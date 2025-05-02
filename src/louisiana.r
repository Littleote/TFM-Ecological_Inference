wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)

library(sf)
library(spdep)
source("plotting.r")
source("utilities.r")

geom <- st_read("../data/US-counties/tl_2020_us_county.shp")

df <- read.csv("../data/Louisiana-2020_presidential_election.csv", skip = 1)
df.income <- read.csv("../data/Louisiana-2020_income_adjusted.csv")
df.age.gender <- read.csv("../data/Louisiana-2020_age-gender.csv")
df.employment <- read.csv("../data/Louisiana-2020_employment.csv")
df.education <- read.csv("../data/Louisiana-2020_education.csv")

true <- list()
data <- list(
  T = (nrow(df) - 4) / 3,
  R = 2,
  C = 2
)

true$count <- array(unlist(df[
  3 + 3 * (1:data$T),
  paste0(c("WHITE", "BLACK", "OTHER")[1:data$R], ".", rep(1:data$C, each = data$R))
]), c(data$T, data$R, data$C))
true$beta <- aperm(apply(true$count, c(1, 2), norm.1), c(2, 3, 1))

data <- c(data, list(
  Y.1 = apply(true$count, c(1, 2), sum),
  Y.2 = apply(true$count, c(1, 3), sum),
  N.1 = apply(true$count, 1, sum),
  N.2 = apply(true$count, 1, sum)
))
if (all(data$N.1 == data$N.2)) {
  data$N <- data$N.1
}
data$f.1 <- t(apply(data$Y.1, 1, norm.1))
data$f.2 <- t(apply(data$Y.2, 1, norm.1))

rand.covariate <- rnorm(data$T)
parse <- function(df, row, stride, offset = 2) {
  as.numeric(gsub(",", "", df[row, stride * 1:data$T - stride + offset]))
}
var.income <- parse(df.income, 12, 8)
var.age <- parse(df.age.gender, 19, 4)
var.gender <- parse(df.age.gender, 93, 4) / parse(df.age.gender, 92, 4)
var.education.white <- parse(df.education, 31, 4) / parse(df.education, 30, 4)
var.education.black <- parse(df.education, 37, 4) / parse(df.education, 36, 4)
var.employment.white.collar <- (parse(df.employment, 2, 12) + parse(df.employment, 4, 12)) / parse(df.employment, 1, 12)
single.covariate <- var.income
multi.covariate <- rbind(
  var.income,
  var.age,
  var.gender,
  var.education.white,
  var.education.black,
  var.employment.white.collar
)
multi.covariate <- whitening::whiten(t(multi.covariate))
geom.louisiana <- geom[geom$STATEFP == 22, ]
parish.idx <- order(geom.louisiana$NAME)
reorder <- sample(data$T)
adjacency <- 1 == nb2mat(poly2nb(geom.louisiana), style = "B")[parish.idx, parish.idx]
rand.adjacency <- adjacency[reorder, reorder]
bounds.args <- list(ylim = c(0, 1))
spatial.args <- list(map = geom.louisiana)

source("logit_covariate.r")
source("cluster.r")
source("spatial.r")
source("covariate.r")
source("wrapper.r")

curry <- function(fun, fun.name, data, true, config = NULL, plotting = NULL, ...) {
  common.params <- list(...)
  return(function(table, param.name, ...) {
    do.call(
      test.model,
      c(list(
        model = fun,
        name = paste(fun.name, param.name, sep = ", "),
        table = table,
        data = data,
        true = true,
        config = config,
        plotting = plotting,
        ...
      ), common.params)
    )
  })
}

table <- NULL

# LOGIT COVARIATE
run <- curry(
  logit.covariate.model, "logit covariate", data, true, BUGS.config,
  plot.call(covariate.plot, geom.louisiana, NULL, bounds.args),
  burn = 30000
)
table <- run(table, "random covariate", X = rand.covariate)
table <- run(table, "single covariate", X = single.covariate)
table <- run(table, "multiple covariate", X = multi.covariate)

# LOGIT COVARIATE + RANDOMNESS
run <- curry(
  logit.covariate.model, "logit covariate + randomness", data, true, BUGS.config,
  plot.call(covariate.plot, geom.louisiana, NULL, bounds.args),
  burn = 50000, free = TRUE
)
table <- run(table, "random covariate", X = rand.covariate)
table <- run(table, "single covariate", X = single.covariate)
table <- run(table, "multiple covariate", X = multi.covariate)

# CLUSTER
run <- curry(
  cluster.model, "cluster", data, true, BUGS.config,
  plot.call(cluster.plot, geom.louisiana, NULL, bounds.args)
)
table <- run(table, "2 clusters", K = 2)
table <- run(table, "3 clusters", K = 3)
table <- run(table, "4 clusters", K = 4)

# COVARIATE
run <- curry(
  covariate.model, "covariate", data, true, BUGS.config,
  plot.call(covariate.plot, geom.louisiana, NULL, bounds.args),
  burn = 30000
)
table <- run(table, "random covariate", X = rand.covariate)
table <- run(table, "single covariate", X = single.covariate)
table <- run(table, "multiple covariate", X = multi.covariate)

# SPATIAL
run <- curry(
  spatial.model, "spatial", data, true, BUGS.config,
  plot.call(spatial.plot, geom.louisiana, spatial.args, bounds.args),
  burn = 30000, BUGS = "WinBUGS"
)
table <- run(table, "random adjacency", ADJ = rand.adjacency)
table <- run(table, "louisiana adjacency", ADJ = adjacency)

# LPHOM
run <- curry(
  lphom.model, "lphom", data, true, NULL,
  plot.call(NULL, geom.louisiana, NULL, bounds.args)
)
table <- run(table, "lclphom", lphom = lphom::lclphom)
table <- run(table, "nslphom", lphom = lphom::nslphom)
table <- run(table, "rslphom", lphom = lphom::rslphom)
table <- run(table, "tslphom", lphom = lphom::tslphom)

# ECOL RxC
run <- curry(
  ecol.model, "ecolRxC", data, true, NULL,
  plot.call(NULL, geom.louisiana, NULL, bounds.args),
  confidence = 0.95
)
table <- run(table, "ecolRxC", ecol = ecolRxC::ecolRxC)

write.csv(table, file = "louisiana_stats.csv")
save.image(file = "louisiana.RData")
