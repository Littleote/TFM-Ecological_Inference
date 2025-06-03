wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)

library(sf)
library(spdep)
source("plotting.r")
source("utilities.r")

DATA <- "Louisiana"
df <- read.csv("../data/Louisiana-2020_presidential_election.csv", skip = 1)
geom <- st_read("../data/US-counties/tl_2020_us_county.shp")
df.income <- read.csv("../data/Louisiana-2020_income_adjusted.csv")
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

parse <- function(df, row, stride, offset = 2) {
  as.numeric(gsub(",", "", df[row, stride * 1:data$T - stride + offset]))
}
var.income <- parse(df.income, 12, 8)
var.education.white <- parse(df.education, 31, 4) / parse(df.education, 30, 4)
var.education.black <- parse(df.education, 37, 4) / parse(df.education, 36, 4)
var.education.gap <- var.education.white - var.education.black
var.both <- whitening::whiten(cbind(var.income, var.education.gap))
geom.louisiana <- geom[geom$STATEFP == 22, ]
parish.idx <- order(geom.louisiana$NAME)
adjacency <- 1 == nb2mat(poly2nb(geom.louisiana), style = "B")[parish.idx, parish.idx]
bounds.args <- list(ylim = c(0, 1))
spatial.args <- list(map = geom.louisiana)

source("logit_covariate.r")
source("cluster.r")
source("spatial.r")
source("covariate.r")
source("wrapper.r")

table <- NULL

# LOGIT COVARIATE
run <- curry(
  logit.covariate.model, "logit covariate", data, true, BUGS.config,
  plot.call(covariate.plot, geom.louisiana, NULL, bounds.args),
  burn = 90000, thin = 20, rand.effects = TRUE
)
table <- run(table, "no covariate", X = NULL)
table <- run(table, "income", X = var.income)
table <- run(table, "education", X = var.education.gap)
table <- run(table, "income and education", X = var.both)

# LOGIT COVARIATE WITHOUT RANDOM EFFECT
run <- curry(
  logit.covariate.model, "logit covariate without random effects", data, true, BUGS.config,
  plot.call(covariate.plot, geom.louisiana, NULL, bounds.args),
  burn = 10000, rand.effects = FALSE
)
table <- run(table, "no covariate", X = NULL)
table <- run(table, "income", X = var.income)
table <- run(table, "education", X = var.education.gap)
table <- run(table, "income and education", X = var.both)

# CLUSTER
run <- curry(
  cluster.model, "cluster", data, true, BUGS.config,
  plot.call(cluster.plot, geom.louisiana, NULL, bounds.args),
  burn = 30000
)
table <- run(table, "1 cluster", K = 1)
table <- run(table, "2 clusters", K = 2)
table <- run(table, "3 clusters", K = 3)
table <- run(table, "4 clusters", K = 4)

# COVARIATE
run <- curry(
  covariate.model, "covariate", data, true, BUGS.config,
  plot.call(covariate.plot, geom.louisiana, NULL, bounds.args),
  burn = 70000, thin = 20
)
table <- run(table, "no covariate", X = NULL)
table <- run(table, "income", X = var.income)
table <- run(table, "education", X = var.education.gap)
table <- run(table, "income and education", X = var.both)

# SPATIAL
run <- curry(
  spatial.model, "spatial", data, true, BUGS.config,
  plot.call(spatial.plot, geom.louisiana, spatial.args, bounds.args),
  burn = 30000, BUGS = "WinBUGS", thin = 20
)
table <- run(table, "no adjacency", ADJ = NULL)
table <- run(table, "louisiana map", ADJ = adjacency)

# LPHOM
run <- curry(
  lphom.model, "lphom", data, true, NULL,
  plot.call(NULL, geom.louisiana, NULL, bounds.args)
)
table <- run(table, "nslphom", lphom = lphom::nslphom)

# ECOL RxC
run <- curry(
  ecol.model, "ecolRxC", data, true, NULL,
  plot.call(NULL, geom.louisiana, NULL, bounds.args),
  confidence = 0.95
)
table <- run(table, "ecolRxC", ecol = ecolRxC::ecolRxC)

write.csv(table, file = "louisiana_stats.csv")
save.image(file = "louisiana.RData")

plot.error(table, glue("../plots/{DATA}/error/full"))
mean.sd.tradeoff(table, 2, 1, glue("../plots/{DATA}/error/full"))
bias.variance.tradeoff(table, 2, 1, glue("../plots/{DATA}/error/full"))
