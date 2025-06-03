wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)

source("utilities.r")
source("plotting.r")
library(ei.Datasets)
DATA <- "NewZealand"

extract <- function(elec, region, threshold = 0.05) {
  index <- which(elec[,2] == region)
  count <- as.matrix(elec[[index, 5]][[1]][,-1])
  R <- nrow(count)
  C <- ncol(count)
  total <- sum(count)
  idx.1 <- which(apply(count, 1, sum) / total > threshold)
  other.1 <- (1:R)[-idx.1]
  idx.2 <- which(apply(count, 2, sum) / total > threshold)
  other.2 <- (1:C)[-idx.2]
  data <- list()
  true <- list()
  Y.1 <- as.data.frame(elec[[index, 3]][[1]][,-(1:2)])
  Y.2 <- as.data.frame(elec[[index, 4]][[1]][,-(1:2)])
  data$R <- length(idx.1) + 1
  data$C <- length(idx.2) + 1
  data$T <- nrow(Y.1)
  data$Y.1 <- Y.1[,idx.1]
  data$Y.1[,"Other"] <- apply(Y.1[,other.1], 1, sum)
  data$Y.1 <- as.matrix(data$Y.1)
  data$Y.2 <- Y.2[,idx.2]
  data$Y.2[,"Other"] <- apply(Y.2[,other.2], 1, sum)
  data$Y.2 <- as.matrix(data$Y.2)
  true$all <- matrix(NA, data$R, data$C)
  true$all[1:length(idx.1), 1:length(idx.2)] <- count[idx.1, idx.2]
  true$all[1:length(idx.1),data$C] <- apply(count[idx.1, other.2], 1, sum)
  true$all[data$R,1:length(idx.2)] <- apply(count[other.1, idx.2], 2, sum)
  true$all[data$R, data$C] <- sum(count[other.1, other.2])
  dimnames(true$all) <- list(names(data$Y.1), names(data$Y.2))
  true$all.beta <- t(apply(info$true$all, 1, norm.1))
  return(list(data = data, true = true))
}

generate.local <- function(data, global, var, noise, force, explainability) {
  sim.data <- runBUGS(
    "JAGS",
    data = c(
      data, list(
        explainability = explainability,
        global.beta = global,
        noise = noise,
        force = force,
        X = var
      )
    ),
    param = "beta",
    model = model.file("generate-local"),
    n.chains = 1,
    purge.chains = 0
  )
  
  sim.beta <- sim.data$BUGSoutput$last.values[[1]]$beta
  sim.counts <- sim.beta
  for (t in 1:data$T) {
    for (r in 1:data$R) {
      y <- data$Y.1[t, r]
      sim.counts[t,r,] <- y * sim.beta[t, r,]
      for(iter in 1:10) {
        exact <- sum(round(sim.counts[t,r,]))
        if (exact == y) break
        sim.counts[t,r,] <- sim.counts[t,r,] - sign(exact - y) * 2 ^ -iter
      }
      sim.counts[t,r,] <- round(sim.counts[t,r,])
    }
  }
  new.true <- list(
    beta = aperm(apply(sim.counts, c(1, 2), norm.1), c(2, 3, 1)),
    count = sim.counts
  )
  new.true$beta[is.nan(new.true$beta)] <- 1 / data$C
  new.data <- list(
    T = data$T,
    R = data$R,
    C = data$C
  )
  
  new.data <- c(new.data, list(
    Y.1 = apply(new.true$count, c(1, 2), sum),
    Y.2 = apply(new.true$count, c(1, 3), sum),
    N.1 = apply(new.true$count, 1, sum),
    N.2 = apply(new.true$count, 1, sum),
    N = apply(new.true$count, 1, sum)
  ))
  new.data$f.1 <- t(apply(new.data$Y.1, 1, norm.1))
  new.data$f.2 <- t(apply(new.data$Y.2, 1, norm.1))
  return(list(true = new.true, data = new.data))
}

info <- extract(ei_NZ_2020, "Auckland Central")

X <- rnorm(info$data$T)
out <- generate.local(info$data, info$true$all.beta, X, noise = 0.2, force = 1e4, explainability = 0.9)
bounds.args <- list(ylim = c(0, 1))
data <- out$data
true <- out$true

source("logit_covariate.r")
source("covariate.r")

table <- NULL

# LOGIT COVARIATE
run <- curry(
  logit.covariate.model, "logit covariate", data, true, BUGS.config,
  plot.call(covariate.plot, NULL, NULL, bounds.args),
  burn = 30000, thin = 10, rand.effects = TRUE
)
table <- run(table, "no covariate", X = NULL)
table <- run(table, "with covariate", X = X)

# COVARIATE
run <- curry(
  covariate.model, "covariate", data, true, BUGS.config,
  plot.call(covariate.plot, NULL, NULL, bounds.args),
  burn = 30000, thin = 10
)
table <- run(table, "no covariate", X = NULL)
table <- run(table, "income", X = X)

write.csv(table, file = "new_zealand_stats.csv")
save.image(file = "new_zealand.RData")

plot.error(table, glue("../plots/{DATA}/error/full"))
mean.sd.tradeoff(table, 2, 1, glue("../plots/{DATA}/error/full"))
bias.variance.tradeoff(table, 2, 1, glue("../plots/{DATA}/error/full"))
