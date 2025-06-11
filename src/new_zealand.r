wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)

source("utilities.r")
source("plotting.r")
library(ei.Datasets)
library(tidyr)
library(glue)
# LOC <- "Auckland Central"; combine = FALSE
LOC <- "Waiariki"; combine = TRUE

extract <- function(elec, region, threshold = 0.05, combine = TRUE) {
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
  if (combine) {
    aux.table <- elec[[index, 3]][[1]]
    aux.table[
      c(mapply(
        function(x) {
          contains("Special Votes", var = x)
        },
        aux.table["Address"])
      ),
      "City"
    ] <- "Special"
    city <- fill(aux.table, City)$City
    Y.1 <- aggregate(Y.1, list(city), sum)[,-1]
    Y.2 <- aggregate(Y.2, list(city), sum)[,-1]
  }
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
  true$all.beta <- t(apply(true$all, 1, norm.1))
  return(list(data = data, true = true))
}

generate.local <- function(data, global, var, noise, explainability, force.exp, force.noise = 0.2 * force.exp, integer.count = TRUE) {
  sim.data <- runBUGS(
    "JAGS",
    data = c(
      data, list(
        explainability = explainability,
        global.beta = global,
        noise = noise,
        noise.2 = noise,
        force.exp = force.exp,
        force.noise = force.noise,
        X = var
      )
    ),
    param = "beta",
    model = model.file("generate-local"),
    n.chains = 1,
    purge.chains = 0
  )
  
  sim.beta <- sim.data$BUGSoutput$last.values[[1]]$beta
  sim.counts <- sim.beta * array(data$Y.1, dim(sim.beta))
  if (integer.count) {
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
  }
  new.true <- list(
    beta = aperm(apply(sim.counts, c(1, 2), norm.1), c(2, 3, 1)),
    count = sim.counts
  )
  new.true$beta[is.nan(new.true$beta)] <- 1 / data$C
  new.true$global.count <- apply(new.true$count, 2:3, sum)
  new.true$global.beta <- apply(new.true$beta, 2:3, mean)
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
  properties <- list(
    explainability = mean(apply(new.true$beta, 2:3, function(Y) { abs(cor(var, Y)) } )),
    noise = mean(apply(new.true$beta, 2:3, var ) / (new.true$global.beta * (1 - new.true$global.beta)))
  )
  return(list(true = new.true, data = new.data, properties = properties))
}

info <- extract(ei_NZ_2020, LOC, combine = combine)
var.X <- rnorm(info$data$T)
noise.range <- rep(c(0.1, 0.5, 0.9), each = 3)
exp.range <- rep(c(0.1, 0.5, 0.9), times = 3)

source("logit_covariate.r")
source("covariate.r")

table <- NULL

for (i in 1:9) {
  print(i)
  out <- generate.local(
    info$data, info$true$all.beta, var.X, noise = noise.range[i],
    force.exp = ifelse(i %% 3 == 0, 1e5, 1e4),
    force.noise = ifelse(i %% 3 == 0, 1e5, 2e3),
    explainability = exp.range[i], integer = TRUE
  )
  bounds.args <- list(ylim = c(0, 1))
  DATA <- glue("NewZealand/{LOC}-{i}")
  data <- out$data
  true <- out$true

  sub.table <- NULL

  # LOGIT COVARIATE
  print("Logit covariate")
  run <- curry(
    logit.covariate.model, "logit covariate", data, true, BUGS.config,
    plot.call(covariate.plot, NULL, NULL, bounds.args),
    burn = 30000, iter = 2000, thin = 10, rand.effects = TRUE
  )
  sub.table <- run(sub.table, "no covariate", X = NULL)
  sub.table <- run(sub.table, "with covariate", X = var.X)

  # COVARIATE
  print("King's covariate")
  run <- curry(
    covariate.model, "covariate", data, true, BUGS.config,
    plot.call(covariate.plot, NULL, NULL, bounds.args),
    burn = 30000, iter = 2000, thin = 10
  )
  sub.table <- run(sub.table, "no covariate", X = NULL)
  sub.table <- run(sub.table, "with covariate", X = var.X)

  sub.table["Noise"] <- out$properties$noise
  sub.table["Explainability"] <- out$properties$explainability
  table <- rbind(table, sub.table)
}

array(table[["Noise"]], c(4, 3, 3))[1,,]
array(table[["Explainability"]], c(4, 3, 3))[1,,]

write.csv(table, file = glue("out/New Zealand-{LOC}_stats.csv"))
save.image(file = glue("out/New Zealand-{LOC}.RData"))