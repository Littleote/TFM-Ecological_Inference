source("utilities.r")
library(glue)

logit.covariate.model <- function(data, X, BUGS = "JAGS", rand.effects = TRUE, iter = 2000, burn = 1000, chain = 4, thin = 10) {
  extra <- list()
  no.args <- TRUE
  if (!is.null(X)) {
    no.args <- FALSE
    extra$P <- length(X) / data$T
    extra$X <- array(X, c(data$T, extra$P))
  }
  data <- c(
    data[c("T", "R", "C", "Y.1", "Y.2", "N.1", "N.2")],
    extra
  )
  prior <- c(
    list(d = array(1, c(data$R, data$C))),
    (if (rand.effects) list(l = 1 + ifelse(no.args, 0, data$P)) else NULL),
    (if (!no.args) list(g = array(1, c(data$P, data$R, data$C))) else NULL)
  )

  sim <- runBUGS(
    BUGS, c(data, prior),
    parameters = c(
      "beta", "theta.2", "delta",
      (if (rand.effects) c("lambda") else NULL),
      (if (!no.args) c("gamma") else NULL)
    ),
    model = model.file(glue("logit_covariate{ifelse(rand.effects, '', '-no_effects')}{ifelse(no.args, '-no_args', '')}")),
    n.iter = iter * thin + burn, n.burnin = burn, n.thin = thin, n.chains = chain
  )

  ref <- sim$BUGSoutput$median
  arr <- sim$BUGSoutput$sims.array
  at <- index(arr)
  return(c(
    list(
      beta = get.array(apply(arr[, , at("beta")], 3, median), dim(ref$beta)),
      bounds = hi.low(sim$BUGSoutput$sims.array, data$T, data$R, data$C),
      delta = get.array(apply(arr[, , at("delta")], 3, median), dim(ref$delta)),
      sim = sim
    ),
    (if (!no.args) list(gamma = get.array(apply(arr[, , at("gamma")], 3, median), dim(ref$gamma))) else NULL)
  ))
}
