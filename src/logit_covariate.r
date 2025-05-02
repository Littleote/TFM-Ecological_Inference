source("utilities.r")

logit.covariate.model <- function(data, X, BUGS = "JAGS", free = FALSE, iter = 2000, burn = 1000, chain = 4, thin = 10) {
  extra <- list()
  extra$P <- length(X) / data$T
  extra$X <- array(X, c(data$T, extra$P))
  data <- c(
    data[c("T", "R", "C", "Y.1", "Y.2", "N.1", "N.2")],
    extra
  )
  prior <- c(
    list(
      d = array(1, c(data$R, data$C)),
      g = array(1, c(data$P, data$R, data$C))
    ),
    (if (free) list(l = 1 + data$P) else NULL)
  )
  sim <- runBUGS(
    BUGS, c(data, prior),
    parameters = c(
      "beta", "theta.2", "delta", "gamma", (if (free) c("lambda") else NULL)
    ),
    model = model.file(ifelse(free, "logit_covariate+free", "logit_covariate")),
    n.iter = iter * thin + burn, n.burnin = burn, n.thin = thin, n.chains = chain
  )

  ref <- sim$BUGSoutput$median
  arr <- sim$BUGSoutput$sims.array
  at <- index(arr)
  return(list(
    beta = get.array(apply(arr[, , at("beta")], 3, median), dim(ref$beta)),
    bounds = hi.low(sim$BUGSoutput$sims.array, data$T, data$R, data$C),
    delta = get.array(apply(arr[, , at("delta")], 3, median), dim(ref$delta)),
    gamma = get.array(apply(arr[, , at("gamma")], 3, median), dim(ref$gamma)),
    sim = sim
  ))
}
