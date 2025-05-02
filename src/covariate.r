source("utilities.r")

error <- function(param, data) {
  # Reconstruct matrices
  delta <- matrix(
    c(param[data$delta.idx], rep(0, data$R)),
    ncol = data$C, nrow = data$R
  )
  gamma <- array(
    c(param[data$gamma.idx], rep(0, data$R * data$P)),
    c(data$P, data$R, data$C)
  )

  # Sum of squared errors
  score <- 0
  for (i in seq(data$T)) {
    beta <- t(apply(exp(delta + apply(gamma * data$X[i, ], c(2, 3), sum)), 1, norm.1))
    pred <- data$f.1[i, ] %*% beta
    score <- score + sum((data$f.2[i, ] - pred)^2)
  }
  return(score)
}

numerical.model <- function(data, X, verbose = FALSE, maxit = 5000) {
  C.minus.1 <- data$C - 1

  # Standardize to prevent numerical issues
  data$P <- length(X) / data$T
  data$X <- array((X - mean(X)) / sd(X), c(data$T, data$P))

  # Optimization parameters
  param <- c(
    matrix(0, ncol = C.minus.1, nrow = data$R), # delta
    rep(matrix(0, ncol = C.minus.1, nrow = data$R), each = data$P) # gamma
  )
  data$delta.idx <- seq(C.minus.1 * data$R)
  data$gamma.idx <- seq(C.minus.1 * data$R * data$P) + C.minus.1 * data$R

  # Run optimization
  best <- optim(param, error, data = data, control = list(
    trace = ifelse(verbose, 1, 0), maxit = maxit
  ))
  if (best$convergence != 0) {
    print(paste(best$convergence, best$message))
    stop()
  }

  # Reconstruct matrices
  delta <- matrix(
    c(best$par[data$delta.idx], rep(0, data$R)),
    ncol = data$C, nrow = data$R
  )
  gamma <- array(
    c(best$par[data$gamma.idx], rep(0, data$R * data$P)),
    c(data$P, data$R, data$C)
  )
  beta <- array(NA, c(data$T, data$R, data$C))
  for (i in seq(data$T)) {
    beta[i, , ] <- t(apply(exp(delta + apply(gamma * data$X[i, ], c(2, 3), sum)), 1, norm.1))
  }
  return(list(beta = beta, delta = delta, gamma = gamma))
}

covariate.model <- function(data, X, BUGS = "JAGS", approx = NULL, iter = 2000, burn = 1000, chain = 4, thin = 10, ...) {
  stopifnot(all(data$N1 == data$N2))

  if (is.null(approx)) {
    approx <- numerical.model(data, X, ...)
  }
  C.minus.1 <- data$C - 1
  # Standardize to prevent numerical issues
  extra <- list()
  extra$P <- length(X) / data$T
  extra$X <- array((X - mean(X)) / sd(X), c(data$T, extra$P))
  data <- c(
    data[c("T", "R", "C", "f.1", "Y.2", "N")],
    extra
  )
  prior <- list(
    l = 1
  )

  # Initial randomness
  S <- 100
  init.fn <- function() {
    list(
      beta = approx$beta,
      delta.aux = approx$delta[, 1:C.minus.1, drop = FALSE],
      gamma.aux = approx$gamma[, , 1:C.minus.1, drop = FALSE],
      lambda = approx$lambda %||% 1
    )
  }
  sim <- runBUGS(
    BUGS, c(data, prior),
    parameters = c("beta", "theta.2", "gamma", "delta", "lambda"),
    model = model.file("covariate"), inits = init.fn,
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
