source("utilities.r")

spatial.model <- function(data, ADJ, BUGS, iter = 2000, burn = 1000, chain = 4, thin = 10) {
  stopifnot(data$C == 2, data$R == 2)
  stopifnot(all(data$N1 == data$N2))
  data <- list(
    T = data$T,
    N = data$N,
    f.1 = data$f.1,
    f.2 = data$f.2,
    ADJ.adj = unlist(apply(ADJ, 1, which)),
    ADJ.num = apply(ADJ, 1, sum),
    ADJ.weights = rep(1, sum(ADJ))
  )

  test <- list(
    T = data$T,
    tau.u = 1,
    ADJ.adj = data$ADJ.adj,
    ADJ.num = data$ADJ.num,
    ADJ.weights = data$ADJ.weights
  )
  test.sim <- runBUGS(
    BUGS,
    DIC = FALSE,
    data = c(test), parameters = c("upsilon"), model = model.file("car_normal"),
    n.iter = 1000, n.burnin = 100, n.thin = 1, n.chains = 1
  )
  car.normal <- test.sim$BUGSoutput$sims.array
  factor <- mean(apply(car.normal, 1, sd))

  # Uninformed prior with TMV assuming beta has flat(0, 1) prior
  #   -> logit(beta) has logistic(0, 1) prior
  TMV <- pi^2 / 3
  a <- 1.583179
  r <- qgamma(0.05, shape = a, rate = 1) * TMV
  prior <- list(
    a.u = rep(a, 2),
    r.u = rep(r * factor, 2),
    a.v = rep(a, 2),
    r.v = rep(r, 2)
  )
  sim <- runBUGS(
    BUGS,
    data = c(data, prior), parameters = c("beta", "theta.2", "delta", "nu", "upsilon", "tau.v", "tau.u"), model = model.file("spatial"),
    n.iter = iter * thin + burn, n.burnin = burn, n.thin = thin, n.chains = chain
  )

  ref <- sim$BUGSoutput$median
  arr <- sim$BUGSoutput$sims.array
  at <- index(arr)
  return(list(
    beta = get.array(apply(arr[, , at("beta")], 3, median), dim(ref$beta)),
    bounds = hi.low(sim$BUGSoutput$sims.array, data$T, 2, 2),
    delta = get.array(apply(arr[, , at("delta")], 3, median), dim(ref$delta)),
    nu = get.array(apply(arr[, , at("nu")], 3, median), dim(ref$nu)),
    upsilon = get.array(apply(arr[, , at("upsilon")], 3, median), dim(ref$upsilon)),
    sim = sim
  ))
}
