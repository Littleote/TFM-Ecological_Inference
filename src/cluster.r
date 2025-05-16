source("utilities.r")

match.clusters <- function(sims.array) {
  arr <- sims.array
  at <- index(arr)
  chain <- dim(arr)[2]
  K <- length(at("omega"))
  # TODO: Add beta into consideration
  if (K == 1) {
    mapping <- array(1, c(K, chain))
  } else {
    mapping <- apply(apply(sims.array[, , at("omega")], c(2, 3), mean), 1, order)
  }
  for (c in seq(chain)) {
    for (cluster in seq(K)) {
      new.cluster <- mapping[cluster, c]
      arr[, c, at("beta", cluster)] <- sims.array[, c, at("beta", new.cluster)]
      arr[, c, at("omega", new.cluster)] <- sims.array[, c, at("omega", new.cluster)]
    }
    arr[, c, at("zeta")] <- factor(
      sims.array[, c, at("zeta")],
      labels = seq(K),
      levels = mapping[, c]
    )
  }
  return(list(array = arr, at = at))
}

stats.mode <- function(x) {
  as.numeric(names(which.max(table(x))))
}

cluster.model <- function(data, K, BUGS = "JAGS", iter = 2000, burn = 1000, chain = 4, thin = 10) {
  data <- c(
    data[c("T", "R", "C", "Y.1", "Y.2", "f.1", "N.1", "N.2")],
    list(K = K)
  )
  file <- if (data$K > 1) model.file("cluster") else model.file("cluster-1_cluster")

  test <- list(
    b = array(1, c(data$K, data$R, data$C)),
    v = array(1, c(data$K)),
    a = array(1, c(data$K, data$R))
  )
  test.sim <- runBUGS(
    BUGS, c(data, test),
    parameters = c("beta", "theta.1", "omega", "zeta"), model = file,
    n.iter = iter * thin + burn, n.burnin = burn, n.thin = thin, n.chains = chain
  )
  test.matching <- match.clusters(test.sim$BUGSoutput$sims.array)
  test.zeta <- apply(test.matching$array[, , test.matching$at("zeta")], 3, stats.mode)
  weight <- 1
  prior <- c(list(
    b = weight * data$C * aperm(array(
      apply(test.matching$array[, , test.matching$at("beta")], 3, median),
      c(data$C, data$R, data$K)
    ), c(3, 2, 1)),
    a = weight * data$R * array(
      apply(
        array(test.matching$array[, , test.matching$at("theta.1")], c(iter, chain, data$R, data$T)),
        3, tapply, rep(test.zeta, each = iter * chain), median
      ),
      c(data$K, data$R)
    )
  ), (if (K > 1) list(v = weight * data$K * apply(test.matching$array[, , test.matching$at("omega")], 3, median)) else NULL))

  sim <- runBUGS(
    BUGS, c(data, prior),
    parameters = c("beta", "theta.1", "theta.2", "omega", "zeta"), model = file,
    n.iter = iter * thin + burn, n.burnin = burn, n.thin = thin, n.chains = chain
  )
  matching <- match.clusters(sim$BUGSoutput$sims.array)
  sim$BUGSoutput$sims.array <- matching$array


  cluster.beta <- aperm(
    array(
      apply(matching$array[, , matching$at("beta")], 3, median),
      c(data$C, data$R, data$K)
    ),
    c(3, 2, 1)
  )
  omega <- if (K > 1) apply(matching$array[, , matching$at("omega")], 3, median) else 1
  zeta <- apply(matching$array[, , matching$at("zeta")], 3, stats.mode)
  beta <- array(NA, c(data$T, data$R, data$C))
  bounds <- array(NA, c(2, data$T, data$R, data$C), list(c("2.5%", "97.5%"), NULL, NULL, NULL))
  cluster.bounds <- hi.low(sim$BUGSoutput$sims.array, data$K, data$R, data$C)
  for (i in seq(data$T)) {
    beta[i, , ] <- cluster.beta[zeta[i], , ]
    bounds[, i, , ] <- cluster.bounds[, zeta[i], , ]
  }
  return(list(
    beta = beta,
    bounds = bounds,
    cluster.beta = cluster.beta,
    omega = omega,
    cluster = zeta,
    sim = sim
  ))
}
