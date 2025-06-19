library(tidyr)

# Set path once
if(!file.exists("path.RData")){
  WinBUGS.path <- rstudioapi::selectDirectory(
    caption = "Select WinBUGS Directory",
    label = "Select",
    path = rstudioapi::getActiveProject()
  )
  save(file = "path.RData", list = c("WinBUGS.path"))
} else {
  load("path.RData")
}

# Helper functions
norm.1 <- function(v) {
  if(sum(v) < 1e-100) {
    0 * v + 1 / length(v)
  } else {
    v / sum(v)
  }
}

model.file <- function(name) {
  normalizePath(paste0("../models/", name, ".txt"))
}

get.array <- function(arr, shape) {
  aperm(array(arr, rev(shape)), length(shape):1)
}

BUGS.config <- function(..., iter = 2000, burn = 1000, chain = 4, thin = 10, BUGS = "JAGS") {
  return(list(Iterations = iter, Burned = burn, Chains = chain, Thinning = thin, BUGS = BUGS))
}

as.global <- function(X) {
  apply(X, 2:3, sum)
}

# Sort BUGS results
"[.BUGS.param" <- function(x, i) {
  class(x) <- "list"
  x <- x[i]
  class(x) <- "BUGS.param"
  x
}

">.BUGS.param" <- function(a, b) {
  a <- a[[1]]
  b <- b[[1]]
  if (a$name != b$name) {
    return(a$name > b$name)
  }
  diff <- which(a$subindex != b$subindex)[1]
  if (is.na(diff)) {
    return(FALSE)
  }
  return((a$subindex > b$subindex)[diff])
}

"<.BUGS.param" <- function(a, b) b > a

"==.BUGS.param" <- function(a, b) ifelse(a > b || b > a, FALSE, TRUE)

# Utility functions
index <- function(names) {
  if (!is.null(dim(names))) {
    names <- dimnames(names)[[3]]
  }
  params <- mapply(
    function(x, idx) {
      list(
        name = x[[1]],
        subindex = if (length(x) == 1) NULL else as.numeric(unlist(x[2:length(x)])),
        idx = idx
      )
    },
    strsplit(names, "\\[|,|\\]"), seq(length(names)),
    SIMPLIFY = FALSE
  )
  class(params) <- "BUGS.param"
  params <- sort(params)
  function(target, subindex = NULL) {
    unlist(mapply(function(x) {
      if (x$name == target && all(subindex == x$subindex[1:length(subindex)] | subindex == -1)) x$idx
    }, params))
  }
}

runBUGS <- function(BUGS, data, inits = NULL, parameters, model, n.chains, purge.chains = 2, ...) {
  force(model)
  total.chains <- n.chains + purge.chains
  if (BUGS == "JAGS") {
    library(R2jags)
    sim <- R2jags::jags.parallel(
      data, inits, parameters, model,
      n.chains = total.chains, ...,
    )
  } else if (BUGS == "WinBUGS") {
    library(R2WinBUGS)
    BUGSoutput <- R2WinBUGS::bugs(
      data, inits, parameters, model,
      n.chains = total.chains, ...,
      bugs.directory = WinBUGS.path
    )
    sim <- list(
      BUGSoutput = BUGSoutput,
      parameters.to.save = parameters,
      model.file = model,
      n.iter = BUGSoutput$n.iter,
      DIC = BUGSoutput$is.DIC
    )
  } else {
    stop(paste(BUGS, "not allowed, used JAGS or WinBUGS instead."))
  }
  if (purge.chains == 0) {
    return(sim)
  }
  keep <- sim$BUGSoutput$sims.array %>%
    apply(., 2:3, sort) %>%
    apply(., 3, function(x) {
      as.matrix(dist(t(x), method = "man"))
    }) %>%
    apply(., 1, sum) %>%
    array(., c(total.chains, total.chains)) %>%
    apply(., 1, sum) %>%
    order(.) %>%
    .[1:n.chains]
  sim$BUGSoutput$sims.array <- sim$BUGSoutput$sims.array[, keep, ]
  return(sim)
}

hi.low <- function(arr, N, R, C) {
  search <- index(arr)
  bounds <- aperm(array(
    apply(arr[, , search("beta")], 3, quantile, c(0.025, 0.975)),
    c(2, C, R, N),
    list(c("2.5%", "97.5%"), paste("column", 1:C), paste("row", 1:R), NULL)
  ), c(1, 4, 3, 2))
  return(bounds)
}

bias.deviation <- function(out, true, data) {
  R <- data$R
  C <- data$C
  bias <- array(dim = c(R, C))
  deviation <- array(dim = c(R, C))
  for (r in 1:R) {
    for (c in 1:C) {
      error <- out$beta[, r, c] - true$beta[, r, c]
      weight <- data$Y.1[, r]
      bias[r, c] <- sum(weight * error) / sum(weight)
      deviation[r, c] <- sqrt(sum(weight * (error - bias[r, c])^2) / (sum(weight) - 1))
    }
  }

  return(list(
    bias = bias,
    deviation = deviation
  ))
}

get.error <- function(out, true) {
  og.count <- as.global(out$count)
  tg.count <- as.global(true$count)
  list(
    local = list(
      MAE = sum(abs(out$count - true$count)) / sum(true$count),
      MSE = sqrt(sum((out$count - true$count)^2)) / sum(true$count)
    ),
    global = list(
      MAE = sum(abs(og.count - tg.count)) / sum(tg.count),
      MSE = sqrt(sum((og.count - tg.count)^2)) / sum(tg.count)
    )
  )
}

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

test.model <- function(model, name, table = NULL, data, true, ..., config = NULL, plotting = NULL, verbose = FALSE) {
  execution.time <- system.time(out <- model(data, ...))
  out$count <- out$count %||% out$beta * array(data$Y.1, dim(true$beta))

  if (!is.null(plotting)) {
    plotting(name, out, true)
  }

  bias.sd <- bias.deviation(out, true, data)
  error <- get.error(out, true)
  median.bound.size <- apply(out$bounds[2, , , ] - out$bounds[1, , , ], c(2, 3), median)
  in.bounds <- apply(out$bounds[1, , , ] < true$beta & true$beta < out$bounds[2, , , ], c(2, 3), mean)

  elems <- c(
    unlist(bias.sd),
    out$beta[1,,],
    out$count[1,,],
    median.bound.size,
    in.bounds,
    execution.time[3],
    unlist(error$local),
    unlist(error$global)
  )


  if (is.null(table)) {
    table <- data.frame(as.list(elems), row.names = name)
    naming <- paste0("[", 1:data$R, ", ", rep(1:data$C, each = data$R), "]")
    dimnames(table)[[2]] <- c(
      paste(
        rep(
          c(
            "Bias",
            "Deviation",
            "Example Beta",
            "Example Counts",
            "Bounds Size",
            "Accuracy"
          ),
          each = length(naming)
        ),
        naming
      ),
      "Execution Time",
      paste(rep(c("Local", "Global"), each=2), c("MAE", "MSE"))
    )
  } else {
    table[name, 1:length(elems)] <- elems
  }
  table[name, "Model"] <- name
  table[name, "MCMC"] <- ifelse(is.null(out$sim), "False", "True")
  if (!is.null(config)) {
    extra <- config(out, true, data, ...)
    for (var.name in names(extra)) {
      table[name, var.name] <- extra[[var.name]]
    }
  }

  if (verbose) {
    df <- data.frame(list(
      beta.bias = c(beta.bias),
      count.bias = c(count.bias),
      beta.deviation = c(beta.deviation),
      count.deviation = c(count.deviation)
    ), row.names = paste(seq(data$R), rep(seq(data$C), each = data$R)))
    print(df)
  }

  return(table)
}
