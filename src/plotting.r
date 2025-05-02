library(glue)

WIDTH <- 900
HEIGHT <- 600
RED <- "#ff7f7f"
ORANGE <- "#ffbf7f"


plot.call <- function(model.plot, map, model.args, bounds.args) {
  function(name, out, true) {
    directory <- paste(c("../plots", strsplit(name, ", ")[[1]]), collapse = "/")
    bounds.args <- bounds.args %||% list()
    bounds.args$directory <- directory
    bounds.args$out <- out
    bounds.args$true <- true
    do.call(bounds.plot, bounds.args)

    if (!is.null(model.plot)) {
      model.args <- model.args %||% list()
      model.args$directory <- directory
      model.args$out <- out
      do.call(model.plot, model.args)
    }

    beta.map.plot(out, true, map, directory)
  }
}

bounds.plot <- function(out, true, ylim = NULL, directory = NULL, ...) {
  # Override save functions if not saving
  if (is.null(directory)) {
    png <- function(...) {}
    dev.off <- function(...) {}
  } else {
    if (!dir.exists(directory)) dir.create(directory, recursive = TRUE)
  }

  T <- dim(true$beta)[1]
  R <- dim(true$beta)[2]
  C <- dim(true$beta)[3]
  calc.ylim <- is.null(ylim)
  png(file = glue("{directory}/beta_bounds.png"), width = WIDTH, height = HEIGHT)
  par(mfrow = c(R, C), mar = c(2, 2, 2, 0) + 0.1)
  for (r in 1:R) {
    for (c in 1:C) {
      if (calc.ylim) {
        ylim <- range(c(true$beta, out$bounds))
      }
      match <- out$bounds[1, , r, c] < true$beta[, r, c] & true$beta[, r, c] < out$bounds[2, , r, c]
      plot(1, xlim = c(1, T), ylim = ylim, main = bquote(beta[.(r) ~ .(c)]), type = "n", ...)
      for (t in 1:T) {
        lines(c(t, t), out$bounds[, t, r, c], col = "blue", lwd = 3, lend = 3)
      }
      points(out$beta[, r, c], col = "black", pch = 18)
      points(true$beta[, r, c], col = ifelse(match, "green", "red"), pch = 20)
    }
  }
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
  dev.off()
}

bias.variance.tradeoff <- function(table, R, C) {
  groups <- factor(mapply(
    function(x) {
      strsplit(x, ", ")[[1]][1]
    },
    dimnames(table)[[1]]
  ))
  par(mar = c(4, 4, 2, 0) + 0.1)
  for (r in 1:R) {
    for (c in 1:C) {
      for (kind in c("Individual", "Region")) {
        bias.2 <- table[, glue("{kind} bias [{r}, {c}]")]^2
        variance <- table[, glue("{kind} deviation [{r}, {c}]")]^2
        plot.tradeoff(
          bias.2, variance, groups,
          mirror = F,
          xlab = bquote(Bias^2 * .("[") * beta[.(r) ~ .(c)] * .("]")),
          ylab = bquote(Var * .("[") * beta[.(r) ~ .(c)] * .("]")),
          main = bquote(.(kind) * .(" bias-variance for ") * beta[.(r) ~ .(c)])
        )
      }
    }
  }
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
}

mean.sd.tradeoff <- function(table, R, C) {
  groups <- factor(mapply(
    function(x) {
      strsplit(x, ", ")[[1]][1]
    },
    dimnames(table)[[1]]
  ))
  par(mar = c(4, 4, 2, 0) + 0.1)
  for (r in 1:R) {
    for (c in 1:C) {
      for (kind in c("Individual", "Region")) {
        mean <- table[, glue("{kind} bias [{r}, {c}]")]
        sd <- table[, glue("{kind} deviation [{r}, {c}]")]
        plot.tradeoff(
          mean, sd, groups,
          mirror = T,
          xlab = bquote(Bias * .("[") * beta[.(r) ~ .(c)] * .("]")),
          ylab = bquote(sigma * .("[") * beta[.(r) ~ .(c)] * .("]")),
          main = bquote(.(kind) * .(" bias-deviation for ") * beta[.(r) ~ .(c)])
        )
      }
    }
  }
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
}

plot.tradeoff <- function(x, y, groups, mirror, chain = TRUE, ...) {
  # Plot margins
  margin <- 1.1
  xlim <- margin * max(abs(x))
  ylim <- margin * max(y)
  plot(
    1,
    xlim = c(ifelse(mirror, -xlim, 0), xlim), ylim = c(0, ylim),
    type = "n", xaxs = "i", yaxs = "i", ...
  )

  # Bias regions
  polygon(c(0, xlim, xlim, 0), c(0, xlim * 1 / qnorm(0.90), 0, 0), col = ORANGE, border = NA)
  polygon(c(0, xlim, xlim, 0), c(0, xlim * 1 / qnorm(0.95), 0, 0), col = RED, border = NA)
  if (mirror) {
    polygon(c(0, -xlim, -xlim, 0), c(0, xlim * 1 / qnorm(0.90), 0, 0), col = ORANGE, border = NA)
    polygon(c(0, -xlim, -xlim, 0), c(0, xlim * 1 / qnorm(0.95), 0, 0), col = RED, border = NA)
    abline(v = 0, lty = 4)
  }


  components <- tapply(1:length(groups), groups, function(x) {
    x
  })
  subindex <- unlist(components)
  subindex[subindex] <- unlist(tapply(
    groups, groups, function(x) {
      1:length(x)
    }
  ))
  palette <- rainbow(length(levels(groups)))

  points(x, y, pch = 14 + subindex, col = palette[groups])
  if (chain) {
    for (comp in 1:length(components)) {
      n <- length(components[[comp]]) - 1
      if (n >= 1) {
        idx <- components[[comp]]
        arrows(
          x[idx[1:n]], y[idx[1:n]],
          x[idx[1 + 1:n]], y[idx[1 + 1:n]],
          length = 0.1, col = palette[comp]
        )
      }
    }
  }
  legend("topright", legend = levels(groups), fill = palette, border = NA)
}

plot.error <- function(table) {
  groups <- factor(mapply(
    function(x) {
      strsplit(x, ", ")[[1]][1]
    },
    row.names(table)
  ))
  config <- mapply(
    function(x) {
      strsplit(x, ", ")[[1]][2]
    },
    row.names(table)
  )

  palette <- rainbow(length(levels(groups)))
  par(mar = c(8, 2, 2, 0) + 0.1)
  for (err in c("EI", "EPW", "EQ")) {
    barplot(table[, err], names.arg = config, col = palette[groups], las = 3, main = err)
    legend("topright", legend = levels(groups), fill = palette, border = NA)
  }
  par(mar = c(2, 2, 2, 0) + 0.1)
}

BUGS.trace <- function(arr, var, search = NULL, subindex = NULL, layout = NULL, ylim = NULL, ...) {
  # Extract dimensions and indexes
  iter <- dim(arr)[1]
  chain <- dim(arr)[2]
  var.names <- dimnames(arr)[[3]]
  search <- search %||% index(var.names)
  idx <- search(var, subindex)
  stopifnot(length(idx) > 0)

  # Porces NULL
  calc.ylim <- is.null(ylim)
  layout <- layout %||% rep(ceiling(sqrt(length(idx))), 2)

  # Plot traces
  par(mfrow = layout, mar = c(2, 2, 2, 0) + 0.1)
  for (i in idx) {
    if (calc.ylim) {
      ylim <- range(arr[, , i])
    }
    plot(1, xlim = c(0, iter), ylim = ylim, type = "n", main = var.names[i], ...)
    for (c in seq(chain)) {
      lines(arr[, c, i], col = rainbow(chain)[c])
    }
  }
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
}

covariate.plot <- function(out, t.seq = NULL, directory = NULL) {
  # Override save functions if not saving
  if (is.null(directory)) {
    png <- function(...) {}
    dev.off <- function(...) {}
  } else {
    if (!dir.exists(directory)) dir.create(directory, recursive = TRUE)
  }

  # Extract variables
  arr <- out$sim$BUGSout$sims.array
  var.names <- dimnames(arr)[[3]]
  search <- index(var.names)
  T <- dim(out$beta)[1]
  R <- dim(out$beta)[2]
  C <- dim(out$beta)[3]
  P <- dim(out$gamma)[1]
  t.seq <- t.seq %||% seq(T)
  p.seq <- seq(P)

  # Generate plots
  for (t in t.seq) {
    png(file = glue("{directory}/beta_{t}.png"), width = WIDTH, height = HEIGHT)
    BUGS.trace(arr, var = "beta", search = search, subindex = c(t))
    dev.off()
    png(file = glue("{directory}/theta(2)_{t}.png"), width = WIDTH, height = HEIGHT)
    BUGS.trace(arr, var = "theta.2", search = search, subindex = c(t), layout = c(C, 1))
    dev.off()
  }
  png(file = glue("{directory}/delta.png"), width = WIDTH, height = HEIGHT)
  BUGS.trace(arr, var = "delta", search = search)
  dev.off()

  for (p in p.seq) {
    png(file = glue("{directory}/gamma_{p}.png"), width = WIDTH, height = HEIGHT)
    BUGS.trace(arr, var = "gamma", search = search, subindex = c(p))
    dev.off()
  }
  if (!is.null(search("lambda"))) {
    png(file = glue("{directory}/lambda.png"), width = WIDTH, height = HEIGHT)
    BUGS.trace(arr, var = "lambda", search = search)
    dev.off()
  }
}

cluster.plot <- function(out, t.seq = NULL, directory = NULL) {
  # Override save functions if not saving
  if (is.null(directory)) {
    png <- function(...) {}
    dev.off <- function(...) {}
  } else {
    if (!dir.exists(directory)) dir.create(directory, recursive = TRUE)
  }

  # Extract variables
  arr <- out$sim$BUGSout$sims.array
  var.names <- dimnames(arr)[[3]]
  search <- index(var.names)
  T <- dim(out$beta)[1]
  R <- dim(out$beta)[2]
  C <- dim(out$beta)[3]
  K <- dim(out$cluster.beta)[1]
  t.seq <- t.seq %||% seq(T)
  k.seq <- seq(K)

  # Generate plots
  for (k in k.seq) {
    png(file = glue("{directory}/beta_{k}.png"), width = WIDTH, height = HEIGHT)
    BUGS.trace(arr, var = "beta", search = search, subindex = c(k))
    dev.off()
  }
  png(file = glue("{directory}/omega.png"), width = WIDTH, height = HEIGHT)
  BUGS.trace(arr, var = "omega", search = search, layout = c(K, 1))
  dev.off()

  for (t in t.seq) {
    png(file = glue("{directory}/zeta_{t}.png"), width = WIDTH, height = HEIGHT)
    BUGS.trace(arr, var = "zeta", search = search, subindex = c(t))
    dev.off()
    png(file = glue("{directory}/theta(1)_{t}.png"), width = WIDTH, height = HEIGHT)
    BUGS.trace(arr, var = "theta.1", search = search, subindex = c(t), layout = c(R, 1))
    dev.off()
    png(file = glue("{directory}/theta(2)_{t}.png"), width = WIDTH, height = HEIGHT)
    BUGS.trace(arr, var = "theta.2", search = search, subindex = c(t), layout = c(C, 1))
    dev.off()
  }
}

spatial.plot <- function(out, t.seq = NULL, directory = NULL, map = NULL) {
  # Override save functions if not saving
  if (is.null(directory)) {
    png <- function(...) {}
    dev.off <- function(...) {}
  } else {
    if (!dir.exists(directory)) dir.create(directory, recursive = TRUE)
  }

  # Extract variables
  arr <- out$sim$BUGSout$sims.array
  var.names <- dimnames(arr)[[3]]
  search <- index(var.names)
  T <- dim(out$beta)[1]
  t.seq <- t.seq %||% seq(T)

  # Generate plots
  for (t in t.seq) {
    png(file = glue("{directory}/beta_{t}.png"), width = WIDTH, height = HEIGHT)
    BUGS.trace(arr, var = "beta", search = search, subindex = c(t))
    dev.off()
    png(file = glue("{directory}/theta(2)_{t}.png"), width = WIDTH, height = HEIGHT)
    BUGS.trace(arr, var = "theta.2", search = search, subindex = c(t), layout = c(2, 1))
    dev.off()
    png(file = glue("{directory}/nu_{t}.png"), width = WIDTH, height = HEIGHT)
    BUGS.trace(arr, var = "nu", search = search, subindex = c(-1, t), layout = c(2, 1))
    dev.off()
    png(file = glue("{directory}/upsilon_{t}.png"), width = WIDTH, height = HEIGHT)
    BUGS.trace(arr, var = "upsilon", search = search, subindex = c(-1, t), layout = c(1, 1))
    dev.off()
  }
  png(file = glue("{directory}/delta.png"), width = WIDTH, height = HEIGHT)
  BUGS.trace(arr, var = "delta", search = search, layout = c(2, 1))
  dev.off()
  png(file = glue("{directory}/tau(nu).png"), width = WIDTH, height = HEIGHT)
  BUGS.trace(arr, var = "tau.v", search = search, layout = c(2, 1))
  dev.off()
  png(file = glue("{directory}/tau(upsilon).png"), width = WIDTH, height = HEIGHT)
  BUGS.trace(arr, var = "tau.u", search = search, layout = c(1, 1))
  dev.off()

  if (!is.null(map)) {
    png(file = glue("{directory}/upsilon_map.png"), width = WIDTH, height = HEIGHT)
    upsilon <- apply(arr[, , search("upsilon")], 3, median)
    map.plot(map, upsilon, "upsilon")
    dev.off()
    for (i in 1:2) {
      png(file = glue("{directory}/nu[{i}]_map.png"), width = WIDTH, height = HEIGHT)
      nu <- apply(arr[, , search("nu", i)], 3, median)
      map.plot(map, nu, glue("nu[{i}]"))
      dev.off()
    }
  }
}

map.plot <- function(map, var, var.name, breaks = NULL, pal = NULL, centered = TRUE) {
  if (is.null(breaks)) {
    lim <- max(abs(range(var)))
    step <- floor(log(lim / 5, 10))
    top <- round(lim, -step)
    breaks <- seq(ifelse(centered, -top, 0), top, by = 10^step)
  }
  pal <- pal %||% ifelse(centered, colorspace::diverge_hsv, colorspace::sequential_hcl)
  map$var <- var
  plot(map["var"], main = var.name, pal = pal, breaks = breaks)
}

beta.map.plot <- function(out, true, map, directory = NULL) {
  # Override save functions if not saving
  if (is.null(directory)) {
    png <- function(...) {}
    dev.off <- function(...) {}
  } else {
    if (!dir.exists(directory)) dir.create(directory, recursive = TRUE)
  }

  # Generate plots
  R <- dim(out$beta)[2]
  C <- dim(out$beta)[2]
  for (r in 1:R) {
    for (c in 1:C) {
      png(file = glue("{directory}/beta[{r},{c}]_map.png"), width = WIDTH, height = HEIGHT)
      map.plot(map, out$beta[, r, c], glue("beta[{r},{c}]"), seq(0, 1, 0.01))
      dev.off()
      png(file = glue("{directory}/beta[{r},{c}]_error_map.png"), width = WIDTH, height = HEIGHT)
      map.plot(map, abs(out$beta[, r, c] - true$beta[, r, c]), glue("beta[{r},{c}] error"), centered = FALSE)
      dev.off()
    }
  }
}
