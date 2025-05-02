lphom.model <- function(data, lphom, ...) {
  out <- lphom(data$Y.1, data$Y.2, ...)
  beta <- aperm(out$VTM.prop.units, c(3, 1:2))
  bounds <- aperm(
    array(
      rbind(
        c(out$deterministic.bounds$lower.units),
        c(out$deterministic.bounds$upper.units)
      ),
      c(2, data$R, data$C, data$T)
    ),
    c(1, 4, 2:3)
  )
  return(list(beta = beta, bounds = bounds))
}

ecol.model <- function(data, ecol, ...) {
  out <- ecol(data$Y.1, data$Y.2, ...)
  beta <- aperm(out$VTM.units, c(3, 1:2))
  bounds <- aperm(
    array(
      rbind(
        c(out$VTM.lower.units),
        c(out$VTM.upper.units)
      ),
      c(2, data$R, data$C, data$T)
    ),
    c(1, 4, 2:3)
  )
  return(list(beta = beta, bounds = bounds))
}
