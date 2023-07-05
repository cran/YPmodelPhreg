fn_uhx <- function(b, od, owz, jh = 1E-8, ...) {

  k <- ncol(owz) - 1
  mb <- as.numeric(b) + cbind(0, diag(k+2) * jh)
  u <- fn_ux(mb, od, owz)

  # p <- k + 2 #length(b)
  h <- matrix(NA, k + 2, k + 2)
  for (i in 1:(k + 2)) {
    h[,i] <- (u[,i+1] - u[, 1])/jh
  }

  ep <- 0.0001; cr <- 0

  while ((min(abs(svd(h)$d)) < 1.e-4) & (cr<100)) {
    ep <- 2 * ep
    h <- h + ep * diag(k + 2)
  }

  res <- list()
  res$u <- u[,1]
  res$h <- h

  return(res)

}
