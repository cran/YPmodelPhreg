fn_uh <- function(b, od, oz, jh = 1E-8, ...) {

  mb <- as.numeric(b) + cbind(0, diag(2) %x% t(c(jh, -jh)))
  u <- fn_u(mb, od, oz)
  h <- rbind((u[,2] - u[,3])/jh/2, (u[,4] - u[,5])/jh/2)
  h[1,2] <- h[1,2]/2 + h[2,1]/2
  h[2,1] <- h[1,2]

  ep <- 0.0001; cr <- 0
  kc <- nrow(h)

  while ((min(abs(svd(h)$d)) < 1.e-4) & (cr<100)) {
    ep <- 2 * ep
    h <- h + ep * diag(kc)
  }

  res <- list()
  res$u <- u[,1]
  res$h <- h

  return(res)

}
