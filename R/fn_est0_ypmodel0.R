fn_est0_ypmodel0 <- function(oy, od, oz, b0, ...) {

  n <- length(oy)
  cpor <- mean(od)
  u0 <- cpor * oy/oy[n]
  u0 <- u0/(1 - u0)
  qx <- oy
  dd <- which(od != 0)
  tk <- length(dd)
  yt <- oy[dd]

  nb <- b0
  ob <- nb + 0.1
  jh <- 1E-8
  m <- 3

  nm <- 50*(cpor < .05) + 25*(cpor >= .05)*(cpor < .1) +
    15*(cpor >= .1)*(cpor < .2) + 10*(cpor >= .2)*(cpor < .5) + 7*(cpor >= .5)

  ittm <- 100
  cn <- 0
  stoprule <- FALSE
  while (stoprule == FALSE & cn < ittm) {

    b <- as.numeric(nb) + cbind(0, diag(jh, 2))

    S0 <- matrix(1, m, n)
    S11 <- matrix(1, m, n)
    S12 <- matrix(1, m, n)

    for (ii in 1:m) {
      for (jj in 1:tk) {
        tb <- b[, ii]
        ki <- dd[jj]
        expbz <- exp(tb[1] * oz + tb[2] * qx[ki] * oz)
        S0[ii, ki] <- sum((oy >= oy[ki]) * expbz)
        S11[ii, ki] <- sum((oy >= oy[ki]) * expbz * oz)
        S12[ii, ki] <- sum((oy >= oy[ki]) * qx[ki] * expbz * oz)
      }
    }

    u1 <- u2 <- rep(1, m)
    for (numi in 1:m) {
      u1[numi] <- sum(od * (-oz + (S11[numi,]/S0[numi,])))
      u2[numi] <- sum(od * (-oz * qx + (S12[numi,]/S0[numi,])))
    }

    u <- rbind(u1, u2)
    us_ob <- sum(u[,1]^2)

    h <- diag(2)
    for (i in 1:2) {
      h[,i] <- (u[,i+1] - u[,1])/jh
    }

    ep <- 0.0001
    cr <- 0

    while ((min(abs(svd(h)$d)) < 1.e-4) & (cr<100)) {
      ep <- 2 * ep
      h <- h + ep * diag(2) ## kc?
    }

    po <- -as.numeric(solve(h) %*% u[, 1])

    cn <- cn + 1
    t <- 1
    bc <- nb + t*po
    idb <- sum((abs(bc) > log(nm)))
    tb <- 0

    while ((idb > 0)&(tb < 20)) {
      tb <- tb + 1
      t <- t * 0.618
      bc <-  nb + t*po
      idb <- sum((abs(bc) > log(nm)))
    }

    nb <- bc
    stoprule <- sum(abs(round(1000*ob) - round(1000*nb))) == 0

    ob_store <- ob
    ob <- nb

    us_nb <- sum(u[,1]^2)

    # os <- min(sum(u^2))
    # of <- os
  }

  if(us_nb > us_ob) {
    minb <- ob_store
    os <- us_ob
  } else {
    minb <- nb
    os <- us_nb
  }

  res <- list()
  res$nb <- nb
  res$ob <- ob_store
  res$minb <- minb
  res$os <- os

  return(res)

}
