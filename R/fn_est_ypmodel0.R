fn_est_ypmodel0 <- function(oy, od, oz, b0, ...) {

  n <- length(oy)
  cpor <- mean(od)
  u0 <- cpor * oy/oy[n]
  u0 <- u0/(1 - u0)
  dd <- which(od != 0)
  yt <- oy[dd]
  qx <- oy

  nm0 <- 7
  nm <- 50*(cpor < .05) + 25*(cpor >= .05)*(cpor < .1) +
    15*(cpor >= .1)*(cpor < .2) + 10*(cpor >= .2)*(cpor < .5) + 7*(cpor >= .5)

  ftd <- which(od * oz != 0)
  ntd <- length(ftd)
  fcd <- which(od * (1 - oz) != 0)
  ncd <- length(fcd)
  nt3 <- sum(yt <= min(c(oy[ftd[ntd-1]], oy[fcd[ncd-1]])))

  nt1 <- max(c(ftd[ceiling(ntd/2)], fcd[ceiling(ncd/2)]))

  nb <- b0
  tem <- u0

  tem <- u0[dd[nt3]]
  hrt <- exp(b0[1] + b0[2]*qx)
  if (tem < 0.01 - 1 + hrt[dd[nt3]] * exp(-b0[1])) {
    tem = 0.01 -1 + hrt[dd[nt3]] * exp(-b0[1])
  }

  nb[2] <- -log((1 + tem - hrt[dd[nt3]]/hrt[nt1])/tem/hrt[dd[nt3]])

  upb <- log(c(nm0, nm)) - 0.01
  lowb <- -log(c(nm0, nm)) + 0.01
  nb <- apply(rbind(nb, upb), 2, min)
  nb <- apply(rbind(nb, lowb), 2, max)

  # ob <- nb + 0.1
  ob <- nb

  of <- Inf
  ittm <- 100
  cn <- 0
  stoprule <- FALSE
  while (stoprule == FALSE & cn < ittm) {

    cn <- cn + 1
    os_store <- of
    ob_store <- ob

    uh_ob <- fn_uh(ob, od, oz)
    u <- uh_ob$u
    h <- uh_ob$h
    po <- as.numeric(-solve(h) %*% u)
    of <- sum(u^2)

    # testing
    t <- 1
    bc <- ob + t*po

    idb <- sum((abs(bc) > log(c(nm0, nm))))
    tb <- 0

    while ((idb > 0)&(tb < 20)) {
      tb <- tb + 1
      t <- t * 0.3
      bc <-  ob + t*po
      idb <- sum((abs(bc) > log(c(nm0, nm))))
    }

    uh_bc <- fn_uh(bc, od, oz)
    os <- sum(uh_bc$u^2)

    if (os < of) {
      nb <- bc
      of <- os
    } else {
      bc_set <- matrix(NA, 2, 4)
      for (ib in 1:4) {
        bc_set[, ib] <- ob + t *  ib * po/5
      }
      bc_u <- fn_u(bc_set, od, oz)
      bc_os <- colSums(bc_u^2)
      I <- which.min(bc_os)
      nb <- bc_set[,I]
      of <- bc_os[I]
    }

    stoprule <- sum(abs(round(1000*ob) - round(1000*nb))) == 0

    of <- os
    ob <- nb

  }

  if(os_store < os) {
    minb <- ob_store
    os <- os_store
  } else {
    minb <- nb
    os <- os
  }

  res <- list()
  res$nb <- nb
  res$ob <- ob_store
  res$minb <- minb
  res$os <- os

  return(res)

}
