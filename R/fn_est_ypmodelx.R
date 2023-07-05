fn_est_ypmodelx <- function(oy, od, owz, b0, ...) {

  n <- length(oy)
  cpor <- mean(od)
  dd <- which(od != 0)

  kc <- ncol(owz)
  k <- kc - 1
  m <- 2*(k+2)
  tem <- diag(k + 2)
  tec <- tem[,2]

  nm0 <- 7
  nm <- 50*(cpor < .05) + 25*(cpor >= .05)*(cpor < .1) +
    15*(cpor >= .1)*(cpor < .2) + 10*(cpor >= .2)*(cpor < .5) + 7*(cpor >= .5)

  # ob <- b0 + 0.1
  ob <- b0

  of <- Inf
  ittm <- 100
  cn <- 0
  stoprule <- FALSE
  while (stoprule == FALSE & cn < ittm) {

    cn <- cn + 1
    os_store <- of
    ob_store <- ob

    uh_ob <- fn_uhx(ob, od, owz)
    u <- uh_ob$u
    h <- uh_ob$h
    po <- as.numeric(-solve(h) %*% u)
    of <- sum(u^2)

    # testing
    t <- 1
    bc <- ob + t*po

    idb <- sum((abs(bc) > log(nm0 * rep(1, k+2) + (nm - nm0) * tec)))
    tb <- 0

    while ((idb > 0)&(tb < 20)) {
      tb <- tb + 1
      t <- t * 0.3
      bc <-  ob + t*po
      idb <- sum((abs(bc) > log(nm0 * rep(1, k+2) + (nm - nm0) * tec)))
    }

    uh_bc <- fn_uhx(bc, od, owz)
    os <- sum(uh_bc$u^2)

    if (os < of) {
      nb <- bc
      of <- os
    } else {
      bc_set <- matrix(NA, k+2, m)
      for (ib in 1:m) {
        bc_set[, ib] <- ob + t *  ib * po/(m + 1)
      }
      bc_u <- fn_ux(bc_set, od, owz)
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
