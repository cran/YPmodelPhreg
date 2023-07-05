fn_ux <- function(mb, od, owz, ...) {

  n <- nrow(owz)
  k <- ncol(owz) - 1
  oz <- owz[, 1]
  ox <- owz[, -1]

  if (is.null(nrow(mb))) {

    gama1 <- as.numeric(exp(-oz * mb[1] - ox %*% mb[3:(k+2)]))
    gama2 <- as.numeric(exp(-oz * mb[2] - ox %*% mb[3:(k+2)]))

    sm <- n:1
    dlam <- od/sm

    lam2 <- cumsum(dlam * gama2)
    p <- exp(-lam2)
    pl <- c(1, p[1:(n-1)])
    ru <- cumsum(pl * dlam * gama1)/p

    deni <- gama1 + gama2 * ru

    lik <- sum(od * log(deni) + log(deni/gama1)/gama2)

    u <- rep(NA, 5)
    u[1] <- sum((oz*od) * (-gama1/deni) + oz * ru/deni)
    u[2] <- sum((oz*od) * (-gama2*ru/deni) - oz * ru/deni + oz * log(deni/gama1)/gama2)

    for (i in 1:k) {
      u[i+2] <- sum(-ox[,i] * od + ox[,i] * (log(deni/gama1)/gama2))
    }

  } else {

    m <- ncol(mb)
    gama1 <- exp(-oz %x% t(mb[1,]) - ox %*% mb[3:(k+2), ])
    gama2 <- exp(-oz %x% t(mb[2,]) - ox %*% mb[3:(k+2), ])

    sm <- n:1
    dlam <- (od/sm) %x% t(rep(1, m))

    lam2 <- apply(dlam * gama2, 2, cumsum)
    p <- exp(-lam2)
    pl <- rbind(matrix(1, ncol = m), p[1:(n-1), ])
    ru <- apply(pl * dlam * gama1, 2, cumsum)/p

    deni <- gama1 + gama2 * ru

    lik <- od %*% log(deni) + colSums(log(deni/gama1)/gama2)

    u <- matrix(NA, nrow = k + 2, ncol = m)
    u[1, ] <- (oz*od) %*% (-gama1/deni) + t(oz) %*% (ru/deni)
    u[2, ] <- (oz*od) %*% (-gama2*ru/deni) - t(oz) %*% (ru/deni) + t(oz) %*% (log(deni/gama1)/gama2)

    for (i in 1:k) {
      u[i+2, ] <- -sum(ox[,i] * od) + t(ox[,i]) %*% (log(deni/gama1)/gama2)
    }

  }

  u

  return(u)
}

