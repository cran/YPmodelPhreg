fn_u <- function(mb, od, oz, ...) {

  n <- length(oz)

  if (is.null(nrow(mb))) {

    gama1 <- exp(-oz * mb[1])
    gama2 <- exp(-oz * mb[2])

    sm <- n:1
    dlam <- od/sm

    gaw2 <- dlam[1:n] * gama2
    gaw1 <- dlam[1:n] * gama1

    p <- exp(-cumsum(gaw2))
    pl <- c(1, p[1:(n-1)])
    ru <- cumsum(pl * gaw1)/p
    rul <- c(0, ru[1:(n-1)])

    deni <- gama1 + gama2 * ru
    denil <- gama1 + gama2 * rul

    lik <- sum(od * log(deni) + log(deni/gama1)/gama2)

    u1 <- sum((oz*od) * (-gama1/deni) + oz * ru/deni)
    u2 <- sum((oz*od) * (-gama2*ru/deni) - oz * ru/deni + oz * log(deni/gama1)/gama2)

    u <- c(u1, u2)

  } else {

    m <- ncol(mb)
    gama1 <- exp(-oz %x% t(mb[1,]))
    gama2 <- exp(-oz %x% t(mb[2,]))

    sm <- n:1
    dlam <- (od/sm) %x% t(rep(1, m))

    gaw2 <- dlam[1:n,] * gama2
    gaw1 <- dlam[1:n,] * gama1

    p <- exp(-apply(gaw2, 2, cumsum))
    pl <- rbind(matrix(1, ncol = m), p[1:(n-1), ])
    ru <- apply(pl * gaw1, 2, cumsum)/p
    rul <- rbind(matrix(0, ncol = m), ru[1:(n-1), ])

    deni <- gama1 + gama2 * ru
    denil <- gama1 + gama2 * rul

    lik <- od %*% log(deni) + colSums(log(deni/gama1)/gama2)

    u1 <- (oz*od) %*% (-gama1/deni) + t(oz) %*% (ru/deni)
    u2 <- (oz*od) %*% (-gama2*ru/deni) - t(oz) %*% (ru/deni) + t(oz) %*% (log(deni/gama1)/gama2)

    u <- rbind(u1, u2)

  }

  return(u)
}

