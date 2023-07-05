ypreg.default <- function(data, alpha = 0.05, time.hr = NULL, L = NULL, U = NULL, repnum = 5000, tau = NULL, ...) {

  match.call()

  y <- as.numeric(data[, 1])
  d <- as.numeric(data[, 2])
  z <- as.numeric(data[, 3])
  x <- as.matrix(data[, -c(1:3)])

  if (ncol(x) == 0) stop("Please supply at least one X variable to data")
  xnames <- colnames(x)

  n <- length(y)
  p <- ncol(x)

  I <- order(y)
  oy <- y[I]
  od <- d[I]
  oz <- z[I]
  ox <- x[I,]
  owz <- cbind(oz, ox)

  k <- ncol(owz)

  mean.x <- colMeans(x)
  sd.x <- apply(x, 2, sd)

  if (any(c(mean.x > 1, sd.x > 1.5))) stop("Please standadize covariates X first")

  # cox regression
  fit_coxph <- coxph(Surv(oy, od) ~ owz)
  best_cox <- fit_coxph$coefficients

  res0_yp0 <- fn_est0_ypmodel0(oy, od, oz, b0 = c(0, 0))
  ini_b0 <- res0_yp0$minb

  res_yp0 <- fn_est_ypmodel0(oy, od, oz, ini_b0)
  best_b0 <- res_yp0$minb

  ini_b <- c(best_b0, best_cox[-1])
  res_ypx <- fn_est_ypmodelx(oy, od, owz, ini_b)
  best_ypx <- res_ypx$minb
  res_summ <- fn_summ(best_ypx, oy, od, owz, alpha, L, U, h = 10E-8, repnum, tau)

  if (!is.null(time.hr)) {
    res_hrci <- fn_hrci(time.hr, res_summ)
  } else {
    res_hrci <- NULL
  }

  results <- c()
  results$input.data <- data
  results$fit_coxph <- fit_coxph
  results$ini_b0 <- ini_b0
  results$best_b0 <- best_b0
  results$best_ypx <- best_ypx
  results$res_summ <- res_summ
  results$res_hrci <- res_hrci
  results$tau <- tau
  results$xnames <- xnames
  results$p <- p
  results$n <- n
  results$time.hr <- time.hr
  results$alpha <- alpha

  class(results) <- "ypreg"

  results

}

