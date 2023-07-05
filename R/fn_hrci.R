fn_hrci <- function(time.hr, x, ...) {

  hr <- x$hr
  yt <- x$yt
  nwt <- x$nwt
  za <- x$zalpha

  n <- x$n

  if (all(min(time.hr) < yt) | all(max(time.hr) > yt)) stop("Please set time.hr to be values within the observed data")

  pp <- length(time.hr)

  all_hrci <- NULL
  for (i in 1:pp) {
    nk1 <- sum(yt < time.hr[i])
    hrci <- c(yt[nk1], hr[nk1], exp(-za/sqrt(n)*nwt[nk1]/hr[nk1]) * hr[nk1],
               exp(za/sqrt(n)*nwt[nk1]/hr[nk1])*hr[nk1])
    all_hrci <- rbind(all_hrci, hrci)
  }

  all_hrci <- matrix(all_hrci, ncol = 4)
  colnames(all_hrci) <- c("y", "hr", "lowerIC", "upperIC")
  rownames(all_hrci) <- paste("t=", time.hr, sep = "")

  all_hrci

}

