print.ypreg <- function(x, ...) {

  ci_level <- round((1 - x$alpha) * 100, 2)
  digit <- paste("%.", max(3, getOption("digits") - 4), "f", sep = "")
  digitp <- paste("%.", max(3, getOption("digits") - 3), "f", sep = "")

  display_p <- function(sv) ifelse(as.numeric(sv) < 0.0001, "<0.0001", sv)

  xnames <- x$xnames
  n <- x$n
  p <- x$p
  alpha <- x$alpha

  # tab1
  best <- x$best_ypx
  exbcip <- x$res_summ$exbcip
  expb <- exbcip[,1]
  lower.expb <- exbcip[,2]
  upper.expb <- exbcip[,3]
  SE <- exbcip[,4]
  z <- exbcip[,5]
  pvalue <- exbcip[,6]

  best <- sprintf(digit, best)
  expb <- sprintf(digit, expb)
  lower.expb <- sprintf(digit, lower.expb)
  upper.expb <- sprintf(digit, upper.expb)
  expbCI <-  paste("(", lower.expb, ", ", upper.expb, ")", sep = "")
  SE <- sprintf(digit, SE)
  z <- sprintf(digit, z)
  pvalue <- display_p(sprintf(digitp, pvalue))

  tab1 <- data.frame(best, SE, expb, expbCI, z, pvalue)
  colnames(tab1) <- c("coef", "se(coef)", "exp(coef)", paste((1 - alpha) * 100, "% CI.exp(coef)", sep = ""), "z", "p")
  rownames(tab1) <- c("b1", "b2", xnames)

  #tab2
  avhrcip <- x$res_summ$avhrcip
  avhr <- sprintf(digit, avhrcip[1])
  SE.avhr <- sprintf(digit, avhrcip[2])
  lower.avhr <- sprintf(digit, avhrcip[3])
  upper.avhr <- sprintf(digit, avhrcip[4])
  avhrCI <-  paste("(", lower.avhr, ", ", upper.avhr, ")", sep = "")
  z.avhr <- sprintf(digit, avhrcip[5])
  pvalue.avhr <- display_p(sprintf(digitp, avhrcip[6]))

  data.frame(avhr, avhrCI)
  tab2 <- data.frame(avhr, SE.avhr, avhrCI, z.avhr, pvalue.avhr)
  colnames(tab2) <- c("avhr", "se(avhr)", paste((1 - alpha) * 100, "% CI", sep = ""), "z", "p*")
  rownames(tab2) <- ""

  time.hr <- x$time.hr
  if (!is.null(time.hr)) {
    hrci <- x$res_hrci
    tab3 <- NULL
    for (u in 1:nrow(hrci)) {
      # y.hr <- sprintf(digit, hrci[u, 1])
      hr <- sprintf(digit, hrci[u, 2])
      lowerIC.hr <- sprintf(digit, hrci[u, 3])
      upperIC.hr <- sprintf(digit, hrci[u, 4])
      CI.hr <-  paste("(", lowerIC.hr, ", ", upperIC.hr, ")", sep = "")
      tab3 <- rbind(tab3, c(hr, CI.hr))
    }
    tab3 <- matrix(tab3, ncol = 2)
    tab3 <- as.data.frame(tab3)
    colnames(tab3) <- c("hr", paste((1 - alpha) * 100, "% CI", sep = ""))
    rownames(tab3) <- paste("time=", time.hr, sep = "")
  }

  cat("I. -------------- < Estimation Results > ----------------------\n\n")
  print(tab1)
  cat("---------------------------------------------------------------\n\n")

  cat("II. ------------- < The average hazard ratio > ----------------\n\n")
  print(tab2)
  cat("\n")
  cat("* two-sided p-value\n")
  cat("---------------------------------------------------------------\n\n")
  cat("III. ------------ < Hazard ratios at time t > -----------------\n\n")
  if (!is.null(time.hr)) {
    print(tab3)
  } else {
    cat(" No time point is entered.\n")
    cat(" Please input some time points into the 'time.hr' argument.\n")
  }
  cat("---------------------------------------------------------------\n")


}
