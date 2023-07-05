plot.ypreg <- function(x, legend.loc = "bottomleft", lty = 1, lty.pci = 3, lty.sci = 2, col = "black", col.pci = "blue3", col.sci = "firebrick3",  lwd = 2, lwd.pci = 2, lwd.sci = 2, ...) {

  match.call()
  dotlist <- list(...)

  digit <- paste("%.", max(2, getOption("digits") - 5), "f", sep = "")

  hr <- x$res_summ$hr
  yt <- x$res_summ$yt
  nwt <- x$res_summ$nwt
  ca22 <- x$res_summ$ca22
  zalpha <- x$res_summ$zalpha
  ld2 <- x$res_summ$ld2
  ud2 <- x$res_summ$ud2
  Lr <- x$res_summ$Lr
  Ur <- x$res_summ$Ur
  Lv <- sprintf(digit, Lr)
  Uv <- sprintf(digit, Ur)

  n <- x$n
  alpha <- x$alpha

  # simultaneous CI (log-transformed)
  upp22 <- exp(ca22/sqrt(n)*nwt/hr) * hr;
  low22 <- exp(-ca22/sqrt(n)*nwt/hr) * hr;

  # point-wise CI (log-transformed)
  upp3 <- exp(zalpha/sqrt(n) * nwt/hr) * hr
  low3 <- exp(-zalpha/sqrt(n) * nwt/hr) * hr

  ld3 <- sum(yt <= min(yt))
  ud3 <- sum(yt <= max(yt))

  xmax <- round(max(yt), 2)
  ymax <- round(max(c(upp22, upp3)) * 1.2, 2)

  ci.text <- paste(round((1 - alpha) * 100, 2), "% Pointwise and Simultaneous Confidence Intervals", sep = "")
  pci.shorttext <- paste(round((1 - alpha) * 100, 2), "% pointwise CI", sep = "")
  sci.shorttext <- paste(round((1 - alpha) * 100, 2), "% simultaneous CI", " on [L, U] = [", Lv, ", ", Uv, "]", sep = "")
  main.text <- paste("Hazard Ratio with ", ci.text, collapse = "")

  dotlist$lty <- lty
  dotlist$lwd <- lwd
  dotlist$col <- col
  if(is.null(dotlist$xlab)) dotlist$xlab = "Time"
  if(is.null(dotlist$ylab)) dotlist$ylab = "Hazard Ratio"
  if(is.null(dotlist$main)) dotlist$main = main.text
  if(is.null(dotlist$xlim)) dotlist$xlim = c(0, xmax)
  if(is.null(dotlist$ylim)) dotlist$ylim = c(0, ymax)
  if(is.null(dotlist$type)) dotlist$type = 'l'
  if(is.null(dotlist$axes)) dotlist$axes = FALSE
  if(is.null(dotlist$xaxs)) dotlist$xaxs = "i"
  if(is.null(dotlist$yaxs)) dotlist$yaxs = "i"
  dotlist$x = yt
  dotlist$y = hr

  do.call(what = plot, args = dotlist)
  polygon(x = c(yt[ld2:ud2], yt[ud2:ld2]),
          y = c(low22[ld2:ud2], upp22[ud2:ld2]), lty = 0, col = rgb(0.70,0.30,0.20, 0.1))

  abline(h = 1, lty = 2, col = "black")

  lines(yt[ld3:ud3], upp3[ld3:ud3], type = 'l', lty = lty.pci, lwd = lwd.pci, col = col.pci)
  lines(yt[ld3:ud3], low3[ld3:ud3], type = 'l', lty = lty.pci, lwd = lwd.pci, col = col.pci)
  lines(yt[ld2:ud2], upp22[ld2:ud2], type = 'l', lty = lty.sci, lwd = lwd.sci, col = col.sci)
  lines(yt[ld2:ud2], low22[ld2:ud2], type = 'l', lty = lty.sci, lwd = lwd.sci, col = col.sci)

  x.max <- dotlist$xlim[2]
  y.max <- dotlist$ylim[2]

  if(dotlist$axes == FALSE) {
    x.axis0 <- signif(seq(from = 0, to = x.max, length.out = 6),2)
    x.axis <- c(x.axis0[x.axis0 <= signif(x.max, 2)], signif(x.max, 2))
    y.axis0 <- signif(seq(from = 0, to = y.max, length.out = 6),2)
    y.axis <- c(y.axis0[y.axis0 <= signif(y.max, 2)], signif(y.max, 2))
    axis(1, x.axis)
    axis(2, c(1, y.axis))
  }

  legend(legend.loc, legend = c("hazard ratio", pci.shorttext, sci.shorttext),
         col = c(dotlist$col, col.pci, col.sci), lty = c(lty, lty.pci, lty.sci),
         lwd = c(lwd, lwd.pci, lwd.sci), box.lty = 0)

}
