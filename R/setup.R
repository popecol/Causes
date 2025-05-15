
# Setup and a set of helper functions.

library(RColorBrewer)

# high-contrast colour scheme

# Colors
col <- c('#004488', '#DDAA33', '#BB5566')

col2 <- adjustcolor(col, alpha = 1/2)
col4 <- adjustcolor(col, alpha = 1/4)
col10 <- adjustcolor(col, alpha = 1/10)
col_shaded <- adjustcolor(col, alpha = 1/35)




get_names <- function(obj) {
  # Returns a vector of desired variable names.
  # obj: an object inheriting from "gam" or "character".

  exclude <- c("plot_id", "observer_id", "fyear", "x", "y", "x,y")

  if(inherits(obj, "gam")) {
    nm <- names(obj[["var.summary"]])
  }
  else {
    if(inherits(obj, "character")) {
      nm <- sub(".*\\((.*)\\).*", "\\1", obj)
    }
  }

  nm <- nm[is.na(match(nm, exclude))]
  return(nm)
}



first_toupper <- function(x) {
  # Capitalises the first letter in a string.
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  return(x)
}



quiet <- function(x) {
  # Helper fun: returns function output without printed messages
  # Author: Hadley Wickham

  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}



part_res <- function(obj, ...) {
  UseMethod("part_res", obj)
}



part_res.gam <- function(obj, out = NULL, ...)
{
  # Calculates components needed to produce partial dependence plots.
  # Arguments:
  #   obj:  fitted 'gam' or 'bam' object
  #   out:  vector of variables to omit
  #   ...:  additional arguments passed to 'plot.gam'
  # Value: named list

  require(mgcv)
  pg <- plot.gam(obj, select = 0, n = 100, ...)
  nam <- names(obj[["var.summary"]])
  id <- !nam %in% out
  nam <- nam[id]
  pg <- pg[id]
  names(pg) <- nam
  return(pg)
}



pres_plot <- function(obj, trans = I, scale = -1, ...) {
  # Plots partial response curves for selected variables of 'gam' objects.
  # Arguments:
  #     obj:  fitted 'gam' or 'bam' object
  #   trans: function to apply before plotting
  #   scale: -1 (default): the same y-axis scale for each plot; 0: different y axis for each plot
  #   ...:  additional arguments passed to 'plot'

  # variables to exclude from plotting:
  out <- c("plot_id", "observer_id", "fyear", "x", "y")

  pr <- part_res(obj, out) # partial residuals
  pr <- lapply(pr, function(x) append(x, list(lci = x$fit - 1.96 * x$se, uci = x$fit + 1.96 * x$se)))
  pr <- lapply(pr, function(x) append(x, list(xlab_spaced = first_toupper(gsub("_", " ", x$xlab)))))

  ylim <- range(sapply(pr, function(x) range(x$lci, x$uci)))
  if(scale == -1)
    pr <- lapply(pr, function(x) append(x, list(ylim = ylim))) else
      pr <- lapply(pr, function(x) append(x, list(ylim = range(x$lci, x$uci))))

  all_names <- names(pr)
  # id <- ord(all_names)
  # all_names <- all_names[id]

  rows <- floor(sqrt(length(all_names)))
  cols <- ceiling(length(all_names) / rows)
  op <- par(mfrow = c(rows, cols), mar = c(5, 3, 1, 1))

  for (i in all_names)
  {
    with(pr[[i]], plot(x, trans(fit), type = "n", ylim = trans(ylim), xlab = xlab_spaced, ylab = "", ...))
    with(pr[[i]], polygon(c(x, rev(x)), c(trans(lci), trans(rev(uci))), border = NA, col = col4[2]))
    with(pr[[i]], lines(x, trans(fit), col = col[2], lty = 1, lwd = 2))
  }
  par(op)
}



prevalence <- function(x) sum(x > 0) / length(x)



threshold <- function(x, y)
{
  # Assignes zeros to all values in x that are not greater than y.
  x[x <= y] <- 0
  x
}



hist_comp <- function(x, y, xname = NA, cex.lab = 1.3, ...)
{
  # Plots a histogram of x and adds an outline of y as a reference.
  # Can be used to visually compare histograms.

  # Arguments:
  #     x: a vector of values for which the histogram is desired
  #     y: a vector of values for the reference histogram
  # xname: variable name (optional)
  #   ...: arguments passed to 'hist'

  hy <- hist(y, plot = FALSE, ...)
  mid <- hy$mids
  mid <- mid - (mid[2] - mid[1]) / 2
  ydens <- hy$density
  mid <- c(0, mid); ydens <- c(0, ydens)

  hx <- hist(x, plot = FALSE, ...)
  # hx <- hist(x, plot = FALSE)
  xdens <- hx$density
  hx$breaks = hy$breaks
  if(!is.na(xname)) hx$xname <- xname

  plot(hx, freq = FALSE, ylim = range(xdens, ydens), col = col4[2], border = col2[2], main = "", ylab = "Probability density", cex.lab = cex.lab)
  lines(mid, ydens, type = "s", lwd = 2, col = col2[1])
}
