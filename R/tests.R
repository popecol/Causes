
# A set of functions to evaluate the models.

library(RColorBrewer)
# display.brewer.all(4, type = "qual")

test <- function(true, sim, p = 0.05) {
  var_left <- !is.na(sim) # variables remained after preliminary screening (based on correlations with log-density)
  signif <-  (sim < p) & var_left # significant variables in the GAM model (at the 'p' level)
  # correct <- true == signif
  # a <- data.frame(true = true, sim = sim[, 1], var_left = var_left[, 1], signif = signif[, 1], correct = correct[, 1])
  signif
}


evaluate_model <- function(true, sim, p = 0.05) {
  # Evaluates a model (the efficiency of variable selection)
  
  reference <- factor(true != 0, levels = c(TRUE, FALSE)) # conversion to factor
  
  signif <- as.data.frame((sim < p) & (!is.na(sim))) # variables selected as significant at the p-level.
  signif[] <- lapply(signif, \(x) factor(x, levels = c(TRUE, FALSE))) # conversion to factor
  tab <- lapply(signif, \(x) table(reference, x)) # confronting with the truth
  
  # [1, 1] true positives
  # [2, 1] false positives
  # [1, 2] false negatives
  # [2, 2] true negatives
  
  # Accuracy = correct classification rate
  # Probability of ((detecting an effect if there is one) OR (detecting no effect if there is none)).
  acc <- sapply(tab, \(x) 100 * (x[1, 1] + x[2, 2]) / sum(x))
  
  # Sensitivity (true positive rate = power of a test)
  # Probability of detecting an effect if there is one. 
  # Probability of correct classification of a factor that has an influence on population changes.
  tpr <- sapply(tab, \(x) 100 * (x[1, 1] / (x[1, 1] + x[1, 2])))
  
  # Specificity (true negative rate)
  # Probability of detecting no effect if there is none.
  # Probability of correctly classifying a factor that has no influence on population changes.
  tnr <- sapply(tab, \(x) 100 * (x[2, 2] / (x[2, 2] + x[2, 1])))
  
  # Balanced accuracy (accuracy estimate for imbalanced data)
  ba <- (tpr + tnr) / 2
  
  data.frame(acc, tpr, tnr, ba)
}


plt_multi <- function(m, measure = c("acc", "tpr", "tnr", "ba"), bw = 3, dy = 0.1, yaxis = FALSE) {
  # Plots the results of the model evaluation for many species
  # Arguments:
  #          m: list of lists of model evaluations for each species (with elements created by the 'evaluate_model' function)
  #    measure: specific measure (see: 'evaluate_model')
  #         bw: bandwidth for density distribution estimation
  #         dy: scaling factor defining the vertical spacing of the density plots
  #      yaxis: plot y-axis labels?
  
  # Invisibly returns a list of the kernel density distributions.
  
  measure <- match.arg(measure) 
  k <- length(m[[i]])
  
  for(i in 1:length(m)) {
    p <- lapply(m[[i]], \(x) x[[measure]])
    d <- lapply(p, \(x) density(x, bw = bw, from = 0, to = 100))
    y <- 1:k * dy
    xlab <- switch(measure,
                   acc = "Accuracy",
                   tpr = "Sensitivity",
                   tnr = "Specificity",
                   ba  = "Balanced accuracy")
    
    for(a in 1:length(d)) {
      xd <- d[[a]][["x"]]  
      yd <- d[[a]][["y"]]
      yd[1] <- yd[length(yd)] <- 0
      yd <- yd + (a - 1) * dy
      
      if(i == 1 && a == 1) {  
        plot(xd, yd, type = "n", ylim = c(0, k * dy - (dy / k * 1.5)), yaxt = "n", xlab = xlab, ylab = "", cex.lab = 1.8, xlim = c(0, 100)) 
        grid(ny = NA)
      }
      
      polygon(xd, yd, col = col4[a], border = NA)
      yd[i] <- yd[length(yd)] <- NA
      lines(xd, yd, col = "black")
    }
    
    if(yaxis)
      axis(2, at = y - dy, labels = lab, las = 1, cex.axis = 1.5)
    else
      axis(2, at = y - dy, labels = FALSE)
    
    axis(4, at = y - dy, labels = FALSE)
  }
  return(invisible(d))
  
}  



plt <- function(m, measure = c("acc", "tpr", "tnr", "ba"), bw = 3, dy = 0.1, yaxis = FALSE, plot.me = TRUE) {
  # Plots the results of the model evaluation.
  # Arguments:
  #          m: list of model evaluations (with elements created by the 'evaluate_model' function)
  #    measure: specific measure (see: 'evaluate_model')
  #         bw: bandwidth for density distribution estimation
  #         dy: scaling factor defining the vertical spacing of the density plots
  #      yaxix: plot y-axis labels?
  #    plot.me: plot result?
  # Invisibly returns a list of the kernel density distributions.
  
  measure <- match.arg(measure)
  p <- lapply(m, \(x) x[[measure]])
  d <- lapply(p, \(x) density(x, bw = bw, from = 0, to = 100))
  
  if(plot.me) {
    k <- length(m)
    y <- 1:k * dy
    xlab <- switch(measure,
                   acc = "Accuracy",
                   tpr = "Sensitivity",
                   tnr = "Specificity",
                   ba  = "Balanced accuracy")
    
    xd <- d[[1]][["x"]];  yd <- d[[1]][["y"]]
    yd[1] <- yd[length(yd)] <- 0
    plot(xd, yd, type = "n", ylim = c(0, k * dy - (dy / k * 1.5)), yaxt = "n", xlab = xlab, ylab = "", cex.lab = 0.3, xlim = c(0, 100))
    grid(ny = NA)
    polygon(xd, yd, col = col[1], border = NA)
    yd[1] <- yd[length(yd)] <- NA
    lines(xd, yd, col = sil)
    for (i in 2:k) {
      xd <- d[[i]][["x"]]
      yd <- d[[i]][["y"]]
      yd[1] <- yd[length(yd)] <- 0
      yd <- yd + (i - 1) * dy
      polygon(xd, yd, col = col[i], border = NA)
      yd[1] <- yd[length(yd)] <- NA
      lines(xd, yd, col = sil)
      if(yaxis)
        axis(2, at = y - dy, labels = lab, las = 1, cex.axis = 1.5, cex.lab = 0.3)
      else
        axis(2, at = y - dy, labels = FALSE)
      axis(4, at = y - dy, labels = FALSE)
    }
  }
  return(invisible(d))
}


plot_ctr <- function(m, measure = c("acc", "tpr", "tnr", "ba"), yaxis = FALSE) {
  # Plots contrasts (pairwise differences between measures of variable selection efficiency).
  # Arguments:
  #          m: list of model evaluations (with elements created by the 'evaluate_model' function)
  #    measure: specific measure (see: 'evaluate_model')
  
  between <- function(x, lower, upper) x > lower & x < upper
  
  measure <- match.arg(measure)
  p <- sapply(m, \(x) x[[measure]])
  pairs <- combn(length(m), 2)
  
  id <- outer(lab_short, lab_short, "paste", sep = " - ")
  idx <- lower.tri(id)
  ctr <- id[idx]
  k <- length(ctr)
  
  ct <- matrix(NA, nrow(p), k)
  for (i in 1:k) {
    v1 <- pairs[1, i]
    v2 <- pairs[2, i]
    ct[, i] <- p[, v2] - p[, v1]
  }
  ct <- data.frame(ct)
  names(ct) <- ctr
  
  mea <- apply(ct, 2, mean)
  ci_50 <- apply(ct, 2, dhdi, credMass = 0.5)
  ci_95 <- apply(ct, 2, dhdi, credMass = 0.95)
  ci_99 <- apply(ct, 2, dhdi, credMass = 0.99)
  
  pd <- data.frame(t(rbind(mea, ci_50, ci_95, ci_99)))
  names(pd) <- c("m", "lci_50", "uci_50", "lci_95", "uci_95", "lci_99", "uci_99")
  pd <- transform(pd, zero = as.numeric(between(0, lci_95, uci_95)) + 1)
  
  xlab <- switch(measure,
                 acc = "accuracy",
                 tpr = "sensitivity",
                 tnr = "specificity",
                 ba  = "balanced accuracy")
  # xlab <- paste("Difference in", xlab)
  y <- 1:k
  
  plot(y ~ m, pd, type = "n", xlim = range(pd), ylim = c(0.7, k + 0.3), yaxt = "n", xlab = xlab, ylab = "", cex.lab = 1.8)
  grid(ny = NA)
  abline(v = 0, col = "grey30")
  segments(pd$lci_95, y, pd$uci_95, y, col = col2[1:2][pd$zero], lwd = 4)
  segments(pd$lci_50, y, pd$uci_50, y, col = col2[1:2][pd$zero], lwd = 6)
  points(pd$m, y, pch = 16, cex = 2.1, col = "white")
  points(pd$m, y, pch = 16, cex = 1.7, col = pal[1:2][pd$zero])
  points(pd$m, y + 0.3, pch = "*", cex = 2, col = c(1, NA)[pd$zero])
  # points(pd$m, y + 0.3, pch = "**", cex = 2, col = c(1, NA)[pd$z])
  
  if(yaxis)
    axis(2, at = y, labels = ctr, las = 1, cex.axis = 1.5)
  else
    axis(2, at = y, labels = FALSE)
  axis(4, at = y, labels = FALSE)
  mtext("Difference in:", side = 1, outer = TRUE, cex = 1.2, adj = -0.25, padj = -2)
  
  return(invisible(pd))
}


dhdi <- function(x, credMass = 0.95, ...) {
  # Calculate the highest density interval for a kernel density distribution.
  d <- density(x, ...)
  ci <- suppressWarnings(HDInterval::hdi(d, credMass = credMass))
}
