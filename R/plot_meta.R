
plot_meta <- function(fitted_model, ran = FALSE, points = FALSE, ...) {
  # Plots results of multi-species simulations.
  # Whiskers represent 95% confidence intervals. Lower and upper hinges 
  # are 50% confidence intervals and the middle line is the predicted mean. 
  
  # Arguments:
  # fitted_model: a model fitted using 'glmmTMB'
  #          ran: show random effects of species?
  #       points: show raw data points?
  #         ... : arguments passed to 'plot' and 'ci' functions
  
  require(glmmTMB)
  
  ci <- function(fit, se, level = 0.95, type = c("link", "response"), large_se = FALSE, ...) {
    # Calculates confidence intervals from the fitted values and standard errors.
    # If all values are identical for a given combination of factors, then
    # the error bar ranges from 0 to 1. If lagrge_se = FALSE, it will not be drawn.
    
    if(!large_se)
      se <- ifelse(se > 10, 0, se)
    type <- match.arg(type)
    trans <- switch(type, link = I, response = plogis)
    q <- qnorm((1 - level) / 2, lower.tail = FALSE)
    lci <- fit - q * se
    uci <- fit + q * se
    df <- cbind(lci, uci)
    trans(df)
  }
  
  data <- fitted_model[["frame"]]
  mod <- levels(data$mod)
  sel <- levels(data$sel)
  newdat <- expand.grid(mod = mod, sel = sel, species = NA, w = 20)
  pre <- predict(fitted_model, newdat, se = TRUE)
  fit <- as.vector(pre$fit); se <- as.vector(pre$se.fit)
  if(ran) {
    # Add the random effect variance (i.e. for a species) to the model prediction variance.
    v <- as.numeric(VarCorr(fitted_model)[["cond"]][["species"]])
    se <- sqrt(se^2 + v)
  }
  ci95 <- ci(fit, se, 0.95, "response", ...)
  ci50 <- ci(fit, se, 0.5, "response", ...)
  newdat <- cbind(newdat[, 1:2], plogis(fit), ci50, ci95)
  names(newdat) <- c("mod", "sel", "fit", "lci50", "uci50", "lci95", "uci95")
  
  f <- formula(paste(names(data)[1], "~ mod"))
  
  x <- seq(length(mod))
  w <- rep(0.2, 2); w[1] <- -w[1]
  xlim <- c(min(x) - 0.6, max(x) + 0.6)
  
  # if(points) {
  #   ylim <- range(data[1], 0.5)
  # } else {
  #   ylim <- range(newdat[c("lci95", "uci95")], 0.5)
  # }
  ylim <- c(0, 1)
  
  cex.lab <- 1.8
  
  plot(x, x, xlim = xlim, ylim = ylim, type = "n", xaxt = "n", xlab = "", cex.lab = cex.lab, cex.axis = 1.5, ...)
  axis(1, at = seq(mod), labels = mod, las = 2, cex.axis = 1.5)
  mtext("Method", 1, 5, cex = cex.lab * 0.7)
  abline(h = 0.5, col = adjustcolor("black", alpha = 0.3))
  
  for (i in seq(length(sel))) {
    xw <- x + w[i]
    nd <- newdat[newdat$sel == sel[i], ]
    dat <- data[data$sel == sel[i], ]
    const <- names(which(tapply(dat[[1]], dat$mod, sd) == 0))
    if(length(const) > 0) {
      idx <- dat$mod %in% const
      # idx[1:100] <- FALSE
      # summary(idx)
      # dat[idx, 1] <- NA
    }
    
    if(points) {
      require(beeswarm)
      bs <- beeswarm(f, dat, method = "swarm", spacing = 0.01, col = col_shaded[i], do.plot = FALSE)
      bs$x <- bs$x + w[i]
      points(bs, col = col_shaded[i])
    }
    
    # segments(xw, nd$lci95, xw, nd$uci95, lwd = 6)
    # segments(xw, nd$lci50, xw, nd$uci50, lwd = 12)
    # points(xw, nd$fit, pch = "-")
    # points(xw, nd$fit, pch = 20, col = "white")
    
    bp <- boxplot(f, dat, plot = FALSE)
    bp[["stats"]][1, ] <- nd[["lci95"]]
    bp[["stats"]][2, ] <- nd[["lci50"]]
    bp[["stats"]][3, ] <- nd[["fit"]]
    bp[["stats"]][4, ] <- nd[["uci50"]]
    bp[["stats"]][5, ] <- nd[["uci95"]]
    bxp(bp, add = TRUE, boxwex = 0.25, staplewex = 0, boxfill = col4[i], medlwd = 4, 
        medcol = col[i], whisklty = 1, xaxt = "n", yaxt = "n", at = xw, outline = FALSE)
  }
  return(invisible(newdat))
}
