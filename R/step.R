
# Backward model selection for 'gam' objects.

backward <- function(obj, p = 0.05, nvar = NULL, re.test = FALSE, verbose = TRUE, gc = TRUE, ...)
{
  cond <- function(x) {
    if(is.null(nvar))
      ret <- any(x$p.value > p) & nrow(x) > 0
    else
      ret <- nrow(x) > nvar
    return(ret)
  }
  
  op <- options("warn")
  options(warn = -1)
  
  su <- summary(obj, re.test = re.test)
  A <- to.remove(su)
  
  while (cond(A))
  {
    idx <- which.max(A$p.value)
    pval <- prettyNum(A$p.value, digits  = 4)[idx]
    dv <- rownames(A)[idx]
    if(verbose) cat("Removing", dv, "     p-value =", pval, "\n")
    dterm <- paste0(". ~ . - ", paste0(substr(dv, 1, nchar(dv) - 1), ", k = k, bs = 'cr')"))
    newf <- update(formula(obj), dterm)
    if(gc) gc(verbose = FALSE)
    obj <- update(obj, newf, ...)
    (su <- summary(obj, re.test = re.test))
    if (su$m != 0) A <- to.remove(su) else A <- data.frame()
    if(nrow(A) == 0) obj <- update(obj, formula(VE ~ 1)) 
  }
  
  if(verbose) {
    if(is.null(nvar))
      cat("All p-values below", p, "\n")
    else
      cat("No. of variables retained in the final model: ", nvar, "\n")
  }
  
  options(op)
  return(obj)
}


to.remove <- function(su) {
  # Helper function extracting p-values from a 'summary.gam' object, 
  # excluding random factors, time and space.
  
  A <- data.frame(su$s.table)
  A <- A[!grepl("(plot_id|observer_id|year|fyear|(x,y))", rownames(A)), ]
  A <- transform(A, edf = round(edf, 4), p.value = round(p.value, 4))
}
