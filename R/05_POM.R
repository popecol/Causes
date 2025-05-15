
# Fitting a process-oriented model to the VE data


# Setup -------------------------------------------------------------------

library(mgcv)
library(dplyr)
library(pbapply)
library(parallel)

source("R/setup.R")
source("R/set_species.R")
source("R/tests.R")

ncores <- 5L
cl <- makeCluster(ncores)

clusterEvalQ(cl, {
  library(mgcv)
})


# Loading data and models --------------------------------------------------------------------

species_info()
path <- paste0("data/", species)

# The VE data
load(paste0(path, "/ve_prim.RData"))
ve <- ve_prim; rm(ve_prim)
ve_data <- as.data.frame(ve)

# The full model
load(paste0(path, "/gamm_full_inter.RData"))
# summary(full, re.test = FALSE)

# Selected model
load(paste0(path, "/gamm_selected_inter.RData"))
# summary(fit, re.test = FALSE)

# The real data
load(paste0(path, "/data_inter.RData"))
# summary(data)

# Predictions from the selected model
load(paste0(path, "/pred_data.RData"))


# Preparing data ----------------------------------------------------------

# Aligning data
idx <- match(pred_data[["id_year"]], data[["id_year"]])
data <- data[idx, ]

# Variable names
v <- get_names(full)
v <- v[!v %in% "log_dens_lag"]
# v <- v[ord(v)]
v <- v[!v %in% "year"]


# Fitting a model ---------------------------------------------------------

R <- 100 # no. of replicates of the experiment
k <- 6
gamma <- 4
f <- update(formula(full), VE ~ .)

fdat <- data
clusterExport(cl, c("fdat", "ve_data", "get_names", "v", "k", "f"))

fun <- function(i, p = 0) {
  fdat$VE <<- ve_data[, i]

  # sim_full <- gam(f, fdat, family = tw())
  sim_full <- bam(f, fdat, family = tw(), discrete = TRUE)
  # summary(sim_full); plot(sim_full, res = TRUE, pages = 1, pch = 21)

  if(p > 0)
    sim_fit <- bam(f, fdat, family = tw(), discrete = TRUE, gamma = gamma, select = TRUE)
  else
    sim_fit <- sim_full

  # summary(sim_fit); plot(sim_fit, res = TRUE, pages = 1, pch = 21)
  su <- summary(sim_fit)[["s.table"]]
  if(!is.null(su)) {
    idx <- match(v, get_names(rownames(su)))
    ret <- su[idx, ][, "p-value"]
  } else {
    ret <- matrix(NA, nrow = length(v), ncol = 1)
  }
  return(ret)
}


sim <- pblapply(1:R, fun, p = 0, cl = cl)
sim <- simplify2array(sim, higher = FALSE)
rownames(sim) <- v
save(sim, file = paste0(path, "/sim_POM_0.RData"))

# sim <- pblapply(1:R, fun, p = 0.05, cl = cl) # To u mnie nie działa.
sim <- pblapply(1:R, fun, p = 0.05) # A to działa.
sim <- simplify2array(sim, higher = FALSE)
rownames(sim) <- v
save(sim, file = paste0(path, "/sim_POM_005.RData"))


# Evaluation ----------------------------------------------------------

# Variables in the generating model
ref <- summary(fit, re.test = FALSE)[["s.table"]]
idx <- match(v, get_names(rownames(ref)))
ref <- ref[idx, ]
rownames(ref) <- v
true <- !is.na(ref[, 4])

M <- evaluate_model(true, sim)
summary(M)

# Clean up  -----------------------------------------------------------
stopCluster(cl)
