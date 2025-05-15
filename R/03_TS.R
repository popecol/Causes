
# Fitting the Time Series model to the VE data


# Setup -------------------------------------------------------------------

library(mgcv)
library(parallel)
library(dplyr)
library(pbapply)

source("R/setup.R")
# source("R/step.R")
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

# The full model
load(paste0(path, "/gamm_full_inter.RData"))
# summary(full, re.test = FALSE)
# pres_plot(full, scale = 0)

# Selected model
load(paste0(path, "/gamm_selected_inter.RData"))
# summary(fit, re.test = FALSE)
# pres_plot(fit, cex.lab = 1.2)
pres_plot(fit, scale = 0, cex.lab = 1.2)

# The real data
load(paste0(path, "/data_inter.RData"))
# summary(data)

# Predictions from the selected model
load(paste0(path, "/pred_data.RData"))


# Preparing data ----------------------------------------------------------

# Aligning data
idx <- match(pred_data[["id_year"]], data[["id_year"]])
data <- data[idx, ]

# Combining data
ve_data <- data.frame(year = pred_data$year, ve)

# Averaging VE data over years
ve_data <- ve_data %>%
  group_by(year) %>%
  summarize_all(mean, na.rm = TRUE) %>%
  arrange(year)

ve_data <- select(ve_data, -year)
ve_data <- as.data.frame(ve_data)

# Preparing environmental data
v  <- names(full$var.summary)
v <- v[!v %in% c("log_dens_lag", "plot_id", "observer_id", "fyear", "x", "y")]
# v <- v[ord(v)]
data <- data[, v]

# Averaging environmental data over years
env_data <- data %>%
  group_by(year) %>%
  summarize_all(mean, na.rm = TRUE) %>%
  arrange(year)

env_data <- select(env_data, -year)
env_data <- as.data.frame(env_data)
summary(env_data)

# Variable names
v <- names(env_data)

# Fitting the model to replicated VE data ----------------------------------

R <- 100 # no. of replicates of the experiment
k <- 6
gamma <- 4
f <- formula(paste0("VE ~ ", paste0("s(", v, ", k = k, bs = 'cr')", collapse = " + ")))

fdat <- env_data
clusterExport(cl, c("fdat", "ve_data", "get_names", "v", "k", "gamma", "f"))

fun <- function(i, p = 0) {
  fdat$VE <<- ve_data[, i]

  sim_full <- gam(f, fdat, family = tw())
  # sim_full <- bam(f, fdat, family = tw(), discrete = TRUE)
  # summary(sim_full); plot(sim_full, res = TRUE, pages = 1, pch = 21)

  if(p > 0)
    sim_fit <- gam(f, fdat, family = tw(), gamma = gamma, select = TRUE)
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


# sim <- pblapply(1:2, fun, p = 0) # test
sim <- pblapply(1:R, fun, p = 0, cl = cl)
sim <- simplify2array(sim, higher = FALSE)
rownames(sim) <- v
save(sim, file = paste0(path, "/sim_TS_0.RData"))

# sim <- pblapply(1:2, fun, p = 0.05) # test
sim <- pblapply(1:R, fun, p = 0.05, cl = cl)
sim <- simplify2array(sim, higher = FALSE)
rownames(sim) <- v
save(sim, file = paste0(path, "/sim_TS_005.RData"))


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
