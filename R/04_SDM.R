
# Fitting a Species Distribution Model to the VE data

# Setup -------------------------------------------------------------------

library(mgcv)
library(parallel)
library(dplyr)
library(pbapply)

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
dim(ve)

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
ve_data <- data.frame(plot_id = pred_data$plot_id, ve)

# Averaging VE data over sites
ve_data <- ve_data %>%
  group_by(plot_id) %>%
  summarize_all(mean, na.rm = TRUE) %>%
  arrange(plot_id)

ve_data <- select(ve_data, -plot_id)
ve_data <- as.data.frame(ve_data)
ve_data <- round(ve_data, 2)

# Preparing environmental data
v  <- names(full$var.summary)
v <- v[!v %in% c("log_dens_lag", "year", "observer_id", "fyear", "x", "y")]
data <- data[, v]
# apply(data, 2, \(x) length(unique(x)))

# Averaging environmental data over sites
env_data <- data %>%
  group_by(plot_id) %>%
  summarize_all(mean, na.rm = TRUE) %>%
  arrange(plot_id)

env_data <- select(env_data, -plot_id)
env_data <- as.data.frame(env_data)
env_data <- signif(env_data, 2)
summary(env_data)
# apply(env_data, 2, \(x) length(unique(x)))

# Variable names
v <- names(env_data)


# Fitting the model to VE data --------------------------------------------

R <- 100 # no. of replicates of the experiment
k <- 6
gamma <- 4
f <- formula(paste0("VE ~ ", paste0("s(", v, ", k = k, bs = 'cr')", collapse = " + ")))

fdat <- env_data
clusterExport(cl, c("fdat", "ve_data", "get_names", "v", "k", "gamma", "f"))

# control <- list(trace = TRUE)

fun <- function(i, p = 0) {
  fdat$VE <<- ve_data[, i]

  sim_full <- bam(f, fdat, family = tw(), discrete = TRUE)
  # summary(sim_full); plot(sim_full, res = TRUE, pages = 1, pch = 21)

  if(p > 0) {
    sim_fit <- try(bam(f, fdat, family = tw(), discrete = TRUE, gamma = gamma, select = TRUE))
    if(inherits(sim_fit, "try-error"))
      sim_fit <- bam(f, fdat, family = tw(), discrete = TRUE, gamma = gamma, select = TRUE, method = "NCV")
  } else
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
save(sim, file = paste0(path, "/sim_SDM_0.RData"))

# sim <- pblapply(1:R, fun, p = 0.05, cl = NULL)
sim <- pblapply(1:R, fun, p = 0.05, cl = cl)
sim <- simplify2array(sim, higher = FALSE)
rownames(sim) <- v
save(sim, file = paste0(path, "/sim_SDM_005.RData"))


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
