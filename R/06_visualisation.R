
# Visualisation of models performance

# Setup -------------------------------------------------------------------

library(HDInterval)

source("R/setup.R")
source("R/tests.R")

# Assigning labels  --------------------------------

lab <- c("Time Series", "Species\nDistribution Model", "Process-Oriented\n Model")
lab_short <- c("TS", "SDM", "POM")
species <- c( "KS", "FH", "SA", "SU", "SB", "SC", "Z", "KT", "E", "PP")

# Loading simulation results ----------------------------------------------

m0 <- m5 <- list()

for (i in species){

      path <- paste0("data/", i)
      
      # The full model
      load(paste0(path, "/gamm_full_inter.RData"))
      # summary(full, re.test = FALSE)
      
      # Variable names
      v <- get_names(full)
      v <- v[!v %in% "log_dens_lag"]
      # v <- v[ord(v)]
      v <- v[!v %in% "year"]
      
      # The selected model (model used for generating a VS)
      load(paste0(path, "/gamm_selected_inter.RData"))
      
      # Variables in the generating model
      ref <- summary(fit, re.test = FALSE)[["s.table"]]
      idx <- match(v, get_names(rownames(ref)))
      ref <- ref[idx, ]
      rownames(ref) <- v
      true <- !is.na(ref[, 4])
      
# Time series models
load(file = paste0(path, "/sim_TS_0.RData")); TSM_0 <- sim
load(file = paste0(path, "/sim_TS_005.RData")); TSM_5 <- sim

# Spatial distribution models
load(file = paste0(path, "/sim_SDM_0.RData")); SDM_0 <- sim
load(file = paste0(path, "/sim_SDM_005.RData")); SDM_5 <- sim

# Process-based models
load(file = paste0(path, "/sim_POM_0.RData")); PBM_0 <- sim
load(file = paste0(path, "/sim_POM_005.RData")); PBM_5 <- sim

rm(sim)


sims_0 <- list(TSM_0, SDM_0, PBM_0)
sims_5 <- list(TSM_5, SDM_5, PBM_5)
s0 <- lapply(sims_0, \(x) evaluate_model(true = true, sim = x))
s5 <- lapply(sims_5, \(x) evaluate_model(true = true, sim = x))

names(s0) <- names(s5) <- lab

  
  if(length(m0) == 0){
    m0 <- list(s0)
    m5 <- list(s5)
  } else {
    m0 <- c(m0, list(s0))
    m5 <- c(m5, list(s5))
  }
}

names(m0) <- names(m5) <- species
rm(s0,s5)

save(m0, file = "data/model0.RData")
save(m5, file = "data/model5.RData")

# Drawing plots (Figure S4 and S5 Supplementary Material ) -----------------------------

load("data/model0.RData")
load("data/model5.RData")

# No selection 
op <- par(mfrow = c(1, 3), mar = c(5, 1, 1, 1), oma = c(0, 10, 0, 0))
plt_multi(m0, "ba", dy = 0.12, bw = 5, yaxis = T)
plt_multi(m0, "tpr", dy = 0.12, bw = 5)
plt_multi(m0, "tnr", dy = 0.12, bw = 5)


# Variable selection 
op <- par(mfrow = c(1, 3), mar = c(5, 1, 1, 1), oma = c(0, 10, 0, 0))
plt_multi(m5, "ba", dy = 0.12, bw = 5, yaxis = T)
plt_multi(m5, "tpr", dy = 0.12, bw = 5)
plt_multi(m5, "tnr", dy = 0.12, bw = 5)

