
# Fitting the spatiotemporal GAMM with an interacting species


# Setup -------------------------------------------------------------------

library(mgcv)

source("R/setup.R")
source("R/step.R")
source("R/set_species.R")


# The data ----------------------------------------------------------------

# Print species IDs
species_info(print_inter_species = TRUE)

# Interacting species
path <- paste0("data/", inter_species)
load(paste0(path, "/data.RData"))
inter <- data[c("idy", "dens", "log_dens", "log_dens_lag")]
names(inter)[-1] <- paste(inter_species_name, names(inter)[-1], sep = "_")

rm(data)

# Focal species
path <- paste0("data/", species)
load(paste0(path, "/data.RData"))

# Merging
data <- merge(data, inter, by = "idy")
summary(data)
rm(inter)

# List of common variables
# dput(names(data))
vc <- c("idy", "id_year", "id", "plot_id", "observer_id", "x", "y", "year",
        "fyear", "dens", "log_dens")

# List of predictors
v <- c("Natural_vegetation", "Deciduous_forest", "Coniferous_forest",
       "Mixed_forest", "Shrub_herbaceous", "Grassland", "Urban", "Water",
       "Tmin_spring_lag", "Tmin_spring", "Tmin_autumn_lag", "Tmin_winter",
       "Tmax_summer_lag", "Tmax_autumn_lag", "Precip_spring_lag", "Precip_spring",
       "Precip_summer_lag", "Precip_autumn_lag", "Precip_winter",
       paste0(inter_species_name, "_log_dens"), "log_dens_lag")


# GAMM formula ------------------------------------------------------------

f <- formula(paste0("dens ~ ", paste0("s(", v, ", k = k, bs = 'cr')", collapse = " + ")))

(f01 <- update(f, ". ~ . + s(plot_id, bs = 're') + s(observer_id, bs = 're') + s(fyear, bs = 're') + s(year, k = 10, bs = 'gp', m = c(3, 0.5)) + s(x, y, k = 30, bs = 'gp')"))


# GAMM fitting ------------------------------------------------------------

k <- 6
gamma <- 4
nthreads <- 8L

# Full model
full <- bam(f01, data, family = tw(), discrete = TRUE, nthreads = nthreads, gamma = gamma)
summary(full, re.test = FALSE)
pres_plot(full, scale = 0)

# Backward elimination
fit <- backward(full, nvar = 11)

# The final model
summary(fit)

# Some diagnostics:
op <- par(mfrow = c(2, 2)); gam.check(fit, rep = 100); par(op)
data$res <- residuals(fit)

# Partial response curves
pres_plot(fit, cex.lab = 1.4)
pres_plot(fit, cex.lab = 1.4, scale = 0)


# Saving models -----------------------------------------------------------

save(full, file = paste0(path, "/gamm_full_inter.RData"))
save(fit, file = paste0(path, "/gamm_selected_inter.RData"))
save(data, file = paste0(path, "/data_inter.RData"))
