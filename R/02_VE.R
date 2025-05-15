
# Generating replicated VE data


# Setup -------------------------------------------------------------------

library(mgcv)
library(data.table)

source("R/setup.R")
source("R/sim.R")
source("R/set_species.R")


# Load data and models -----------------------------------------------------
species_info()

path <- paste0("data/", species)

# The selected model
load(paste0(path, "/gamm_selected_inter.RData"))
summary(fit, re.test = FALSE)

# The real data
load(paste0(path, "/data_inter.RData"))
observed <- data$dens


# Prepare data for simulations --------------------------------------------

data[["eta"]] <- predict(fit, type = "link", exclude = c("s(plot_id)", "s(observer_id"))
pred_data <- data.table(data[c("id_year", "plot_id", "observer_id", "year", "eta")], key = "id_year")


# Simulate from the model ------------------------------------------------
# First run, to fit calibration functions.

vse <- rep_VSE(n = 1000, fit, pred_data, type = "none")

ve <- vse[, 2, ]


# Fitting Q-Q maps --------------------------------------------------------

# Mapping simulated VE data using observed data as a target:
h <- qq_map(ve, observed, m = 1)
summary(h)
xq <- seq(0, max(ve), len = 100)
yq <- predict(h, list(x = xq)); summary(yq)
plot(xq, yq, type = "l", xlab = "VE data", ylab = "Observed data"); grid(); abline(a = 0, b = 1, lty = 2)

# Calibration function
f = function(x) {
  hh <- predict(h, list(x = x))
  q <- quantile(hh, 1 - prevalence(observed))
  threshold(hh, q)
}

ve_prim <- ve
ve_prim[] <- pbapply(ve, 2, f)
dim(ve)


# Comparing simulated VE and observed data
ma <- max(max(log1p(observed)), max(log1p(ve_prim))) + 0.3
breaks <- seq(0, ma, 1/3)
op <- par(mfrow = c(2, 2), mar = c(4.5, 4.5, 2, 2))
hist_comp(log1p(ve_prim), log1p(observed), breaks = breaks, xname = "Population log-density")
rep_qqplot(observed, ve_prim) # Q-Q plot
rep_hist(observed, ve_prim, mean, breaks = 20, main = "", xlab = "Mean")
rep_hist(observed, ve_prim, sd, breaks = 20, main = "", xlab = "Standard deviation")
par(op)


# Saving simulation results -----------------------------------------------

save(pred_data, file = paste0(path, "/pred_data.RData"))
save(ve, file = paste0(path, "/ve.RData"))
save(ve_prim, file = paste0(path, "/ve_prim.RData"))

