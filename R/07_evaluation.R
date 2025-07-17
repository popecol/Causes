
# Evaluation of model performance

# Setup -------------------------------------------------------------------
 
library(readr)
library(glmmTMB)

source("R/setup.R")
source("R/tests.R")
source("R/plot_meta.R")


# Combining data from multiple species and models -------------------------

combine <- function(models, species) {
  
  cat(species, "\n")
  out <- data.frame()
  
  path <- paste0("data/", species)
  load(paste0(path, "/gamm_full_inter.RData"))
  
  v <- get_names(full)
  v <- v[!v %in% "log_dens_lag"]
  v <- v[!v %in% "year"]
  
  # The selected model (model used for generating a VS)
  load(paste0(path, "/gamm_selected_inter.RData"))
  
  # Variables in the generating model
  ref <- summary(fit, re.test = FALSE)[["s.table"]]
  idx <- match(v, get_names(rownames(ref)))
  ref <- ref[idx, ]
  rownames(ref) <- v
  true <- !is.na(ref[, 4])
  
  for (i in 1:length(models)) {
    
    model <- models[i]
    load(paste0(path, "/sim_", model, ".RData"))
    model_eval <- evaluate_model(true, sim)
    model_eval$species <- species
    model_eval$model <- model
    out <- rbind(out, model_eval)
    rm(sim)
  }
  return(out)
}


models <- c("TS_0", "TS_005", "SDM_0", "SDM_005", "POM_0", "POM_005")
species_list <- c("KS", "FH", "SA", "SU", "SB", "SC", "Z", "KT", "PP", "E")

results_list <- lapply(species_list, function(x) combine(models, x))
names(results_list) <- species_list

measures <- do.call(rbind, results_list)
measures[1:4] <- measures[1:4] / 100
measures <- transform(measures, species = factor(species), model = factor(model, levels = models))
summary(measures)

save(measures, file = "data/measures.RData")


# Fitting a GLMM ----------------------------------------------------------------

library(emmeans)

load("data/measures.RData")

m <- strsplit(as.character(measures$model), "_", fixed = TRUE)
measures$mod <- factor(sapply(m, \(x) x[[1]]), levels = c("TS", "SDM", "POM"))
sel <- factor(sapply(m, \(x) x[[2]]))
levels(sel) <- c("withouth_var_selection", "with_var_selection")
measures$sel <- sel
summary(measures)

measures$w <- rep(20, nrow(measures))


# ACC
hist(measures$acc)
plot(acc ~ mod, measures, las = 2)
plot(acc ~ sel, measures, las = 2)
plot(acc ~ species, measures)

m_acc <- glmmTMB(acc ~ mod * sel + (1 | species), data = measures, family = binomial, weights = w)
summary(m_acc)
car::Anova(m_acc)glmmTMB:::Anova.II.glmmTMB(m_acc)

emmeans(m_acc, pairwise ~ mod | sel, type = "response")
emmeans(m_acc, pairwise ~ sel | mod, type = "response")
emmip(m_acc, ~ sel | mod, type = "response", CIs = TRUE)
emmip(m_acc, ~ mod | sel, type = "response", CIs = TRUE)


# TPR
hist(measures$tpr)
plot(tpr ~ model, measures, las = 2)
plot(tpr ~ species, measures)

m_tpr <- glmmTMB(tpr ~ mod * sel + (1 | species), data = measures, family = binomial, weights = w)
summary(m_tpr)
car::Anova(m_tpr)

plot(tpr ~ sel, measures)

emmeans(m_tpr, pairwise ~ mod | sel, type = "response")
emmeans(m_tpr, pairwise ~ sel | mod, type = "response")
emmip(m_tpr, ~ sel | mod, type = "response", CIs = TRUE)
emmip(m_tpr, ~ mod | sel, type = "response", CIs = TRUE)


# TNR
hist(measures$tnr)
plot(tnr ~ model, measures, las = 2)
plot(tnr ~ species, measures)

m_tnr <- glmmTMB(tnr ~ mod * sel + (1 | species), data = measures, family = binomial, weights = w)
summary(m_tnr)
car::Anova(m_tnr)

emmeans(m_tnr, pairwise ~ mod | sel, type = "response")
emmeans(m_tnr, pairwise ~ sel | mod, type = "response")
emmip(m_tnr, ~ sel | mod, type = "response", CIs = TRUE)
emmip(m_tnr, ~ mod | sel, type = "response", CIs = TRUE)



# Graph (Figure 2, main article)  ------------------------------------------------------------

op <- par(mfrow = c(1, 3), mar = c(6.5, 5, 2, 2))
plot_meta(m_acc, ran = F, ylab = "Accuracy")
plot_meta(m_tpr, ran = F, ylab = "True positive rate")
plot_meta(m_tnr, ran = F, ylab = "True negative rate")
par(op)


