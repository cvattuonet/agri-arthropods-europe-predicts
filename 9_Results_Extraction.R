#---------------------------------------------------------------------------------------------------------------------------------------#
# CODE DETAILS                                                                                                                          #
#                                                                                                                                       #
# Author: Catalina Vattuone  (cvattuonet@gmail.com)                                                                                     #        
# Date of latest update: 19-08-2026                                                                                                     #
# Main use: Study on "Assessing how agricultural practices and landscape composition shape differences in arthropod biodiversity among  #  
#           habitat types in Europe using the PREDICTS database". Internship at INRAE                                                   #
#                                                                                                                                       #
# Content: Extraction code to obtain excel fill with results from all models. It includes:                                              # 
#         - Model information (family, formula, AIC, BIC, logLik, df.resid)                                                             #
#         - Raw estimates for all models and percentage of change compare to reference for average abundance and richness models        #
#         - DHARMa and check_model plots                                                                                                #
#         - Values are presented as one sheet per model. With sheets for all models with delta AIC < 2 for the gaussian abundance model # 
#           and for the richness model. For the rest of the abundance models results are presented for the best-AIC model               #
#---------------------------------------------------------------------------------------------------------------------------------------#

#LIBRARIES  ------------------------------------------------------------
pacman::p_load(readr, dplyr, stringr, purrr, openxlsx, MuMIn, lme4, lmerTest, glmmTMB, broom.mixed, performance, DHARMa)

#LOAD DATABASES   ------------------------------------------------------------
sites_df <- read_delim("Intermediate_dataset/sites_for_abundance_models.csv")
sites_df$Habitat_category <- relevel(factor(sites_df$Habitat_category), ref = "Natural vegetation")

sites_df_richness <- read_delim("Intermediate_dataset/sites_for_richness_models.csv")
sites_df_richness$Habitat_category <- relevel(factor(sites_df_richness$Habitat_category), ref = "Natural vegetation")

#LOAD MODELS------------------------------------------------------------
#abundance
global_gaussian <- readRDS("Models_results/global_gaussian.rds")
dredge_gaussian <- readRDS("Models_results/dredge_gaussian.rds")

global_poisson <- readRDS("Models_results/global_poisson.rds")
dredge_poisson <- readRDS("Models_results/dredge_poisson.rds")

global_nb <- readRDS("Models_results/global_nb.rds")
dredge_nb <- readRDS("Models_results/dredge_nb.rds")

global_tweedie <- readRDS("Models_results/global_tweedie.rds")
dredge_tweedie <- readRDS("Models_results/dredge_tweedie.rds")

global_zinb <- readRDS("Models_results/global_zinb.rds")
dredge_zinb <- readRDS("Models_results/dredge_zinb.rds")

#richness
global_zinb_richness<- readRDS("Models_results/global_zinb_richness.rds")
dredge_zinb_richness <- readRDS("Models_results/dredge_zinb_richness.rds")

##LASSO MODELS ------------------------------------------------------------
lasso_ab_optimal  <- readRDS("Models_results/lasso_ab_optimal.rds")
lasso_ab_500      <- readRDS("Models_results/lasso_ab_500.rds")
lasso_rch_optimal <- readRDS("Models_results/lasso_rch_optimal.rds")

#DIRECTORY------------------------------------------------------------
out_dir <- "Models_results/Final/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

#FUNCTIONS ------------------------------------------------------------
#pulls family/formula/fit-statistics, works for both lmer and glmmTMB objects
extract_model_info <- function(model) {
  fam <- tryCatch(family(model)$family, error = function(e) "gaussian (lmer)")
  list(
    family     = fam,
    formula    = paste(deparse(formula(model), width.cutoff = 500), collapse = " "),
    AIC        = AIC(model),
    BIC        = BIC(model),
    logLik     = as.numeric(logLik(model)),
    neg2logLik = -2 * as.numeric(logLik(model)),
    df_resid   = tryCatch(df.residual(model), error = function(e) NA)
  )
}

#random effects variance/SD table - lmer and glmmTMB expose this differently
get_varcorr_df <- function(model) {
  if (inherits(model, "merMod")) {
    as.data.frame(VarCorr(model)) %>%
      transmute(Group = grp, Term = ifelse(is.na(var1), "(Intercept)", var1),
                Variance = vcov, Std.Dev = sdcor)
  } else {
    broom.mixed::tidy(model, effects = "ran_pars", scales = "vcov") %>%
      transmute(Group = group, Term = term, Variance = estimate, Std.Dev = sqrt(estimate))
  }
}

#number of observations + number of groups per random effect
get_group_counts <- function(model) {
  ngrps <- if (inherits(model, "merMod")) summary(model)$ngrps else summary(model)$ngrps$cond
  bind_rows(
    data.frame(Group = "Observations", N = nobs(model)),
    data.frame(Group = names(ngrps), N = as.numeric(ngrps))
  )
}

#fixed effects estimate/SE - lmer and glmmTMB store this in different places
get_fixed_effects <- function(model) {
  if (inherits(model, "merMod")) {
    coefs <- as.data.frame(summary(lmerTest::as_lmerModLmerTest(model))$coefficients)
    coefs <- coefs %>% rename(p_value = `Pr(>|t|)`) %>% select(Estimate, `Std. Error`, p_value)
  } else {
    coefs <- as.data.frame(summary(model)$coefficients$cond)
    coefs <- coefs %>% rename(p_value = `Pr(>|z|)`) %>% select(Estimate, `Std. Error`, p_value)
  }
  coefs$Term <- rownames(coefs)
  coefs %>% select(Term, Estimate, `Std. Error`, p_value) %>% `rownames<-`(NULL)
}

#builds one full sheet for one model
build_model_sheet <- function(wb, sheet_name, model, n_sim = 12500) {
  addWorksheet(wb, sheet_name)
  setColWidths(wb, sheet_name, cols = 1:4, widths = c(48, 12, 12, 12))
  cat("Building sheet:", sheet_name, "\n")
  
  info   <- extract_model_info(model)
  vc_df  <- get_varcorr_df(model)
  grp_df <- get_group_counts(model)
  fe_df  <- get_fixed_effects(model)
  r2_val <- tryCatch(performance::r2(model), error = function(e) list(R2_conditional = NA, R2_marginal = NA))
  
  r <- 1
  writeData(wb, sheet_name, "Family", startRow = r, startCol = 1); writeData(wb, sheet_name, info$family, startRow = r, startCol = 2)
  r <- r + 1
  writeData(wb, sheet_name, "Formula", startRow = r, startCol = 1); writeData(wb, sheet_name, info$formula, startRow = r, startCol = 2)
  r <- r + 2
  
  writeData(wb, sheet_name, data.frame(
    Metric = c("AIC", "BIC", "logLik", "-2*logLik", "df.resid"),
    Value  = c(info$AIC, info$BIC, info$logLik, info$neg2logLik, info$df_resid)
  ), startRow = r, startCol = 1)
  r <- r + 8
  
  writeData(wb, sheet_name, "Random effects (variance / SD)", startRow = r, startCol = 1)
  r <- r + 1
  writeData(wb, sheet_name, vc_df, startRow = r, startCol = 1)
  r <- r + nrow(vc_df) + 3
  
  writeData(wb, sheet_name, "N observations / groups", startRow = r, startCol = 1)
  r <- r + 1
  writeData(wb, sheet_name, grp_df, startRow = r, startCol = 1)
  r <- r + nrow(grp_df) + 3
  
  writeData(wb, sheet_name, "Fixed effects", startRow = r, startCol = 1)
  r <- r + 1
  writeData(wb, sheet_name, fe_df, startRow = r, startCol = 1)
  r <- r + nrow(fe_df) + 3
  
  writeData(wb, sheet_name, data.frame(
    Metric = c("R2 conditional", "R2 marginal"),
    Value  = c(r2_val$R2_conditional, r2_val$R2_marginal)
  ), startRow = r, startCol = 1)
  r <- r + 6
  
  #Cook's distance - only supported natively for lmer (merMod) objects
  if (inherits(model, "merMod")) {
    cooks_max <- max(cooks.distance(model), na.rm = TRUE)
    writeData(wb, sheet_name, data.frame(
      Metric = "Max Cook's distance", Value = cooks_max
    ), startRow = r, startCol = 1)
    r <- r + 3
  }
  
  #DHARMa plot
  set.seed(123)
  sim <- DHARMa::simulateResiduals(model, n = n_sim)
  tmp_dharma <- tempfile(fileext = ".png")
  png(tmp_dharma, width = 1600, height = 800, res = 150)
  plot(sim)
  dev.off()
  writeData(wb, sheet_name, paste0("DHARMa residual diagnostics (n = ", n_sim, " simulations)"), startRow = r, startCol = 1)
  insertImage(wb, sheet_name, tmp_dharma, startRow = r + 1, startCol = 1, width = 9, height = 4.5)
  r <- r + 24
  
  #check_model plot
  cm_plot <- tryCatch(performance::check_model(model), error = function(e) NULL)
  if (!is.null(cm_plot)) {
    tmp_cm <- tempfile(fileext = ".png")
    png(tmp_cm, width = 2000, height = 2800, res = 150)
    print(cm_plot)
    dev.off()
    writeData(wb, sheet_name, "check_model() diagnostics", startRow = r, startCol = 1)
    insertImage(wb, sheet_name, tmp_cm, startRow = r + 1, startCol = 1, width = 10, height = 14)
  }
}

build_avg_model_sheet <- function(wb, sheet_name, dredge_obj, models_list, habitat_df, landscape_df, interactions_df, delta_cutoff = 2) {
  addWorksheet(wb, sheet_name)
  setColWidths(wb, sheet_name, cols = 1:8, widths = c(48, 12, 12, 12, 12, 12, 12, 12))
  cat("Building sheet:", sheet_name, "\n")
  
  n_models  <- length(models_list)
  best_aicc <- dredge_obj$AICc[1]
  
  r2_vals <- purrr::map_dfr(models_list, function(m) {
    r2 <- tryCatch(performance::r2(m), error = function(e) list(R2_conditional = NA, R2_marginal = NA))
    data.frame(R2_conditional = as.numeric(r2$R2_conditional), R2_marginal = as.numeric(r2$R2_marginal))
  })
  
  r <- 1
  writeData(wb, sheet_name, "Model-averaged results", startRow = r, startCol = 1)
  r <- r + 1
  writeData(wb, sheet_name, paste0("Averaged over ", n_models, " models with delta AICc < ", delta_cutoff), startRow = r, startCol = 1)
  r <- r + 2
  
  writeData(wb, sheet_name, data.frame(
    Metric = c("Best model AICc", "N models averaged", "R2 conditional (mean)", "R2 conditional (SD)", "R2 marginal (mean)", "R2 marginal (SD)"),
    Value  = c(best_aicc, n_models,
               mean(r2_vals$R2_conditional, na.rm = TRUE), sd(r2_vals$R2_conditional, na.rm = TRUE),
               mean(r2_vals$R2_marginal, na.rm = TRUE), sd(r2_vals$R2_marginal, na.rm = TRUE))
  ), startRow = r, startCol = 1)
  r <- r + 9
  
  writeData(wb, sheet_name, "Habitat category effects (% change vs. Natural vegetation)", startRow = r, startCol = 1)
  r <- r + 1
  writeData(wb, sheet_name, habitat_df, startRow = r, startCol = 1)
  r <- r + nrow(habitat_df) + 3
  
  writeData(wb, sheet_name, "Landscape composition / practice main effects (% change)", startRow = r, startCol = 1)
  r <- r + 1
  writeData(wb, sheet_name, landscape_df, startRow = r, startCol = 1)
  r <- r + nrow(landscape_df) + 3
  
  writeData(wb, sheet_name, "Simple slopes by habitat (variables whose effect depends on habitat)", startRow = r, startCol = 1)
  r <- r + 1
  writeData(wb, sheet_name, interactions_df, startRow = r, startCol = 1)
  r <- r + nrow(interactions_df) + 3
}

#GET BEST  MODELS ------------------------------------------------------------
gaussian_models <- get.models(dredge_gaussian, subset = delta < 2)
cat("Gaussian delta<2 models found:", length(gaussian_models), "\n")   #confirm this is 10
gaussian_models <- lapply(gaussian_models, function(m) update(m, REML = TRUE))

poisson_best <- get.models(dredge_poisson, subset = 1)[[1]]
negbin_best  <- get.models(dredge_nb,  subset = 1)[[1]]
tweedie_best <- get.models(dredge_tweedie,  subset = 1)[[1]]
zinb_best    <- get.models(dredge_zinb,  subset = 1)[[1]]

richness_models <- get.models(dredge_zinb_richness, subset = delta < 2)
cat("Richness (ZINB) delta<2 models found:", length(richness_models), "\n")

#CREATE SHEET NAMES AND DESCRIPTIONS ------------------------------------------------------------
#for abundance
gaussian_sheet_names <- paste0("AB_G_", seq_along(gaussian_models))
gaussian_descriptions <- paste0("Abundance - Gaussian model, rank ", seq_along(gaussian_models),
                                " of delta AIC < 2 (", length(gaussian_models), " models total)")

sheet_index <- data.frame(
  Sheet = c(gaussian_sheet_names, "AB_P", "AB_NB", "AB_T", "AB_ZINB"),
  Description = c(gaussian_descriptions,
                  "Abundance - Poisson model, best model (lowest AICc)",
                  "Abundance - Negative Binomial model, best model (lowest AICc)",
                  "Abundance - Tweedie model, best model (lowest AICc)",
                  "Abundance - Zero-inflated Negative Binomial model, best model (lowest AICc)"))

#for richness
richness_sheet_names <- paste0("R_ZINB_", seq_along(richness_models))
richness_descriptions <- paste0("Taxa richness - Zero-inflated Negative Binomial model, rank ", seq_along(richness_models),
                                " of delta AIC < 2 (", length(richness_models), " models total)")
sheet_index <- bind_rows(sheet_index, data.frame(Sheet = richness_sheet_names, Description = richness_descriptions))


#for average models
sheet_index <- bind_rows(sheet_index, data.frame(Sheet = c("AB_AVG", "R_AVG"),
                                                 Description = c("Abundance - model-averaged results (habitat, landscape, and interaction effects) across delta AICc < 2 model set",
                                                                 "Taxa richness - model-averaged results (habitat, landscape, and interaction effects) across delta AICc < 2 model set")))


#BUILD WORKBOOK ------------------------------------------------------------
wb <- createWorkbook()

addWorksheet(wb, "Index")
writeData(wb, "Index", "Abundance model results - sheet index", startRow = 1, startCol = 1)
writeData(wb, "Index", sheet_index, startRow = 3, startCol = 1)

for (i in seq_along(gaussian_models)) {
  build_model_sheet(wb, gaussian_sheet_names[i], gaussian_models[[i]])
}

build_model_sheet(wb, "AB_P",    poisson_best)
build_model_sheet(wb, "AB_NB",   negbin_best)
build_model_sheet(wb, "AB_T",    tweedie_best)
build_model_sheet(wb, "AB_ZINB", zinb_best)

saveWorkbook(wb, paste0(out_dir, "Model_Results.xlsx"), overwrite = TRUE)

#ADDING RICHNESS MODELS ------------------------------------------------------------
writeData(wb, "Index", sheet_index, startRow = 3, startCol = 1)

for (i in seq_along(richness_models)) {
  build_model_sheet(wb, richness_sheet_names[i], richness_models[[i]])
}

saveWorkbook(wb, paste0(out_dir, "Model_Results.xlsx"), overwrite = TRUE)

##MAIN RESULTS ------------------------------------------------------------
get_coefs_vcov_avg <- function(dredge_obj, delta_cutoff = 2) {
  n_models <- nrow(subset(dredge_obj, delta < delta_cutoff))
  cat("Averaging", n_models, "models with delta <", delta_cutoff, "\n")
  
  avg <- model.avg(dredge_obj, subset = delta < delta_cutoff, fit = TRUE)
  b  <- coef(avg, full = TRUE)
  vc <- as.matrix(vcov(avg, full = TRUE))
  
  clean <- function(x) sub("^cond\\((.*)\\)$", "\\1", x)
  names(b) <- clean(names(b))
  dimnames(vc) <- list(clean(rownames(vc)), clean(colnames(vc)))
  
  keep <- !grepl("^zi\\(|^zi~", names(b))
  list(b = b[keep], vc = vc[keep, keep, drop = FALSE])
}

cv_m12  <- get_coefs_vcov_avg(dredge_gaussian)
cv_m15b <- get_coefs_vcov_avg(dredge_zinb_richness)

land_use_levels <- c("Natural vegetation", "Semi-natural vegetation", "Pasture", "Cropland")

#split a Habitat_category term/level string into land_use + intensity
parse_habitat_level <- function(x) {
  data.frame(intensity = factor(case_when(
    str_detect(x, "Minimal") ~ "Minimal",
    str_detect(x, "Light")   ~ "Light",
    str_detect(x, "Intense") ~ "Intense",
    TRUE ~ "Reference"), levels = c("Reference", "Minimal", "Light", "Intense")),
    land_use = str_remove(x, "_(Minimal|Light|Intense) use$") )
}

count_habitat_sites <- function(data) {
  counts <- data %>% count(Habitat_category, name = "n")
  bind_cols(counts, parse_habitat_level(as.character(counts$Habitat_category))) %>%
    select(land_use, intensity, n)
}

get_pvalues <- function(dredge_obj, delta_cutoff = 2) {
  avg <- model.avg(dredge_obj, subset = delta < delta_cutoff, fit = TRUE)
  pmat <- as.data.frame(summary(avg)$coefmat.full)
  pmat$term <- sub("^cond\\((.*)\\)$", "\\1", rownames(pmat))
  pmat %>% select(term, p_value = `Pr(>|z|)`)
}

get_habitat_pctdiff <- function(cv) {
  terms <- grep("^Habitat_category[^:]+$", names(cv$b), value = TRUE)
  df <- data.frame(term = terms, estimate = as.numeric(cv$b[terms]), se = sqrt(diag(cv$vc)[terms])) %>%
    mutate(level = str_remove(term, "^Habitat_category")) %>%
    bind_cols(parse_habitat_level(.$level)) %>%
    mutate(conf.low = estimate - 1.96 * se,
           conf.high = estimate + 1.96 * se,
           pct_diff = (exp(estimate) - 1) * 100,
           pct_low  = (exp(conf.low) - 1) * 100,
           pct_high = (exp(conf.high) - 1) * 100)
  
  bind_rows(data.frame(term = NA_character_, land_use = "Natural vegetation", intensity = "Reference", pct_diff = 0, pct_low = 0, pct_high = 0),
            df %>% select(term, land_use, intensity, pct_diff, pct_low, pct_high) ) %>%
    mutate(land_use = factor(land_use, levels = land_use_levels))
}

pvals_m12  <- get_pvalues(dredge_gaussian)
pvals_m15b <- get_pvalues(dredge_zinb_richness)

habitat_m12  <- get_habitat_pctdiff(cv_m12)  %>%
  left_join(count_habitat_sites(sites_df), by = c("land_use", "intensity")) %>%
  mutate(land_use = factor(land_use, levels = land_use_levels))

habitat_m15b <- get_habitat_pctdiff(cv_m15b) %>%
  left_join(count_habitat_sites(sites_df_richness), by = c("land_use", "intensity")) %>%
  mutate(land_use = factor(land_use, levels = land_use_levels))

habitat_m12  <- habitat_m12  %>% left_join(pvals_m12,  by = "term")
habitat_m15b <- habitat_m15b %>% left_join(pvals_m15b, by = "term")

## LANDSCAPE MAIN EFFECTS - join p-values ------------------------------------------------------------
predictor_labels <- c(
  urban_pct_1000m_mod                  = "Urban cover",
  water_pct_1000m_mod                  = "Water cover",
  natural_pct_1000m_mod                = "Natural habitat cover",
  pasture_pct_1000m_mod                = "Pasture cover",
  pesticide_1000m                      = "Pesticide use",
  crop_diversity_machefer_alpha_1000m  = "Crop diversity",
  pct_very_small_and_small_field_1000m = "Small field prevalence",
  pct_very_large_and_large_field_1000m = "Large field prevalence",
  climate_mean_temp_c_1000m            = "Mean temperature")

get_landscape_pctdiff <- function(cv) {
  terms <- setdiff(names(cv$b), c("(Intercept)", "(Int)", grep("Habitat_category", names(cv$b), value = TRUE)))
  terms <- terms[!str_detect(terms, ":")]
  
  data.frame(term = terms, estimate = as.numeric(cv$b[terms]), se = sqrt(diag(cv$vc)[terms])) %>%
    mutate(conf.low = estimate - 1.96 * se, conf.high = estimate + 1.96 * se,
           pct_diff = (exp(estimate) - 1) * 100,
           pct_low  = (exp(conf.low) - 1) * 100,
           pct_high = (exp(conf.high) - 1) * 100,
           label = recode(term, !!!predictor_labels),
           depends_on_habitat = purrr::map_lgl(term, ~ any(grepl(
             paste0("(^Habitat_category.*:", .x, "$)|(^", .x, ":Habitat_category.*$)"), names(cv$b)))))
}

landscape_m12  <- get_landscape_pctdiff(cv_m12)
landscape_m15b <- get_landscape_pctdiff(cv_m15b)

landscape_m12  <- landscape_m12  %>% left_join(pvals_m12,  by = "term")
landscape_m15b <- landscape_m15b %>% left_join(pvals_m15b, by = "term")

##SIMPLE SLOPES BY HABITAT - recompute with p-values ------------------------------------------------------------
compute_simple_slopes <- function(cv, predictor) {
  b <- cv$b; vc <- cv$vc
  if (!(predictor %in% names(b))) return(NULL)
  int_terms <- grep(paste0("(^Habitat_category.*:", predictor, "$)|(^", predictor, ":Habitat_category.*$)"), names(b), value = TRUE)
  if (length(int_terms) == 0) return(NULL)
  levels_found <- str_remove(str_extract(int_terms, "Habitat_category[^:]+"), "^Habitat_category")
  ref_row <- data.frame(level = "Natural vegetation", slope = as.numeric(b[[predictor]]), se = sqrt(vc[predictor, predictor]))
  
  other_rows <- purrr::map2_dfr(int_terms, levels_found, function(term, lvl) {
    slope_val <- as.numeric(b[[predictor]] + b[[term]])
    var_val   <- vc[predictor, predictor] + vc[term, term] + 2 * vc[predictor, term]
    data.frame(level = lvl, slope = slope_val, se = sqrt(var_val))
  })
  
  bind_rows(ref_row, other_rows) %>%
    mutate(predictor = predictor,
           z = slope / se,
           p_value = 2 * (1 - pnorm(abs(z))),
           conf.low = slope - 1.96 * se, conf.high = slope + 1.96 * se,
           pct_diff = (exp(slope) - 1) * 100,
           pct_low  = (exp(conf.low) - 1) * 100,
           pct_high = (exp(conf.high) - 1) * 100) %>%
    bind_cols(parse_habitat_level(.$level))
}

#rerun so interactions_m12/m15b pick up the new p_value column
interacting_predictors <- c("pesticide_1000m", "crop_diversity_machefer_alpha_1000m", "pct_very_small_and_small_field_1000m", "natural_pct_1000m_mod", "water_pct_1000m_mod")

get_interaction_slopes <- function(cv) {
  df <- purrr::map_dfr(interacting_predictors, ~compute_simple_slopes(cv, .x))
  present <- interacting_predictors[interacting_predictors %in% unique(df$predictor)]
  df %>%
    mutate(land_use = factor(land_use, levels = land_use_levels),
           label = factor(recode(predictor, !!!predictor_labels), levels = unname(predictor_labels[present])))
}

interactions_m12  <- get_interaction_slopes(cv_m12)
interactions_m15b <- get_interaction_slopes(cv_m15b)

#ADDING AVERAGED MODEL SHEETS ------------------------------------------------------------
build_avg_model_sheet( wb, "AB_AVG", dredge_gaussian, gaussian_models,
                       habitat_m12     %>% select(term, land_use, intensity, pct_diff, pct_low, pct_high, p_value, n),
                       landscape_m12   %>% select(label, pct_diff, pct_low, pct_high, p_value, depends_on_habitat, legend_group),
                       interactions_m12 %>% select(label, land_use, intensity, pct_diff, pct_low, pct_high, p_value))

build_avg_model_sheet(wb, "R_AVG", dredge_zinb_richness, richness_models,
                      habitat_m15b     %>% select(term, land_use, intensity, pct_diff, pct_low, pct_high, p_value, n),
                      landscape_m15b   %>% select(label, pct_diff, pct_low, pct_high, p_value, depends_on_habitat, legend_group),
                      interactions_m15b %>% select(label, land_use, intensity, pct_diff, pct_low, pct_high, p_value))

writeData(wb, "Index", sheet_index, startRow = 3, startCol = 1)
saveWorkbook(wb, paste0(out_dir, "Model_Results.xlsx"), overwrite = TRUE)



#ADDING LASSO RESULTS ------------------------------------------------------------
##LASSO HELPERS ------------------------------------------------------------
get_lasso_se <- function(model) {
  se <- as.numeric(model$fixerror)
  names(se) <- names(model$coefficients)  #fixerror follows the same order as coefficients
  se
}

get_lasso_habitat_pctdiff <- function(model, data) {
  b  <- model$coefficients
  se <- get_lasso_se(model)
  terms <- grep("^Habitat_category[^:]+$", names(b), value = TRUE)
  
  df <- data.frame(term = terms, estimate = as.numeric(b[terms]), se = as.numeric(se[terms])) %>%
    mutate(level = str_remove(term, "^Habitat_category")) %>%
    bind_cols(parse_habitat_level(.$level)) %>%
    mutate(dropped   = estimate == 0,   #shrunk to exactly 0 by the LASSO penalty
           conf.low  = estimate - 1.96 * se, conf.high = estimate + 1.96 * se,
           z         = estimate / se,
           p_value   = ifelse(dropped, NA_real_, 2 * (1 - pnorm(abs(z)))),
           pct_diff  = (exp(estimate) - 1) * 100,
           pct_low   = ifelse(dropped, NA_real_, (exp(conf.low) - 1) * 100),
           pct_high  = ifelse(dropped, NA_real_, (exp(conf.high) - 1) * 100))
  
  bind_rows(
    data.frame(term = NA_character_, land_use = "Natural vegetation", intensity = "Reference",
               pct_diff = 0, pct_low = NA_real_, pct_high = NA_real_, p_value = NA_real_, dropped = FALSE),
    df %>% select(term, land_use, intensity, pct_diff, pct_low, pct_high, p_value, dropped)
  ) %>%
    mutate(land_use = factor(land_use, levels = land_use_levels)) %>%
    left_join(count_habitat_sites(data), by = c("land_use", "intensity"))
}

get_lasso_landscape_pctdiff <- function(model) {
  b  <- model$coefficients
  se <- get_lasso_se(model)
  terms <- setdiff(names(b), c("(Intercept)", grep("Habitat_category", names(b), value = TRUE)))
  terms <- terms[!str_detect(terms, ":")]
  
  data.frame(term = terms, estimate = as.numeric(b[terms]), se = as.numeric(se[terms])) %>%
    mutate(dropped   = estimate == 0,
           conf.low  = estimate - 1.96 * se, conf.high = estimate + 1.96 * se,
           z         = estimate / se,
           p_value   = ifelse(dropped, NA_real_, 2 * (1 - pnorm(abs(z)))),
           pct_diff  = (exp(estimate) - 1) * 100,
           pct_low   = ifelse(dropped, NA_real_, (exp(conf.low) - 1) * 100),
           pct_high  = ifelse(dropped, NA_real_, (exp(conf.high) - 1) * 100),
           label = recode(term, !!!predictor_labels),
           depends_on_habitat = purrr::map_lgl(term, ~ any(grepl(
             paste0("(^Habitat_category.*:", .x, "$)|(^", .x, ":Habitat_category.*$)"), names(b)))))
}

compute_lasso_simple_slopes <- function(model, predictor) {
  b  <- model$coefficients
  se <- get_lasso_se(model)
  if (!(predictor %in% names(b))) return(NULL)
  int_terms <- grep(paste0("(^Habitat_category.*:", predictor, "$)|(^", predictor, ":Habitat_category.*$)"), names(b), value = TRUE)
  if (length(int_terms) == 0) return(NULL)
  levels_found <- str_remove(str_extract(int_terms, "Habitat_category[^:]+"), "^Habitat_category")
  
  ref_row <- data.frame(level = "Natural vegetation", slope = as.numeric(b[[predictor]]), se = as.numeric(se[[predictor]]))
  
  #approximation: no covariance term available from glmmLasso, so Var(sum) = Var(predictor) + Var(term) only
  other_rows <- purrr::map2_dfr(int_terms, levels_found, function(term, lvl) {
    slope_val <- as.numeric(b[[predictor]] + b[[term]])
    se_val    <- sqrt(se[[predictor]]^2 + se[[term]]^2)
    data.frame(level = lvl, slope = slope_val, se = se_val)
  })
  
  bind_rows(ref_row, other_rows) %>%
    mutate(predictor = predictor, z = slope / se, p_value = 2 * (1 - pnorm(abs(z))),
           conf.low = slope - 1.96 * se, conf.high = slope + 1.96 * se,
           pct_diff = (exp(slope) - 1) * 100,
           pct_low  = (exp(conf.low) - 1) * 100,
           pct_high = (exp(conf.high) - 1) * 100) %>%
    bind_cols(parse_habitat_level(.$level))
}

get_lasso_interaction_slopes <- function(model) {
  df <- purrr::map_dfr(interacting_predictors, ~compute_lasso_simple_slopes(model, .x))
  if (nrow(df) == 0) return(df)
  present <- interacting_predictors[interacting_predictors %in% unique(df$predictor)]
  df %>%
    mutate(land_use = factor(land_use, levels = land_use_levels),
           label = factor(recode(predictor, !!!predictor_labels), levels = unname(predictor_labels[present])))
}

build_lasso_sheet <- function(wb, sheet_name, model, lambda, habitat_df, landscape_df, interactions_df) {
  if (sheet_name %in% names(wb)) removeWorksheet(wb, sheet_name)
  addWorksheet(wb, sheet_name)
  setColWidths(wb, sheet_name, cols = 1:8, widths = c(48, 12, 12, 12, 12, 12, 12, 12))
  cat("Building sheet:", sheet_name, "\n")
  
  r <- 1
  writeData(wb, sheet_name, "LASSO-penalized model (glmmLasso)", startRow = r, startCol = 1); r <- r + 1
  writeData(wb, sheet_name, paste0("Lambda (penalty): ", lambda), startRow = r, startCol = 1); r <- r + 2
  
  writeData(wb, sheet_name, data.frame(Metric = c("AIC", "BIC"), Value = c(model$aic, model$bic)), startRow = r, startCol = 1)
  r <- r + 5
  
  writeData(wb, sheet_name, "Habitat category effects (% change vs. Natural vegetation)", startRow = r, startCol = 1); r <- r + 1
  writeData(wb, sheet_name, habitat_df, startRow = r, startCol = 1); r <- r + nrow(habitat_df) + 3
  
  writeData(wb, sheet_name, "Landscape composition / practice main effects (% change)", startRow = r, startCol = 1); r <- r + 1
  writeData(wb, sheet_name, landscape_df, startRow = r, startCol = 1); r <- r + nrow(landscape_df) + 3
  
  writeData(wb, sheet_name, "Simple slopes by habitat (variables whose effect depends on habitat)", startRow = r, startCol = 1); r <- r + 1
  writeData(wb, sheet_name, interactions_df, startRow = r, startCol = 1); r <- r + nrow(interactions_df) + 3
  
  writeData(wb, sheet_name,
            "Note: SEs/CIs/p-values use glmmLasso's per-coefficient SE only (no fixed-effect covariance matrix is available); interaction (simple-slope) SEs ignore the covariance between main effect and interaction term and are approximate. Coefficients shrunk to exactly 0 by the LASSO penalty are marked 'dropped'.",
            startRow = r, startCol = 1)
}

##LASSO RESULTS ------------------------------------------------------------
optimal_lambda_bic_ab  <- 1650
optimal_lambda_bic_rch <- 380

habitat_ab_opt  <- get_lasso_habitat_pctdiff(lasso_ab_optimal,  sites_df)
habitat_ab_500  <- get_lasso_habitat_pctdiff(lasso_ab_500,      sites_df)
habitat_rch_opt <- get_lasso_habitat_pctdiff(lasso_rch_optimal, sites_df_richness)

landscape_ab_opt  <- get_lasso_landscape_pctdiff(lasso_ab_optimal)
landscape_ab_500  <- get_lasso_landscape_pctdiff(lasso_ab_500)
landscape_rch_opt <- get_lasso_landscape_pctdiff(lasso_rch_optimal)

interactions_ab_opt  <- get_lasso_interaction_slopes(lasso_ab_optimal)
interactions_ab_500  <- get_lasso_interaction_slopes(lasso_ab_500)
interactions_rch_opt <- get_lasso_interaction_slopes(lasso_rch_optimal)


wb <- loadWorkbook(paste0(out_dir, "Model_Results.xlsx"))
sheet_index <- readWorkbook(paste0(out_dir, "Model_Results.xlsx"), sheet = "Index", startRow = 3)

build_lasso_sheet(wb, "AB_LASSO_OPT", lasso_ab_optimal, optimal_lambda_bic_ab,
                  habitat_ab_opt   %>% select(term, land_use, intensity, pct_diff, pct_low, pct_high, p_value, n, dropped),
                  landscape_ab_opt %>% select(label, pct_diff, pct_low, pct_high, p_value, depends_on_habitat, dropped),
                  interactions_ab_opt %>% select(label, land_use, intensity, pct_diff, pct_low, pct_high, p_value))

build_lasso_sheet(wb, "AB_LASSO_500", lasso_ab_500, 500,   # fix once lambda=500 is actually used above
                  habitat_ab_500   %>% select(term, land_use, intensity, pct_diff, pct_low, pct_high, p_value, n, dropped),
                  landscape_ab_500 %>% select(label, pct_diff, pct_low, pct_high, p_value, depends_on_habitat, dropped),
                  interactions_ab_500 %>% select(label, land_use, intensity, pct_diff, pct_low, pct_high, p_value))

build_lasso_sheet(wb, "R_LASSO_OPT", lasso_rch_optimal, optimal_lambda_bic_rch,
                  habitat_rch_opt   %>% select(term, land_use, intensity, pct_diff, pct_low, pct_high, p_value, n, dropped),
                  landscape_rch_opt %>% select(label, pct_diff, pct_low, pct_high, p_value, depends_on_habitat, dropped),
                  interactions_rch_opt %>% select(label, land_use, intensity, pct_diff, pct_low, pct_high, p_value))

sheet_index <- sheet_index %>% filter(!Sheet %in% c("AB_LASSO_OPT", "AB_LASSO_500", "R_LASSO_OPT"))
sheet_index <- bind_rows(sheet_index, data.frame(
  Sheet = c("AB_LASSO_OPT", "AB_LASSO_500", "R_LASSO_OPT"),
  Description = c("Abundance - glmmLasso model at BIC-optimal lambda",
                  "Abundance - glmmLasso model at lambda = 500 (fixed comparison)",
                  "Taxa richness - glmmLasso model (Poisson) at BIC-optimal lambda")))

writeData(wb, "Index", sheet_index, startRow = 3, startCol = 1)
saveWorkbook(wb, paste0(out_dir, "Model_Results.xlsx"), overwrite = TRUE)