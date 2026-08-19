#---------------------------------------------------------------------------------------------------------------------------------------#
# CODE DETAILS                                                                                                                          #
#                                                                                                                                       #
# Author: Catalina Vattuone  (cvattuonet@gmail.com)                                                                                     #        
# Date of latest update: 19-08-2026                                                                                                     #
# Main use: Study on "Assessing how agricultural practices and landscape composition shape differences in arthropod biodiversity among  #  
#           habitat types in Europe using the PREDICTS database". Internship at INRAE                                                   #
#                                                                                                                                       #
# Content: Lasso models approach for abundance and richness models as complementary analysis                                            #                           
#---------------------------------------------------------------------------------------------------------------------------------------#

#LIBRARIES  ------------------------------------------------------------
pacman::p_load(glmmLasso,readr)

#LOAD DATABASES   ------------------------------------------------------------
sites_df <- read_delim("Intermediate_dataset/sites_for_abundance_models.csv") 
sites_df$Habitat_category <- relevel(factor(sites_df$Habitat_category), ref = "Natural vegetation")

sites_df_richness <- read_delim("Intermediate_dataset/sites_for_richness_models.csv") 
sites_df_richness$Habitat_category <- relevel(factor(sites_df_richness$Habitat_category), ref = "Natural vegetation")

#LASSO ABUNDANCE------------------------------------------------------------
#components
fixed_ab <- log_TA ~ Habitat_category + urban_pct_1000m_mod + water_pct_1000m_mod + 
  natural_pct_1000m_mod + pasture_pct_1000m_mod + crop_pct_1000m + 
  pesticide_1000m + fertilizer_1000m + 
  crop_diversity_shannon_1000m + crop_diversity_machefer_alpha_1000m + 
  pct_very_small_and_small_field_1000m + pct_very_large_and_large_field_1000m + 
  climate_mean_temp_c_1000m + climate_annual_precip_mm_1000m + 
  Habitat_category:pesticide_1000m +
  Habitat_category:fertilizer_1000m +
  Habitat_category:crop_diversity_shannon_1000m +
  Habitat_category:crop_diversity_machefer_alpha_1000m +
  Habitat_category:pct_very_small_and_small_field_1000m +
  Habitat_category:pct_very_large_and_large_field_1000m 

random_ab <- list(SS = ~1, SSB = ~1, Coordinate_ID = ~1, Sampling_method = ~1)

#optimal LAMBDA searching                     
lambda_grid_ab <- seq(from = 50, to = 2000, by = 50)
bic_values_ab <- rep(NA_real_, length(lambda_grid_ab))

for (i in 1:length(lambda_grid_ab)) {
  print(paste("Estimating for lambda", lambda_grid_ab[i]))
  bic_values_ab[i] <- tryCatch({
    invisible(capture.output(
      m_lasso <- glmmLasso(fix = fixed_ab, rnd = random_ab, data = as.data.frame(sites_df), family = gaussian(), lambda = lambda_grid_ab[i],control = glmmLassoControl())))
    m_lasso$bic
  }, error = function(e) {
    cat("  Failed at lambda =", lambda_grid_ab[i], ":", conditionMessage(e), "\n")
    NA_real_
  })
}

optimal_lambda_bic_ab <- lambda_grid_ab[which.min(bic_values_ab)]

#plot
out_dir <- "Figures/Lasso/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

png(paste0(out_dir, "lasso_bic_abundance.png"), width = 1600, height = 1000, res = 150)
par(cex.lab = 1.4, cex.axis = 1.2, cex.main = 1.4, cex = 1.2)
plot(lambda_grid_ab, bic_values_ab, type = "b", pch = 17, col = "darkblue",
     xlab = expression(lambda), ylab = "BIC",
     main = "LASSO Parameter Selection - BIC vs Lambda - Abundance")
abline(v = optimal_lambda_bic_ab, col = "red", lty = 2, lwd = 2)
dev.off()

#final models
lasso_ab_optimal <- glmmLasso(fix = fixed_ab, rnd = random_ab, data = as.data.frame(sites_df),  family = gaussian(), lambda = optimal_lambda_bic_ab,control = glmmLassoControl(bic = TRUE) )
lasso_ab_500 <- glmmLasso(fix = fixed_ab, rnd = random_ab, data = as.data.frame(sites_df),  family = gaussian(), lambda = 500,control = glmmLassoControl(bic = TRUE) )

#LASSO RICHNESS------------------------------------------------------------

#checking how far we are from poisson
mean(sites_df_richness$taxa_richness)
var(sites_df_richness$taxa_richness)

#data is definetely overdispersed so it needs and observation level random effect
sites_df_richness$obs_id <- factor(1:nrow(sites_df_richness))

#components
fixed_rch <- taxa_richness ~ Habitat_category + urban_pct_1000m_mod + water_pct_1000m_mod + 
  natural_pct_1000m_mod + pasture_pct_1000m_mod + crop_pct_1000m + 
  pesticide_1000m + fertilizer_1000m + 
  crop_diversity_shannon_1000m + crop_diversity_machefer_alpha_1000m + 
  pct_very_small_and_small_field_1000m + pct_very_large_and_large_field_1000m + 
  climate_mean_temp_c_1000m + climate_annual_precip_mm_1000m + 
  Habitat_category:pesticide_1000m +
  Habitat_category:fertilizer_1000m +
  Habitat_category:crop_diversity_shannon_1000m +
  Habitat_category:crop_diversity_machefer_alpha_1000m +
  Habitat_category:pct_very_small_and_small_field_1000m +
  Habitat_category:pct_very_large_and_large_field_1000m 

random_rch <- list(SS = ~1, SSB = ~1, Coordinate_ID = ~1, Sampling_method = ~1, obs_id = ~1)

#optimal LAMBDA searching                     
lambda_grid_rch <- seq(from = 0, to = 2000, by = 50) 
bic_values_rch <- numeric(length(lambda_grid_rch))

for (i in 1:length(lambda_grid_rch)) {
  print(paste("Estimating for lambda", lambda_grid_rch[i]))
  bic_values_rch[i] <- tryCatch({
    invisible(capture.output(
      m_lasso <- glmmLasso(fix = fixed_rch, rnd = random_rch, data = as.data.frame(sites_df_richness),family = poisson(), lambda = lambda_grid_rch[i],
                           control = glmmLassoControl())))
    m_lasso$bic
  }, error = function(e) {
    cat("  Failed at lambda =", lambda_grid_rch[i], ":", conditionMessage(e), "\n")
    NA_real_
  })
}

optimal_lambda_bic_rch <- lambda_grid_rch[which.min(bic_values_rch)]

#plot
png(paste0(out_dir, "lasso_bic_richness.png"), width = 1600, height = 1000, res = 150)
par(cex.lab = 1.4, cex.axis = 1.2, cex.main = 1.4, cex = 1.2)
plot(lambda_grid_rch, bic_values_rch, type = "b", pch = 17, col = "darkblue",
     xlab = expression(lambda), ylab = "BIC",
     main = "LASSO Parameter Selection - BIC vs Lambda - Taxa Richness")
abline(v = optimal_lambda_bic_rch, col = "red", lty = 2, lwd = 2)
dev.off()


#second search optimal LAMBDA searching                     
lambda_grid_rch <- seq(from = 300, to = 500, by = 10) 
bic_values_rch <- numeric(length(lambda_grid_rch))

for (i in 1:length(lambda_grid_rch)) {
  print(paste("Estimating for lambda", lambda_grid_rch[i]))
  bic_values_rch[i] <- tryCatch({
    invisible(capture.output(
      m_lasso <- glmmLasso(fix = fixed_rch, rnd = random_rch, data = as.data.frame(sites_df_richness),
                           family = poisson(), lambda = lambda_grid_rch[i],
                           control = glmmLassoControl())))
    m_lasso$bic
  }, error = function(e) {
    cat("  Failed at lambda =", lambda_grid_rch[i], ":", conditionMessage(e), "\n")
    NA_real_
  })
}

optimal_lambda_bic_rch <- lambda_grid_rch[which.min(bic_values_rch)]

#plot
png(paste0(out_dir, "lasso_bic_richness_zoomed.png"), width = 1600, height = 1000, res = 150)
par(cex.lab = 1.4, cex.axis = 1.2, cex.main = 1.4, cex = 1.2)
plot(lambda_grid_rch, bic_values_rch, type = "b", pch = 17, col = "darkblue",
     xlab = expression(lambda), ylab = "BIC",
     main = "LASSO Parameter Selection - BIC vs Lambda - Taxa Richness Zoomed")
abline(v = optimal_lambda_bic_rch, col = "red", lty = 2, lwd = 2)
dev.off()

#final model
lasso_rch_optimal <- glmmLasso(fix = fixed_rch, rnd = random_rch, data = as.data.frame(sites_df_richness),  family = poisson(), lambda = optimal_lambda_bic_rch,control = glmmLassoControl(bic = TRUE) )

#Saving results------------------------------------------------------------
dir.create("Models_results", showWarnings = FALSE)
saveRDS(lasso_ab_optimal,  "Models_results/lasso_ab_optimal.rds")
saveRDS(lasso_ab_500,      "Models_results/lasso_ab_500.rds")
saveRDS(lasso_rch_optimal, "Models_results/lasso_rch_optimal.rds")