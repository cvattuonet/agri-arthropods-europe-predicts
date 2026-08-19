#---------------------------------------------------------------------------------------------------------------------------------------#
# CODE DETAILS                                                                                                                          #
#                                                                                                                                       #
# Author: Catalina Vattuone  (cvattuonet@gmail.com)                                                                                     #        
# Date of latest update: 19-08-2026                                                                                                     #
# Main use: Study on "Assessing how agricultural practices and landscape composition shape differences in arthropod biodiversity among  #  
#           habitat types in Europe using the PREDICTS database". Internship at INRAE                                                   #
#                                                                                                                                       #
# Content: Mixed effects models running: Abundance and - Richness Models with dredge function                                           #                           
#---------------------------------------------------------------------------------------------------------------------------------------#


#LIBRARIES  ------------------------------------------------------------
pacman::p_load(readr,lme4,glmmTMB,tweedie,parallel,MuMIn)

#LOAD   ------------------------------------------------------------
sites_df <- read_delim("Intermediate_dataset/sites_for_abundance_models.csv") 
sites_df$Habitat_category <- relevel(factor(sites_df$Habitat_category), ref = "Natural vegetation")

sites_df_richness <- read_delim("Intermediate_dataset/sites_for_richness_models.csv") 
sites_df_richness$Habitat_category <- relevel(factor(sites_df_richness$Habitat_category), ref = "Natural vegetation")

#COMMENTS -----------------------------------------------------------
#For each case (both abundance and richness), I run first the global model andthen the dredge ranking by AICc. 
#The extraction of the results is done in a separate script


#ABUNDANCE MODELS-----------------------------------------------------------
##GAUSSIAN-----------------------------------------------------------
#given I am running dredge() models and it is swll, I build this to make sure I am using all the cores available in my computer.
detectCores()
my_cluster <- makeCluster(detectCores() - 1)
clusterEvalQ(my_cluster, library(lme4))   #for lme4
clusterExport(my_cluster, varlist = "sites_df")

global_gaussian <- lmer(log_TA ~ Habitat_category + urban_pct_1000m_mod + water_pct_1000m_mod + natural_pct_1000m_mod + pasture_pct_1000m_mod +
                         pesticide_1000m +crop_diversity_machefer_alpha_1000m + pct_very_small_and_small_field_1000m +pct_very_large_and_large_field_1000m +
                         Habitat_category:natural_pct_1000m_mod+Habitat_category:water_pct_1000m_mod +Habitat_category:pesticide_1000m +
                         Habitat_category:crop_diversity_machefer_alpha_1000m + Habitat_category:pct_very_small_and_small_field_1000m +climate_mean_temp_c_1000m +
                         (1 | SS) + (1 | SSB) + (1 | Coordinate_ID) + (1 | Sampling_method), data = sites_df,REML = FALSE,na.action = "na.fail")

dredge_gaussian <- dredge(global_gaussian, cluster = my_cluster, rank = "AICc",
                         subset = dc(pesticide_1000m, `Habitat_category:pesticide_1000m`) &&
                           dc(natural_pct_1000m_mod, `Habitat_category:natural_pct_1000m_mod`) &&
                           dc(water_pct_1000m_mod, `Habitat_category:water_pct_1000m_mod`) &&
                           dc(crop_diversity_machefer_alpha_1000m, `Habitat_category:crop_diversity_machefer_alpha_1000m`) &&
                           dc(pct_very_small_and_small_field_1000m, `Habitat_category:pct_very_small_and_small_field_1000m`), trace = 2)

dir.create("Models_results", showWarnings = FALSE)
saveRDS(global_gaussian, "Models_results/global_gaussian.rds")
saveRDS(dredge_gaussian, "Models_results/dredge_gaussian.rds")

stopCluster(my_cluster)

##POISSON-----------------------------------------------------------
detectCores()
my_cluster_2 <- makeCluster(detectCores() - 1)
clusterEvalQ(my_cluster_2, library(glmmTMB))
clusterExport(my_cluster_2, varlist = "sites_df")

global_poisson <- glmmTMB(TA ~ Habitat_category + urban_pct_1000m_mod + water_pct_1000m_mod + natural_pct_1000m_mod + 
                                pasture_pct_1000m_mod + pesticide_1000m + crop_diversity_machefer_alpha_1000m + climate_mean_temp_c_1000m +
                                pct_very_small_and_small_field_1000m + pct_very_large_and_large_field_1000m +
                                Habitat_category:natural_pct_1000m_mod+Habitat_category:water_pct_1000m_mod +Habitat_category:pesticide_1000m + 
                                Habitat_category:crop_diversity_machefer_alpha_1000m +Habitat_category:pct_very_small_and_small_field_1000m + 
                                offset(log(Rescaled_sampling_effort_mod)) +(1 | SS) + (1 | SSB) + (1 | Coordinate_ID) + (1 | Sampling_method),
                              family = poisson(), data = sites_df, na.action = "na.fail")

dredge_poisson <- dredge(global_poisson , cluster = my_cluster_2, rank = "AICc",
                             subset =  dc(pesticide_1000m, `Habitat_category:pesticide_1000m`) &&
                               dc(natural_pct_1000m_mod, `Habitat_category:natural_pct_1000m_mod`) &&
                               dc(water_pct_1000m_mod, `Habitat_category:water_pct_1000m_mod`) &&
                               dc(crop_diversity_machefer_alpha_1000m, `Habitat_category:crop_diversity_machefer_alpha_1000m`) &&
                               dc(pct_very_small_and_small_field_1000m, `Habitat_category:pct_very_small_and_small_field_1000m`),trace = 2)

saveRDS(global_poisson, "Models_results/global_poisson.rds")
saveRDS(dredge_poisson, "Models_results/dredge_poisson.rds")

##negative binomial-----------------------------------------------------------
global_nb <- glmmTMB(TA ~ Habitat_category + urban_pct_1000m_mod + water_pct_1000m_mod + natural_pct_1000m_mod + 
                               pasture_pct_1000m_mod + pesticide_1000m + crop_diversity_machefer_alpha_1000m + 
                               pct_very_small_and_small_field_1000m + pct_very_large_and_large_field_1000m +
                               Habitat_category:natural_pct_1000m_mod+
                               Habitat_category:water_pct_1000m_mod +
                               Habitat_category:pesticide_1000m + 
                               Habitat_category:crop_diversity_machefer_alpha_1000m +
                               Habitat_category:pct_very_small_and_small_field_1000m + 
                               climate_mean_temp_c_1000m +
                               offset(log(Rescaled_sampling_effort_mod)) +
                               (1 | SS) + (1 | SSB) + (1 | Coordinate_ID) + (1 | Sampling_method),
                             family = nbinom2(), data = sites_df, na.action = "na.fail")

dredge_nb <- dredge(global_nb, cluster = my_cluster_2, rank = "AICc",
                            subset = dc(pesticide_1000m, `Habitat_category:pesticide_1000m`) &&
                              dc(natural_pct_1000m_mod, `Habitat_category:natural_pct_1000m_mod`) &&
                              dc(water_pct_1000m_mod, `Habitat_category:water_pct_1000m_mod`) &&
                              dc(crop_diversity_machefer_alpha_1000m, `Habitat_category:crop_diversity_machefer_alpha_1000m`) &&
                              dc(pct_very_small_and_small_field_1000m, `Habitat_category:pct_very_small_and_small_field_1000m`), trace = 2)

saveRDS(global_nb, "Models_results/global_nb.rds")
saveRDS(dredge_nb, "Models_results/dredge_nb.rds")

##tweedie-----------------------------------------------------------
global_tweedie <- glmmTMB(TA ~ Habitat_category + urban_pct_1000m_mod + water_pct_1000m_mod + natural_pct_1000m_mod + 
                               pasture_pct_1000m_mod + pesticide_1000m + crop_diversity_machefer_alpha_1000m + 
                               pct_very_small_and_small_field_1000m + pct_very_large_and_large_field_1000m +
                               Habitat_category:natural_pct_1000m_mod+
                               Habitat_category:water_pct_1000m_mod +
                               Habitat_category:pesticide_1000m + 
                               Habitat_category:crop_diversity_machefer_alpha_1000m +
                               Habitat_category:pct_very_small_and_small_field_1000m + 
                               climate_mean_temp_c_1000m +
                               offset(log(Rescaled_sampling_effort_mod)) +
                               (1 | SS) + (1 | SSB) + (1 | Coordinate_ID) + (1 | Sampling_method),
                             family = tweedie(link = "log"), data = sites_df, na.action = "na.fail")

dredge_tweedie <- dredge(global_tweedie, cluster = my_cluster_2, rank = "AICc",
                         subset = dc(pesticide_1000m, `Habitat_category:pesticide_1000m`) &&
                           dc(natural_pct_1000m_mod, `Habitat_category:natural_pct_1000m_mod`) &&
                           dc(water_pct_1000m_mod, `Habitat_category:water_pct_1000m_mod`) &&
                           dc(crop_diversity_machefer_alpha_1000m, `Habitat_category:crop_diversity_machefer_alpha_1000m`) &&
                           dc(pct_very_small_and_small_field_1000m, `Habitat_category:pct_very_small_and_small_field_1000m`), trace = 2)

saveRDS(global_tweedie, "Models_results/global_tweedie.rds")
saveRDS(dredge_tweedie, "Models_results/dredge_tweedie.rds")


##Zero inflated negative binomial-----------------------------------------------------------
global_zinb <- glmmTMB(TA ~ Habitat_category + urban_pct_1000m_mod + water_pct_1000m_mod + natural_pct_1000m_mod + 
                               pasture_pct_1000m_mod + pesticide_1000m + crop_diversity_machefer_alpha_1000m + 
                               pct_very_small_and_small_field_1000m + pct_very_large_and_large_field_1000m +
                               Habitat_category:natural_pct_1000m_mod+
                               Habitat_category:water_pct_1000m_mod +
                               Habitat_category:pesticide_1000m + 
                               Habitat_category:crop_diversity_machefer_alpha_1000m +
                               Habitat_category:pct_very_small_and_small_field_1000m + 
                               climate_mean_temp_c_1000m +
                               offset(log(Rescaled_sampling_effort_mod)) +
                               (1 | SS) + (1 | SSB) + (1 | Coordinate_ID) + (1 | Sampling_method),
                             ziformula = ~1, family = nbinom2(), data = sites_df, na.action = "na.fail")

dredge_zinb <- dredge(global_zinb, cluster = my_cluster_2, rank = "AICc",
                            subset = dc(pesticide_1000m, `Habitat_category:pesticide_1000m`) &&
                              dc(natural_pct_1000m_mod, `Habitat_category:natural_pct_1000m_mod`) &&
                              dc(water_pct_1000m_mod, `Habitat_category:water_pct_1000m_mod`) &&
                              dc(crop_diversity_machefer_alpha_1000m, `Habitat_category:crop_diversity_machefer_alpha_1000m`) &&
                              dc(pct_very_small_and_small_field_1000m, `Habitat_category:pct_very_small_and_small_field_1000m`), trace = 2)

saveRDS(global_zinb, "Models_results/global_zinb.rds")
saveRDS(dredge_zinb, "Models_results/dredge_zinb.rds")

stopCluster(my_cluster_2)

#RICHNESS MODELS-----------------------------------------------------------
detectCores()
my_cluster_3 <- makeCluster(detectCores() - 1)
clusterEvalQ(my_cluster_3, library(glmmTMB))
clusterExport(my_cluster_3, varlist = "sites_df_richness")

##Negative binomial richness-----------------------------------------------------------
global_nb_richness <- glmmTMB(taxa_richness ~ Habitat_category + urban_pct_1000m_mod + water_pct_1000m_mod + natural_pct_1000m_mod +
                              pasture_pct_1000m_mod + pesticide_1000m + crop_diversity_machefer_alpha_1000m +
                              pct_very_small_and_small_field_1000m + pct_very_large_and_large_field_1000m +
                              climate_mean_temp_c_1000m +
                              Habitat_category:natural_pct_1000m_mod+
                              Habitat_category:water_pct_1000m_mod +
                              Habitat_category:pesticide_1000m +
                              Habitat_category:crop_diversity_machefer_alpha_1000m +
                              Habitat_category:pct_very_small_and_small_field_1000m  +
                              (1 | SS) + (1 | SSB) + (1 | Coordinate_ID) + (1 | Sampling_method),
                            family = nbinom2(), data = sites_df_richness, na.action = "na.fail")


dredge_nb_richness <- dredge(global_nb_richness , cluster = my_cluster_3, rank = "AICc",
                           subset =
                             dc(pesticide_1000m, `Habitat_category:pesticide_1000m`) &&
                             dc(natural_pct_1000m_mod, `Habitat_category:natural_pct_1000m_mod`) &&
                             dc(water_pct_1000m_mod, `Habitat_category:water_pct_1000m_mod`) &&
                             dc(crop_diversity_machefer_alpha_1000m, `Habitat_category:crop_diversity_machefer_alpha_1000m`) &&
                             dc(pct_very_small_and_small_field_1000m, `Habitat_category:pct_very_small_and_small_field_1000m`), trace = 2)

saveRDS(global_nb_richness, "Models_results/global_nb_richness.rds")
saveRDS(dredge_nb_richness, "Models_results/dredge_nb_richness.rds")

##Zero inflated negative binomial richness-----------------------------------------------------------
global_zinb_richness <- glmmTMB(taxa_richness ~ Habitat_category + urban_pct_1000m_mod + water_pct_1000m_mod + natural_pct_1000m_mod +
                               pasture_pct_1000m_mod + pesticide_1000m + crop_diversity_machefer_alpha_1000m +
                               pct_very_small_and_small_field_1000m + pct_very_large_and_large_field_1000m +
                               climate_mean_temp_c_1000m +
                               Habitat_category:natural_pct_1000m_mod+
                               Habitat_category:water_pct_1000m_mod +
                               Habitat_category:pesticide_1000m +
                               Habitat_category:crop_diversity_machefer_alpha_1000m +
                               Habitat_category:pct_very_small_and_small_field_1000m +
                               (1 | SS) + (1 | SSB) + (1 | Coordinate_ID) + (1 | Sampling_method),
                             ziformula = ~1, family = nbinom2(), data = sites_df_richness, na.action = "na.fail")

dredge_zinb_richness <- dredge(global_zinb_richness , cluster = my_cluster_3, rank = "AICc",
                            subset =
                              dc(pesticide_1000m, `Habitat_category:pesticide_1000m`) &&
                              dc(natural_pct_1000m_mod, `Habitat_category:natural_pct_1000m_mod`) &&
                              dc(water_pct_1000m_mod, `Habitat_category:water_pct_1000m_mod`) &&
                              dc(crop_diversity_machefer_alpha_1000m, `Habitat_category:crop_diversity_machefer_alpha_1000m`) &&
                              dc(pct_very_small_and_small_field_1000m, `Habitat_category:pct_very_small_and_small_field_1000m`), trace = 2)


saveRDS(global_zinb_richness, "Models_results/global_zinb_richness.rds")
saveRDS(dredge_zinb_richness, "Models_results/dredge_zinb_richness.rds")

stopCluster(my_cluster_3)
