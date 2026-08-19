#---------------------------------------------------------------------------------------------------------------------------------------#
# CODE DETAILS                                                                                                                          #
#                                                                                                                                       #
# Author: Catalina Vattuone  (cvattuonet@gmail.com)                                                                                     #        
# Date of latest update: 19-08-2026                                                                                                     #
# Main use: Study on "Assessing how agricultural practices and landscape composition shape differences in arthropod biodiversity among  #  
#           habitat types in Europe using the PREDICTS database". Internship at INRAE                                                   #
#                                                                                                                                       #
# Content: Running Null Models for deciding on random effects.                                                                          #                           
#---------------------------------------------------------------------------------------------------------------------------------------#


#LIBRARIES  ------------------------------------------------------------
pacman::p_load(readr, lme4, glmmTMB)

#LOAD   ------------------------------------------------------------
sites_df <- read_delim("Intermediate_dataset/sites_for_abundance_models.csv") 
sites_df_richness <- read_delim("Intermediate_dataset/sites_for_richness_models.csv") 

#ABUNDANCE NULL MODELS ------------------------------------------------------------
#gaussian
null_1 <- lmer(log_TA ~ 1 + (1 | SS) + (1 | SSB) , data = sites_df, REML = TRUE)
null_2 <- lmer(log_TA ~ 1 + (1 | SS) + (1 | SSB) + (1 | Coordinate_ID), data = sites_df, REML = TRUE)
null_3 <- lmer(log_TA ~ 1 + (1 | SS) + (1 | SSB) + (1 | Coordinate_ID) + (1 | Sample_start_year), data = sites_df, REML = TRUE)
null_4 <- lmer(log_TA ~ 1 + (1 | SS) + (1 | SSB) + (1 | Coordinate_ID) + (1 | Sampling_method),  data = sites_df, REML = TRUE)
null_5 <- lmer(log_TA ~ 1 + (1 | SS) + (1 | SSB) + (1 | Coordinate_ID) + (1 | Sample_start_year) + (1 | Sampling_method),  data = sites_df, REML = TRUE)

AIC(null_1, null_2, null_3, null_4, null_5)

#RICHNESS NULL MODELS ------------------------------------------------------------
#zero inflated negative binomial
rich_null_1 <- glmmTMB(taxa_richness ~ 1 + (1 | SS) + (1 | SSB) , data = sites_df_richness,family = nbinom2, ziformula = ~ 1)
rich_null_2 <- glmmTMB(taxa_richness ~ 1 + (1 | SS) + (1 | SSB) + (1 | Coordinate_ID) ,data = sites_df_richness,family = nbinom2, ziformula = ~ 1)
rich_null_3 <- glmmTMB(taxa_richness ~ 1 + (1 | SS) + (1 | SSB) + (1 | Coordinate_ID) + (1 | Sample_start_year),data = sites_df_richness,family = nbinom2, ziformula = ~ 1)
rich_null_4 <- glmmTMB(taxa_richness ~ 1 + (1 | SS) + (1 | SSB) + (1 | Coordinate_ID) + (1 | Sampling_method),data = sites_df_richness,family = nbinom2, ziformula = ~ 1)
rich_null_5 <- glmmTMB(taxa_richness ~ 1 + (1 | SS) + (1 | SSB) + (1 | Coordinate_ID) +  (1 | Sample_start_year) + (1 | Sampling_method),data = sites_df_richness,family = nbinom2, ziformula = ~ 1)

AIC(rich_null_1, rich_null_2, rich_null_3, rich_null_4, rich_null_5)