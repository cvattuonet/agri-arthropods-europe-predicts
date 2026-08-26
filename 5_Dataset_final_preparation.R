#---------------------------------------------------------------------------------------------------------------------------------------#
# CODE DETAILS                                                                                                                          #
#                                                                                                                                       #
# Author: Catalina Vattuone  (cvattuonet@gmail.com)                                                                                     #        
# Date of latest update: 19-08-2026                                                                                                     #
# Main use: Study on "Assessing how agricultural practices and landscape composition shape differences in arthropod biodiversity among  #  
#           habitat types in Europe using the PREDICTS database". Internship at INRAE                                                   #
#                                                                                                                                       #
# Content: Final preparation of the biodiversity dataset for modelling                                                                  #                           
#---------------------------------------------------------------------------------------------------------------------------------------#


#LIBRARIES  ------------------------------------------------------------
pacman::p_load(readr, dplyr, tidyr)

#LOAD-----------------------------------------------------------
sites_df<- read_delim("Intermediate_dataset/sites_df_predictors.csv") 
observations_df <- read_delim("Intermediate_dataset/observations_mod.csv")

#HABITAT CATEGORY CREATION-----------------------------------------------------------
#creating the category
sites_df <- sites_df %>% mutate( Habitat_category = paste(Predominant_land_use_mod,Use_intensity,sep = "_"))

#adding it to observations
observations_df <-  observations_df %>%
  left_join(sites_df %>% select(SSBS, Habitat_category), by = "SSBS")

#redefine Habitat "Natural vegetation_Light use" and "Natural vegetation_Minimal use"  to "Natural vegetation"
sites_df <- sites_df %>%
  mutate(Habitat_category = ifelse(Habitat_category %in% c("Natural vegetation_Light use", "Natural vegetation_Minimal use"), "Natural vegetation", Habitat_category))

#drop cannot decide sites
sites_df <- sites_df %>% filter(Use_intensity != "Cannot decide")

#CONVERTING TO FACTOR VALUES ------------------------------------------------------------
sites_df <- sites_df %>%
  mutate(Coordinate_ID = as.factor(Lat_Long_mod),
         SS = as.factor(SS),
         SSB = as.factor(SSB),
         SSBS = as.factor(SSBS),
         Sample_start_year = as.factor(Sample_start_year),
         Sampling_method = as.factor(Sampling_method),
         Use_intensity = droplevels(factor(Use_intensity)),
         Predominant_land_use_mod = droplevels(factor(Predominant_land_use_mod)))

#RESPONSE VARIABLES ------------------------------------------------------------
#abundance
#creating the response variable (log (ab+1))
sites_df <- sites_df %>% mutate(log_TA = log(TA_corrected_mod + 1))

#richness
#all the NA values in genus_richness and taxa_richness need to be transform to 0 (they are NA because there were only observations with abundance 0)
sites_df <- sites_df %>%
  mutate(genus_richness = ifelse(is.na(genus_richness), 0, genus_richness),
         taxa_richness = ifelse(is.na(taxa_richness), 0, taxa_richness))

#we need to take out the references that are only for abundance, given what was agreed on the review of primary literature.
references_only_abundance = c("Darvill et al. 2004", "Nielsen et al. 2011")  #Nielsen et al. 2011 is already out (due to the recategorization of land uses)
sites_df_for_richness <- sites_df  %>% #we take out 2 references for the richness assessment given what was agreed on the review of primary literature.
  filter(!Reference %in% references_only_abundance)
  
#CREATING THE DATABASE FOR MODELLING------------------------------------------------------------
categorical_predictors= c("Habitat_category","Predominant_land_use_mod","Use_intensity") 

continuous_predictors<- c("urban_pct_1000m_mod", "water_pct_1000m_mod", "natural_pct_1000m_mod", 
                                "pasture_pct_1000m_mod", "crop_pct_1000m", "pesticide_1000m", "fertilizer_1000m", 
                                "crop_diversity_shannon_1000m", "crop_diversity_machefer_alpha_1000m", 
                                "pct_very_small_and_small_field_1000m", "pct_very_large_and_large_field_1000m", 
                                "climate_mean_temp_c_1000m", "climate_annual_precip_mm_1000m")

random_effect=c("SS", "SSB", "Coordinate_ID","Sampling_method", "SSBS")
response_abundance= c("log_TA", "TA_corrected_mod","TA")
response_richnes= c("taxa_richness", "genus_richness")
others = c("Rescaled_sampling_effort_mod", "Country", "Sample_start_year", "Diversity_metric_unit")

df_1000 <- sites_df %>%select(all_of(c(continuous_predictors, categorical_predictors, random_effect, response_abundance,others)))
df_rich_1000 <- sites_df_for_richness %>% select(all_of(c(continuous_predictors, categorical_predictors, random_effect, response_richnes, others)))

#REESCALING FROM 0 TO 1------------------------------------------------------------
df_1000_scaled <- df_1000
df_1000_scaled[continuous_predictors] <- lapply(df_1000_scaled[continuous_predictors], scale)

df_rich_1000_scaled <- df_rich_1000
df_rich_1000_scaled[continuous_predictors] <- lapply(df_rich_1000_scaled[continuous_predictors], scale)


#CLEANING NA VALUES------------------------------------------------------------
#inspection of NA values
na_df_1000_scaled  <- tibble(column = names(df_1000_scaled), na_count = colSums(is.na(df_1000_scaled))) 
na_df_rich_1000_scaled  <- tibble(column = names(df_rich_1000_scaled), na_count = colSums(is.na(df_rich_1000_scaled))) 
print(na_df_1000_scaled)
print(na_df_rich_1000_scaled)

#clean NA values given all come from crop_diversity_shannon_1000m
df_1000_scaled_clean <- df_1000_scaled %>%  drop_na(crop_diversity_shannon_1000m)
df_rich_1000_scaled_clean <- df_rich_1000_scaled %>%  drop_na(crop_diversity_shannon_1000m)

#DATAFRAME FOR VIZUALIZATION - NON SCALE-----------------------------------------------------------
df_1000_viz <- df_1000 %>% drop_na(crop_diversity_shannon_1000m)

#DEFINE REFERENCE VARIABLE-----------------------------------------------------------
df_1000_scaled_clean$Predominant_land_use_mod <- relevel(factor(df_1000_scaled_clean$Predominant_land_use_mod), ref = "Natural vegetation")
df_1000_scaled_clean$Habitat_category <- relevel(factor(df_1000_scaled_clean$Habitat_category), ref = "Natural vegetation")
df_1000_scaled_clean$Use_intensity <- relevel(factor(df_1000_scaled_clean$Use_intensity), ref = "Minimal use")

df_rich_1000_scaled_clean$Predominant_land_use_mod <- relevel(factor(df_rich_1000_scaled_clean$Predominant_land_use_mod), ref = "Natural vegetation")
df_rich_1000_scaled_clean$Habitat_category <- relevel(factor(df_rich_1000_scaled_clean$Habitat_category), ref = "Natural vegetation")
df_rich_1000_scaled_clean$Use_intensity <- relevel(factor(df_rich_1000_scaled_clean$Use_intensity), ref = "Minimal use")

#TAKE OUT NON INDIVIDUALS VALUES-----------------------------------------------------------
df_1000_scaled_clean <- df_1000_scaled_clean %>% filter(Diversity_metric_unit == "individuals")
df_rich_1000_scaled_clean <- df_rich_1000_scaled_clean %>% filter(Diversity_metric_unit == "individuals")

#SAVE-----------------------------------------------------------
write.csv(df_1000_scaled_clean, "Intermediate_dataset/sites_for_abundance_models.csv", row.names = FALSE)
write.csv(df_rich_1000_scaled_clean, "Intermediate_dataset/sites_for_richness_models.csv", row.names = FALSE)
write.csv(df_1000_viz, "Intermediate_dataset/abundance_sites_for_vizualization_nonscale.csv", row.names = FALSE)