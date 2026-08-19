#---------------------------------------------------------------------------------------------------------------------------------------#
# CODE DETAILS                                                                                                                          #
#                                                                                                                                       #
# Author: Catalina Vattuone  (cvattuonet@gmail.com)                                                                                     #        
# Date of latest update: 19-08-2026                                                                                                     #
# Main use: Study on "Assessing how agricultural practices and landscape composition shape differences in arthropod biodiversity among  #  
#           habitat types in Europe using the PREDICTS database". Internship at INRAE                                                   #
#                                                                                                                                       #
# Content: Biodiversity dataset creation and filtering for arthropods in Europe using the PREDICTS database                             #                           
#---------------------------------------------------------------------------------------------------------------------------------------#

#LIBRARIES  ------------------------------------------------------------
pacman::p_load(readr,  dplyr)

#LOAD  ------------------------------------------------------------
predicts_df  <- readRDS("Databases/PREDICTS/6fa1dedf-c546-41e0-a470-17c4863686b8.rds") #read the original rds file, downloaded from the PREDICTS website

#FILTERING AND ADDING RELEVANT INFORMATION TO THE OBSERVATION LEVEL DATAFRAME-----------------------------------------------------------
#filter by arthropoda phylum
observations_df <- predicts_df %>% filter(Phylum == "Arthropoda")

#filter by Insecta and Arachnida class 
observations_df <- observations_df %>% filter(Class %in% c("Insecta", "Arachnida"))

#filter by diversity metric type and unit
observations_df <- observations_df %>% filter(Diversity_metric_type %in% c("Abundance"))
observations_df <- observations_df %>% filter(Diversity_metric_unit %in% c("individuals"))

#filter by country (Europe and UK)
eu_uk_countries <- c("Austria", "Belgium", "Bulgaria", "Croatia", "Cyprus", "Czech Republic", "Denmark", "Estonia", "Finland", "France", 
                     "Germany", "Greece", "Hungary", "Ireland", "Italy", "Latvia", "Lithuania", "Luxembourg", "Malta", "Netherlands", 
                     "Poland", "Portugal", "Romania", "Slovakia", "Slovenia", "Spain", "Sweden", "United Kingdom")

observations_df <- observations_df %>% filter(Country %in% eu_uk_countries)

#filter by continental Europe
observations_df <- observations_df %>%
  filter(Longitude > -10 & Latitude > 35)

#add a column with the coordinates pair (latitude, longitude) for each site
observations_df <- observations_df %>%
  mutate(Lat_Long = paste(Latitude, Longitude, sep = ","))

#ad a start and finish year based on the Sample_start_earliest and Sample_end_latest columns
observations_df$Sample_start_year <- substr(observations_df$Sample_start_earliest, 1, 4)
observations_df$Sample_end_year <- substr(observations_df$Sample_end_latest, 1, 4)

#add the title of the source based on the Source_ID column
references_df<- read_delim("Databases/PREDICTS/resource.csv", delim = ",")

observations_df <- observations_df %>%
  left_join(references_df %>% select(Source_ID, Title) %>% 
      distinct(Source_ID, .keep_all = TRUE), by = "Source_ID")

#CREATING DATASET BY SITE WITH TOTAL ABUNDANCE-----------------------------------------------------------
sites_df <- observations_df %>%
  group_by(SSBS) %>%
  summarise(Source_ID = paste(unique(Source_ID), collapse = ", "),
            Reference = paste(unique(Reference), collapse = ", "),
            Title = paste(unique(Title), collapse = ", "),
            Study_number = paste(unique(Study_number), collapse = ", "),
            Study_name = paste(unique(Study_name), collapse = ", "),
            Diversity_metric = paste(unique(Diversity_metric), collapse = ", "),
            Diversity_metric_unit = paste(unique(Diversity_metric_unit), collapse = ", "),
            Sampling_method = paste(unique(Sampling_method), collapse = ", "),
            Site_number = paste(unique(Site_number), collapse = ", "),
            Site_name = paste(unique(Site_name), collapse = ", "),
            Block = paste(unique(Block), collapse = ", "),
            SS = paste(unique(SS), collapse = ", "),
            SSS = paste(unique(SSS), collapse = ", "),
            SSB = paste(unique(SSB), collapse = ", "),
            SSBS = paste(unique(SSBS), collapse = ", "),
            Sample_start_earliest = min(Sample_start_earliest, na.rm = TRUE),
            Sample_start_year = min(Sample_start_year, na.rm = TRUE),
            Sample_end_latest = max(Sample_end_latest, na.rm = TRUE),
            Sample_end_year = max(Sample_end_year, na.rm = TRUE),
            Sample_date_resolution = paste(unique(Sample_date_resolution), collapse = ", "),
            Max_linear_extent_metres = ifelse( all(is.na(Max_linear_extent_metres)), NA_real_,max(Max_linear_extent_metres, na.rm = TRUE) ),
            Sampling_effort = paste(unique(Sampling_effort), collapse = ", "),
            Rescaled_sampling_effort = paste(unique(Rescaled_sampling_effort), collapse = ", "),
            Habitat_as_described = paste(unique(Habitat_as_described), collapse = ", "),
            Predominant_land_use = paste(unique(Predominant_land_use), collapse = ", "),
            Source_for_predominant_land_use= paste(unique(Source_for_predominant_land_use), collapse =", "),
            Use_intensity= paste(unique(Use_intensity), collapse =", "),
            Coordinates_method= paste(unique(Coordinates_method), collapse =", "),
            Longitude= paste(unique(Longitude), collapse = ", "),
            Latitude= paste(unique(Latitude), collapse = ", "),
            Lat_Long = paste(unique(Lat_Long), collapse = ", "),
            Country = paste(unique(Country), collapse = ", "),
            Biome = paste(unique(Biome), collapse = ", "),
            Rank = paste(unique(Rank[Measurement > 0 & !is.na(Measurement) & !is.na(Rank) & Rank != ""]), collapse = ", "),
            N_class = n_distinct(Class[Measurement > 0 & !is.na(Measurement) & !is.na(Class) & Class != ""]),
            Class = paste(unique(Class[Measurement > 0 & !is.na(Measurement) & !is.na(Class) & Class != ""]), collapse = ", "),
            N_order = n_distinct(Order[Measurement > 0 & !is.na(Measurement) & !is.na(Order) & Order != ""]),
            Order = paste(unique(Order[Measurement > 0 & !is.na(Measurement) & !is.na(Order) & Order != ""]), collapse = ", "),
            N_family = n_distinct(Family[Measurement > 0 & !is.na(Measurement) & !is.na(Family) & Family != ""]),
            Family = paste(unique(Family[Measurement > 0 & !is.na(Measurement) & !is.na(Family) & Family != ""]), collapse = ", "),
            N_genus = n_distinct(Genus[Measurement > 0 & !is.na(Measurement) & !is.na(Genus) & Genus != ""]),
            Genus = paste(unique(Genus[Measurement > 0 & !is.na(Measurement) & !is.na(Genus) & Genus != ""]), collapse = ", "),
            N_species = n_distinct(Species[Measurement > 0 & !is.na(Measurement) & !is.na(Species) & Species != ""]),
            Species = paste(unique(Species[Measurement > 0 & !is.na(Measurement) & !is.na(Species) & Species != ""]), collapse = ", "),
            TA=sum(Measurement, na.rm=TRUE),
            TA_corrected=sum(Effort_corrected_measurement, na.rm=TRUE),
            .groups = "drop")


#the sites that do not have a SSB code we create one (same code for all sites of the SS)
sites_df <- sites_df %>% mutate(SSB = ifelse(is.na(Block), paste0(SS,"BLOCK1"), SSB))


#ESTIMATING TAXA AND GENUS RICHNESS BY SITE-----------------------------------------------------------
#estimate taxa richness by site
taxa_richness_by_site <- observations_df %>%
  filter(Measurement > 0) %>%
  mutate(Lowest_Taxon = case_when( #"Lowest_Taxon" column picks the most specific name available
    !is.na(Species) ~ Species,
    !is.na(Genus)   ~ Genus,
    !is.na(Family)  ~ Family,
    !is.na(Order)   ~ Order,
    TRUE            ~ Class)) %>%
  group_by(SSBS) %>%
  summarise(taxa_richness= n_distinct(Lowest_Taxon), #counts unique taxonomic entities
            .groups = "drop")

#merge taxa richness and total abundance with the sites dataframe
sites_df <- sites_df %>%
  left_join(taxa_richness_by_site %>% select(SSBS, taxa_richness), by = "SSBS")

#estimate genus richness by site
genus_richness <- observations_df %>%
  filter(Measurement > 0) %>%
  group_by(SSBS) %>%
  summarise(genus_richness= n_distinct(Genus), .groups = "drop")

sites_df <- sites_df %>%
  left_join(genus_richness, by = "SSBS")

#SAVE RESULTS  ------------------------------------------------------------
dir.create("Intermediate_dataset", showWarnings = FALSE)
write.csv(sites_df, "Intermediate_dataset/sites_df.csv",row.names = FALSE)
write.csv(observations_df, "Intermediate_dataset/observations_df.csv", row.names = FALSE)
