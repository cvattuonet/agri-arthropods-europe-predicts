#---------------------------------------------------------------------------------------------------------------------------------------#
# CODE DETAILS                                                                                                                          #
#                                                                                                                                       #
# Author: Catalina Vattuone  (cvattuonet@gmail.com)                                                                                     #        
# Date of latest update: 19-08-2026                                                                                                     #
# Main use: Study on "Assessing how agricultural practices and landscape composition shape differences in arthropod biodiversity among  #  
#           habitat types in Europe using the PREDICTS database". Internship at INRAE                                                   #
#                                                                                                                                       #
# Content: Biodiversity dataset curation and correction based on primary literature review                                              #                           
#---------------------------------------------------------------------------------------------------------------------------------------#

#LIBRARIES  ------------------------------------------------------------
pacman::p_load(readr,dplyr,stringr)

#LOAD  ------------------------------------------------------------
sites_df <- read_delim("Intermediate_dataset/sites_df.csv")
observations_df <- read_delim("Intermediate_dataset/observations_df.csv")

#PREPARATION ------------------------------------------------------------
sites_mod <- sites_df

sites_mod  <- sites_mod  %>%
  mutate(Study_number_mod = Study_number,
         Site_number_mod = Site_number,
         Site_name_mod = Site_name,
         Block_mod = Block,
         Rescaled_sampling_effort_mod = Rescaled_sampling_effort,
         TA_corrected_mod = TA_corrected,
         Sampling_method_mod = Sampling_method,
         Longitude_mod = Longitude, 
         Latitude_mod = Latitude)

sites_mod <- sites_mod %>%
  select(Source_ID, Reference, Title, Study_number, Study_number_mod, Study_name,
         Site_number, Site_number_mod, Site_name, Site_name_mod, Block, Block_mod, Sampling_method, Sampling_method_mod, TA, TA_corrected, TA_corrected_mod, everything())

#CURATION BY REFERENCE ------------------------------------------------------------
## BERG ET AL 2011  ------------------------------------------------------------
target_row <- which(sites_mod$Reference == "Berg et al. 2011" &  sites_mod$Site_name == "HW_Berg_2011_R263_ABO_E")
sites_mod[target_row, "Predominant_land_use"] <- "Pasture"

## BILLETER ET AL 2008  ------------------------------------------------------------
#Fixing the coordinates of the pitfall traps sites. 
#New coordinates in the Longitud_mod and Latitude_mod based on the site number and country, by looking at the flight trap values.  
flight_trap_coords <- sites_mod %>%
  filter(Reference == "Billeter et al. 2008", Sampling_method == "flight trap") %>%
  distinct(Country, Site_number, .keep_all = TRUE) %>%
  select(Country, Site_number, True_Longitude = Longitude_mod, True_Latitude = Latitude_mod)

sites_mod <- sites_mod %>%
  left_join(flight_trap_coords, by = c("Country", "Site_number")) %>%
  mutate(Longitude_mod = if_else( Reference == "Billeter et al. 2008" & Sampling_method == "pit-fall traps" & !is.na(True_Longitude),True_Longitude,Longitude_mod),
    Latitude_mod = if_else(Reference == "Billeter et al. 2008" & Sampling_method == "pit-fall traps" & !is.na(True_Latitude),True_Latitude,Latitude_mod) ) %>%
  select(-True_Longitude, -True_Latitude)

## BLAKE ET AL 2011  ------------------------------------------------------------
#Take out the reference
sites_mod <- sites_mod %>%filter(!(Reference == "Blake et al. 2011" ))

## CARPENTER ET AL. 2012  ------------------------------------------------------------
#Take out the reference
sites_mod <- sites_mod %>%filter(!(Reference == "Carpenter et al. 2012" ))

## DAVIS ET AL 2010  ------------------------------------------------------------
#Take out the reference
sites_mod <- sites_mod %>% filter(Reference != "Davis et al. 2010")

## DIEKOTTER ET AL 2010  -----------------------------------------------------------
#Take out the reference
sites_mod <- sites_mod %>% filter(Reference != "Diekötter et al. 2006")

## GAUBLOME ET AL 2008  ------------------------------------------------------------
#Take out the  sites with inconsistencies
sites_mod <- sites_mod %>%
  filter(!(Reference == "Gaublomme et al. 2008" & Site_name %in% c("Soignes urban", "Soignes suburban", "Soignes rural", "Soignes extra", "Verrewinkel", "Kleet")))


## HANLEY ET AL 2011  ------------------------------------------------------------
#Giving separate study numbers to each year
sites_mod <- sites_mod %>%
  mutate(Study_number_mod = ifelse(Reference == "Hanley et al. 2011" & Sample_start_year == 2007, 1, Study_number_mod))
sites_mod <- sites_mod %>%
  mutate(Study_number_mod = ifelse(Reference == "Hanley et al. 2011" & Sample_start_year == 2008, 2, Study_number_mod))
sites_mod <- sites_mod %>%
  mutate(Study_number_mod = ifelse(Reference == "Hanley et al. 2011" & Sample_start_year == 2009, 3, Study_number_mod))
sites_mod <- sites_mod %>%
  mutate(Study_number_mod = ifelse(Reference == "Hanley et al. 2011" & Sample_start_year == 2010, 4, Study_number_mod))
sites_mod <- sites_mod %>%
  mutate(Study_number_mod = ifelse(Reference == "Hanley et al. 2011" & Site_name == "SS 514 144_08", 1, Study_number_mod))

## HERRMANN ET AL 2007  ------------------------------------------------------------
#Take out the reference
sites_mod <- sites_mod %>% filter(!(Reference == "Herrmann et al. 2007" ))

## KOHLER ET AL 2008  ------------------------------------------------------------
#Take out the reference
sites_mod <- sites_mod %>%filter(!(Reference == "Kohler et al. 2008" ))

## JONSELL 2012  ------------------------------------------------------------
#Giving separate study numbers to each year
sites_mod <- sites_mod %>%
  mutate(Study_number_mod = ifelse(Reference == "Jonsell 2012" & Sample_start_year == 2001, 1, Study_number_mod))
sites_mod <- sites_mod %>%
  mutate(Study_number_mod = ifelse(Reference == "Jonsell 2012" & Sample_start_year == 2004, 2, Study_number_mod))
sites_mod <- sites_mod %>%
  mutate(Study_number_mod = ifelse(Reference == "Jonsell 2012" & Sample_start_year == 2006, 3, Study_number_mod))
sites_mod <- sites_mod %>%
  mutate(Study_number_mod = ifelse(Reference == "Jonsell 2012" & Sample_start_year == 2007, 4, Study_number_mod))
sites_mod <- sites_mod %>%
  mutate(Study_number_mod = ifelse(Reference == "Jonsell 2012" & Sample_start_year == 2008, 5, Study_number_mod))

## LEIGHTON-GOODALL ET AL 2012  ------------------------------------------------------------
#Take out the sites from 2003
sites_mod <- sites_mod %>% filter(!(Reference == "Leighton-Goodall et al. 2012" & Sample_start_year == 2003))

## MEYER ET AL 2007  ------------------------------------------------------------
#Take out Study 1
sites_mod <- sites_mod %>% filter(!(Reference == "Meyer et al. 2007" & Study_number != 1))

## SAMPLING EFFORT RECALCULATION (all studies)  ------------------------------------------------------------
#Fixing the rescaled sampling effort and the effort-corrected abundance values 
#PREDICTS states that the effort-corrected abundance values are the values corrected across sites within a Study by
#dividing the abundance measurement by sampling effort" assuming that sampled abundances increase linearly with sampling effort
#after first rescaling effort values within each Study to a maximum value of one
#We recompute this ourselves for every study  by dividing Sampling_effort 
#by the maximum Sampling_effort within each SS, since the released Rescaled_sampling_effort field shows 
#internal inconsistencies for some sites, traced to an undocumented rescaling mechanism upstream in PREDICTS.

sites_mod <- sites_mod %>%
  mutate(Sampling_effort = as.numeric(Sampling_effort)) %>%
  group_by(SS) %>%
  mutate(Max_Sampling_effort = max(Sampling_effort, na.rm = TRUE),
         Rescaled_sampling_effort_mod = Sampling_effort / Max_Sampling_effort) %>%
  ungroup() %>%
  select(-Max_Sampling_effort)

#And now we modify TA_corrected_mod to be equal to TA / Rescaled_sampling_effort_mod for every site
sites_mod <- sites_mod %>%
  mutate(TA = as.numeric(TA),
         Rescaled_sampling_effort_mod = as.numeric(Rescaled_sampling_effort_mod)) %>%
  mutate(TA_corrected_mod = TA / Rescaled_sampling_effort_mod)

## OSGATHORPE ET AL. 2012  ------------------------------------------------------------
sites_mod <- sites_mod %>%filter(!(SSBS == "AD1_2012__Osgathorpe 2  3" ))
sites_mod <- sites_mod %>%filter(!(SSBS == "AD1_2012__Osgathorpe 2  13" ))

## SMITH ET AL 2008A  ------------------------------------------------------------
#Take out the reference
sites_mod <- sites_mod %>%filter(!(Reference == "Smith et al. 2008" ))

## SMITH 2006  ------------------------------------------------------------
#Fixing problems in "Habitat sites"
smith2006 <- sites_mod %>% filter(Reference == "Smith 2006")
smith2006 <- smith2006 %>% mutate(Site_name_number = str_extract(Site_name, "\\d+$"))

#Importing the corrected coordinates given by the author. The file was created by us based on a file shared by the main author.
smith2006_coordinates_from_author <- read_delim("Databases/Smith2006_coordinates/Smith2006_coodinates_from_author.csv", delim = ";")

smith2006_coordinates_from_author <- smith2006_coordinates_from_author %>%
  rename(CODE = "Code ", SITE_TYPE = "Site type") %>%
  mutate( Longitude_mod = as.numeric(str_extract(Name, "(?<=POINT \\()[^ ]+")),
          Latitude_mod  = as.numeric(str_extract(Name, "(?<= )[0-9.-]+(?=\\))"))) %>%
  select(CODE, SITE_TYPE, HABITAT, Longitude_mod, Latitude_mod)

smith2006_coordinates_from_author <- smith2006_coordinates_from_author %>%
  mutate(CODE = as.numeric(CODE))

smith2006_study2_keys <- smith2006 %>%
  filter(Study_number == 1) %>%
  select(SSBS, Site_name_number) %>% 
  distinct()

author_coords_with_ssbs <- smith2006_coordinates_from_author %>%
  mutate(CODE_chr = as.character(CODE)) %>%
  left_join(smith2006_study2_keys, by = c("CODE_chr" = "Site_name_number")) %>% 
  filter(!is.na(SSBS)) %>%
  select(SSBS, Site_name_number = CODE_chr, Latitude_mod, Longitude_mod) %>% 
  distinct()

sites_mod <- sites_mod %>%
  left_join(author_coords_with_ssbs, by = "SSBS", suffix = c("", "_author")) %>%
  mutate(Latitude_mod  = ifelse(!is.na(Latitude_mod_author), Latitude_mod_author, Latitude_mod),
         Longitude_mod = ifelse(!is.na(Longitude_mod_author), Longitude_mod_author, Longitude_mod)) %>%
  select(-Latitude_mod_author, -Longitude_mod_author)

#Fixing problems in "Trasect-edges sites"
smith2006_edges <- sites_mod %>% filter(Reference == "Smith 2006" & Study_number == 2)

#change Site_name_number for the number at the end of the Site_name column
smith2006_edges <- smith2006_edges %>%
  mutate(Site_name_number = str_extract(Site_name, "\\d+$"))

#convert to number
smith2006_edges <- smith2006_edges %>%
  mutate(Site_name_number = as.numeric(Site_name_number))

#Importing the corrected coordinates given by the author. The file was created by us based on a file shared by the main author.
#The document sent by the author was mislabeling lat and lon so here we exchange them.
smith2006_coordinates_from_author_edges <- read_delim("Databases/Smith2006_coordinates/Smith2006_coordinates_from_author_edges.csv", delim = ";") %>%
  rename(Latitude_author  = Longitude_author,  
         Longitude_author = Latitude_author  )

#by Site_name_number, change Latitude_mod and Longitud_mod in smith2006_edges by Latitude_author and Longitud_author in smith2006_coordinates_from_author_edges
smith2006_edges <- smith2006_edges %>%
  left_join(smith2006_coordinates_from_author_edges, by = "Site_name_number" ) %>% 
  mutate(Latitude_mod  = ifelse(!is.na(Latitude_author), Latitude_author, Latitude_mod),
         Longitude_mod = ifelse(!is.na(Longitude_author), Longitude_author, Longitude_mod)) %>%
  select(-Latitude_author, -Longitude_author)

#Now we modify Latitude_mod and Longitude_mod by the values is smith2006_edges, by SSBS
sites_mod <- sites_mod %>%
  left_join(smith2006_edges %>% select(SSBS, Latitude_mod, Longitude_mod), by = "SSBS", suffix = c("", "_smith2006_edges")) %>%
  mutate(Latitude_mod  = ifelse(!is.na(Latitude_mod_smith2006_edges), Latitude_mod_smith2006_edges, Latitude_mod),
         Longitude_mod = ifelse(!is.na(Longitude_mod_smith2006_edges), Longitude_mod_smith2006_edges, Longitude_mod)) %>% 
  select(-Latitude_mod_smith2006_edges, -Longitude_mod_smith2006_edges)

## WOODCOCK ET AL. 2007  ------------------------------------------------------------
sites_mod <- sites_mod %>% filter(!(Reference == "Woodcock et al. 2007" & !grepl("7", Site_name)))

## ADDITIONAL REMOVE OF SOURCES DUE TO LAND USE CATEGORIES DEFINITION  ------------------------------------------------------------
#There are additional references that is necesary to take out because when we performe lated tha recategorization of habitats (leaving out some), those reference will not have more than one land use-use intensity to compare
sources_to_delete <- c("Verboven et al. 2012",  "Leighton-Goodall et al. 2012", "Noreika 2009",  "Nielsen et al. 2011", "Gaublomme et al. 2008" , "Magura et al. 2010")
sites_mod <- sites_mod %>%filter(!Reference%in% sources_to_delete)

### RECATGORIZATION OF LAND USE------------------------------------------------------------
#creating a new column for the recategorizations, without loosing the original one
sites_mod <- sites_mod %>% 
  mutate (Predominant_land_use_mod = Predominant_land_use)

#Habitats recategorization
sites_mod <- sites_mod %>%
  mutate (Predominant_land_use_mod = ifelse (Predominant_land_use_mod=="Primary vegetation", "Natural vegetation", Predominant_land_use_mod))
sites_mod <- sites_mod %>%
  mutate (Predominant_land_use_mod = ifelse (Predominant_land_use_mod=="Mature secondary vegetation", "Natural vegetation", Predominant_land_use_mod))
sites_mod <- sites_mod %>%
  mutate (Predominant_land_use_mod = ifelse (Predominant_land_use_mod=="Young secondary vegetation", "Semi-natural vegetation", Predominant_land_use_mod))
sites_mod <- sites_mod %>%
  mutate (Predominant_land_use_mod = ifelse (Predominant_land_use_mod=="Intermediate secondary vegetation", "Semi-natural vegetation", Predominant_land_use_mod))
sites_mod <- sites_mod %>%
  mutate (Predominant_land_use_mod = ifelse (Predominant_land_use_mod=="Secondary vegetation (indeterminate age)", "Semi-natural vegetation", Predominant_land_use_mod))

#Leaving out none pasture, cropland or semi-natural vegetation categories
sites_mod <- sites_mod %>% filter(!Predominant_land_use_mod %in%  c("Urban","Plantation forest", "Cannot decide"))


##FILTERING OBSERVATION LEVEL DATASET------------------------------------------------------------
#Producing a new dataframe of observations filtered by the same changes we did to sites
observations_mod <- observations_df %>% filter(SSBS %in% sites_mod$SSBS)

#Creating a new coordinate ID given latitud and longitud where in some cases changed.
sites_mod <- sites_mod %>% 
  mutate(Lat_Long_mod = paste(Latitude_mod, Longitude_mod, sep = ","))

#Adding new columns to the observation dataset
observations_mod <- observations_mod %>% left_join(sites_mod %>% select(SSBS, Longitude_mod, Latitude_mod, Lat_Long_mod), by = "SSBS")

#ESTIMATING TAXA AND GENUS RICHNESS BY SITE-----------------------------------------------------------
#estimate taxa richness by site
taxa_richness_by_site <- observations_mod %>%
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
sites_mod <- sites_mod %>%
  left_join(taxa_richness_by_site %>% select(SSBS, taxa_richness), by = "SSBS")

#estimate genus richness by site. 
genus_richness <- observations_mod %>%
  filter(Measurement > 0) %>%
  group_by(SSBS) %>%
  summarise(genus_richness= n_distinct(Genus), .groups = "drop")

sites_mod <- sites_mod %>%
  left_join(genus_richness, by = "SSBS")

#SAVE  ------------------------------------------------------------
write.csv(sites_mod, "Intermediate_dataset/sites_mod.csv",row.names = FALSE)
write.csv(observations_mod, "Intermediate_dataset/observations_mod.csv",row.names = FALSE)