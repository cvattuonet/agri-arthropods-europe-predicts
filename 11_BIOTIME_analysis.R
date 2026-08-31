#---------------------------------------------------------------------------------------------------------------------------------------#
# CODE DETAILS                                                                                                                          #
#                                                                                                                                       #
# Author: Catalina Vattuone  (cvattuonet@gmail.com)                                                                                     #        
# Date of latest update: 31-08-2026                                                                                                     #
# Main use: Study on "Assessing how agricultural practices and landscape composition shape differences in arthropod biodiversity among  #  
#           habitat types in Europe using the PREDICTS database". Internship at INRAE                                                   #
#                                                                                                                                       #
# Content: Comparison of BioTIME dataset with PREDICTS dataset, and production of combined dataset                                      #                           
#---------------------------------------------------------------------------------------------------------------------------------------#


#LIBRARY ------------------------------------------------------------
pacman::p_load(DBI, RMariaDB, BioTIMEr, RColorBrewer, taxize, rgbif, dplyr, purrr, tibble, sf, rnaturalearth, rnaturalearthdata, ggplot2,  stringr, scales, ggbeeswarm, readr, tidyr)

#STEP 1: DATA LOADING BIOTIME------------------------------------------------------------
#Downloding link https://biotime.st-andrews.ac.uk/getFullDownload.php
#I started by using the rds and csv files that can be downlodad in the same link
#BUT they don't come with all the available information, that is in the SQL file. 
#So, after installing in my terminal mySQL, and downlowing the SQL file available at the download link of the BioTIME database
#I extract all the tables that come in the SQL by using the following code:

#first I access the data
# con <- dbConnect(RMariaDB::MariaDB(), dbname = "biotime", host = "localhost", user = "root", password = "")
# table_names <- dbListTables(con) #here I list all the tables to later extract
# 
# for (tbl in table_names) { #here I extract every table
#   message("Reading: ", tbl)
#   df <- dbReadTable(con, tbl)
#   saveRDS(df, paste0("Databases/BIOTIME_v2/biotime_v2_", tbl, "_2025.rds"))
#   rm(df)
#   gc()
# }
# dbDisconnect(con)
# 
# #checking that it worked
# list.files("Databases/BIOTIME_v2/", pattern = "\\.rds$") #it prints all the created files

#now I upload all the tables that I will need fro the analysis without having to do the SQL extraction again
biotime_obs <- readRDS("Databases/BIOTIME_v2/biotime_v2_allrawdata_2025.rds")
species_table <- readRDS("Databases/BIOTIME_v2/biotime_v2_species_2025.rds")
methods_table <- readRDS("Databases/BIOTIME_v2/biotime_v2_methods_2025.rds")
citation_table <- readRDS("Databases/BIOTIME_v2/biotime_v2_citation1_2025.rds")
datasets_table <- readRDS("Databases/BIOTIME_v2/biotime_v2_datasets_2025.rds")

#here I check whether the keys that I will use to join the tables are unique, because if they are not unique, 
#then the join will create duplicates and the number of rows will increase.
species_table %>% count(newID) %>% filter(n > 1) %>% nrow()      
methods_table %>% count(STUDY_ID) %>% filter(n > 1) %>% nrow()   
citation_table%>% count(STUDY_ID) %>% filter(n > 1) %>% nrow()  
datasets_table %>% count(STUDY_ID) %>% filter(n > 1) %>% nrow()  

#in the case of the citation I obtained 205 so I need to collaps multiple citations per study
citation_table_collapsed <- citation_table %>%
  group_by(STUDY_ID) %>%
  summarise(CITATION_ID = paste(unique(CITATION_ID), collapse = ", "),
            BIB = paste(unique(BIB), collapse = " ;; "),
            .groups = "drop")

#now I join all the tables to the main one
biotime_obs <- biotime_obs %>%
  left_join(species_table %>% select(newID, GENUS, family, order, class, phylum, kingdom, valid_name, resolution, taxon), by = "newID") %>%
  left_join(methods_table %>% select(STUDY_ID, METHODS, SUMMARY_METHODS), by = "STUDY_ID") %>%
  left_join(citation_table_collapsed, by = "STUDY_ID") %>%
  left_join(datasets_table %>% select(STUDY_ID, TITLE, TAXA, ORGANISMS, AB_BIO), by = "STUDY_ID")

#in biotime_obs rename TITLE for TITLE_dataset
biotime_obs <- biotime_obs %>%
  rename(TITLE_dataset = TITLE)

#STEP 1b: DATA LOADING PREDICTS ------------------------------------------------------------
sites_predicts<- read_delim("Intermediate_dataset/abundance_sites_for_vizualization_nonscale.csv")
n_distinct(sites_predicts$Coordinate_ID)
observations_predicts <- read_delim("Intermediate_dataset/observations_mod.csv")
observations_predicts <- observations_predicts %>%
  filter(SSBS %in% sites_predicts$SSBS)

sites_full_info <- read_delim("Intermediate_dataset/sites_mod.csv")  #there is info on the sites that is not available in the abundance_sites_for_vizualization_nonscale file because 
#it was filtered to only include only information relevant for the model 


#STEP 2: OBSERVATIONS ARE FILTERED  ------------------------------------------------------------
biotime_obs <- biotime_obs %>%
  filter(taxon %in% c("Terrestrial/freshwater invertebrates", "Invertebrates", "Marine/freshwater invertebrates", "Mixed")) %>%
  filter(phylum == "Arthropoda") %>%
  filter(class %in% c("Insecta", "Arachnida")) 

#now we start filtering by location, first by a rought filter because the data does not come with the names of the countries
#adding the countries takes time so is better is we cut it like that first
biotime_obs <- biotime_obs %>%
  filter(LONGITUDE > -10 & LATITUDE > 35) 

#now we want to filter by the countries in the EU and UK so we add the country name first 
biotime_sf <- st_as_sf(biotime_obs, coords = c("LONGITUDE", "LATITUDE"), crs = 4326, remove = FALSE)
world <- ne_countries(scale = "medium", returnclass = "sf")
biotime_obs <- st_join(biotime_sf, world[, c("name")], join = st_intersects) %>%
  rename(Country = name) 

#we take out the points in the ocean
biotime_obs <- biotime_obs %>%
  filter(!is.na(Country)) 

#filter by eu - uk countires
unique(biotime_obs$Country) #check what countries are availble
eu_uk_countries <- c("United Kingdom", "Italy", "Germany", "Netherlands", "Finland", "Latvia", "Belgium", #this is not all eu_uk countries but the ones that are in the dat
                     "Romania","Sweden","Spain","Hungary","Austria", "France", "Denmark","Czechia")

biotime_obs <- biotime_obs %>%
  filter(Country %in% eu_uk_countries)


#STEP 3: SSBS CODE ------------------------------------------------------------
##Source_ID ------------------------------------------------------------
#the citation comes in a format that is hard to read so we have to transform it 
#here are some functions to help extract the field of the authors, 
#clean it, get the first lastname of the authors, the year of the publication and the title of the publication 

extract_author_field <- function(bib_string) {
  start_pos <- regexpr("Author\\s*=\\s*\\{", bib_string)
  if (start_pos == -1) return(NA_character_)
  content_start <- start_pos + attr(start_pos, "match.length")
  chars <- strsplit(substring(bib_string, content_start), "")[[1]]
  depth <- 1
  end_idx <- NA
  for (i in seq_along(chars)) { #list of authors is 
    if (chars[i] == "{") depth <- depth + 1
    if (chars[i] == "}") depth <- depth - 1
    if (depth == 0) { end_idx <- i; break }
  }
  if (is.na(end_idx)) return(NA_character_)
  paste(chars[1:(end_idx - 1)], collapse = "")
}

clean_latex <- function(x) {
  x <- gsub("\\{\\\\[\"'`^~vc.]\\s*([a-zA-Z])\\}", "\\1", x)
  x <- gsub("\\{\\\\([a-zA-Z]+)\\}", "\\1", x)
  x <- gsub("[{}]", "", x)
  trimws(x)
}

#here we divide in many cases because the authors data is reported differently in different cases. 
#Extracting the first last name or organization name is no so easy
get_first_author_lastname <- function(author_field) {
  if (is.na(author_field)) return(NA_character_)
  trimmed <- trimws(author_field)
  
  #ORGANIZATION CASE: whole field wrapped in an extra brace pair, e.g. {National Biodiversity Network}
  if (grepl("^\\{.*\\}$", trimmed)) {
    return(clean_latex(gsub("^\\{|\\}$", "", trimmed)))
  }
  
  #PERSON CASE: split on " and ", take the first author
  first_author <- trimws(strsplit(trimmed, "\\s+and\\s+")[[1]][1])
  if (grepl(",", first_author)) {
    lastname <- trimws(strsplit(first_author, ",")[[1]][1])  #"Lastname, Firstname" format
  } else {
    words <- strsplit(first_author, "\\s+")[[1]]
    lastname <- words[length(words)]  #"Firstname Lastname" format -> take last word
  }
  clean_latex(lastname)
}

#year is also reported within the bibtex string so we extract it with the following function
extract_year <- function(bib_string) {
  m <- regmatches(bib_string, regexpr("Year\\s*=\\s*\\{(\\d{4})\\}", bib_string))
  if (length(m) == 0) return(NA_character_)
  gsub(".*\\{(\\d{4})\\}.*", "\\1", m)
}

#same for year
extract_title_field <- function(bib_string) {
  starts <- gregexpr("Title\\s*=\\s*\\{", bib_string)[[1]]
  if (starts[1] == -1) return(character(0))
  match_lens <- attr(starts, "match.length")
  
  extract_one <- function(start_pos, match_len) {
    content_start <- start_pos + match_len
    chars <- strsplit(substring(bib_string, content_start), "")[[1]]
    depth <- 1
    end_idx <- NA
    for (i in seq_along(chars)) {
      if (chars[i] == "{") depth <- depth + 1
      if (chars[i] == "}") depth <- depth - 1
      if (depth == 0) { end_idx <- i; break }
    }
    if (is.na(end_idx)) return(NA_character_)
    paste(chars[1:(end_idx - 1)], collapse = "")
  }
  
  titles <- mapply(extract_one, starts, match_lens)
  clean_latex(titles[!is.na(titles) & titles != ""])
}

citation_table <- citation_table %>%
  mutate(Year_publication = sapply(BIB, extract_year),
         First_author = sapply(BIB, function(b) get_first_author_lastname(extract_author_field(b))))

citation_table <- citation_table %>%
  mutate(Publication_title = sapply(BIB, function(b) {
    t <- extract_title_field(b)
    if (length(t) == 0) NA_character_ else t[1]   # first citation's title, same convention as First_author/Year_publication
  }))

citation_table_collapsed <- citation_table %>%
  group_by(STUDY_ID) %>%
  summarise(CITATION_ID = paste(unique(CITATION_ID), collapse = ", "),
            BIB = paste(unique(BIB), collapse = " ;; "),
            Year_publication = first(Year_publication), #I use first because sometimes there is multiple citations for one study
            First_author = first(First_author),
            Publication_title = first(Publication_title),
            .groups = "drop") %>%
  mutate(Source_ID = paste(Year_publication, First_author, sep = "_"))

#now I add it to my dataframe
biotime_obs <- biotime_obs %>%
  left_join(citation_table_collapsed %>% select(STUDY_ID, Source_ID, Year_publication, First_author, Publication_title), by = "STUDY_ID")

unique(biotime_obs$Source_ID) #30 studies with 2 weird values:  NA_Langenhaun" and NA, The first one correspond does not have a year in the raw citation data

#but looking into the DOI I found it was published in 2020. So it is hereby corrected
## LANGENHAUN STUDIES (STUDY_ID 625, 626) - BIB missing Year field; confirmed via DOI (10.18728/536.0) that publication year is 2020

citation_table_collapsed <- citation_table_collapsed %>%
  mutate(Year_publication = ifelse(STUDY_ID %in% c(625, 626), "2020", Year_publication),
    Source_ID = ifelse(STUDY_ID %in% c(625, 626), paste(Year_publication, First_author, sep = "_"), Source_ID))

biotime_obs <- biotime_obs %>%
  select(-Source_ID, -Year_publication, -First_author, -Publication_title) %>%
  left_join(citation_table_collapsed %>% select(STUDY_ID, Source_ID, Year_publication, First_author, Publication_title), by = "STUDY_ID")

#now I analyse the missing study which oes have a title but not a citation
missing_citation_study <- biotime_obs %>% filter(is.na(BIB)) %>% distinct(STUDY_ID, TITLE_dataset)
missing_citation_study #in fact there are 3 cases where no citation is available

#I will change in the contact table
contacts_table <- readRDS("Databases/BIOTIME_v2/biotime_v2_contacts_2025.rds")
contacts_table %>% filter(STUDY_ID %in% missing_citation_study$STUDY_ID) %>% select(STUDY_ID, DATA_SOURCE, WEB_LINK) #it seems it non published data, so nothing to cite

#I will give them a Source_ID of unpublish data to avoid problems, using the code 
citation_table_collapsed <- citation_table_collapsed %>%
  bind_rows(tibble(STUDY_ID = c(454, 683, 688),Source_ID = paste0("Unpublished_STUDY", c(454, 683, 688))))

biotime_obs <- biotime_obs %>%
  select(-Source_ID, -Year_publication, -First_author, -Publication_title) %>%
  left_join(citation_table_collapsed %>% select(STUDY_ID, Source_ID, Year_publication, First_author, Publication_title), by = "STUDY_ID")
unique(biotime_obs$Source_ID)


##Block ------------------------------------------------------------
#In Biotime there is no blocks so we assign all the same value
biotime_obs <- biotime_obs %>%
  mutate(Block = 1)

##Study_number ------------------------------------------------------------
#In PREDICTS we give different study numbers to different sampling methods withtin a reference so we inspect if each reference has a different sampling method
biotime_obs <- st_drop_geometry(biotime_obs)
class(biotime_obs)

biotime_obs %>% #yes, in fact 2 cases have more methods, but the problem could be just synthaxis so I will look deeper at those two
  distinct(Source_ID, STUDY_ID, SUMMARY_METHODS) %>%
  group_by(Source_ID) %>%
  summarise(n_distinct_studyid = n_distinct(STUDY_ID),
            n_distinct_methods = n_distinct(SUMMARY_METHODS),
            methods = paste(unique(SUMMARY_METHODS), collapse = " | "),
            .groups = "drop") %>%
  arrange(desc(n_distinct_methods))


biotime_obs <- biotime_obs %>%
  group_by(Source_ID) %>%
  mutate(Study_number = dense_rank(SUMMARY_METHODS)) %>%
  ungroup()

unique(biotime_obs$Study_number)

##Site number ------------------------------------------------------------
biotime_obs <- biotime_obs %>%
  mutate(Coordinate_ID = paste(LATITUDE, LONGITUDE, sep = "_"))

biotime_obs <- biotime_obs %>%
  group_by(Source_ID, Study_number) %>%
  mutate(Site_number = match(Coordinate_ID, unique(Coordinate_ID))) %>%
  ungroup()

##Concatenation ------------------------------------------------------------
biotime_obs <- biotime_obs %>%
  mutate(SS   = paste(Source_ID, Study_number, sep = "_"),
    SSB  = paste(SS, Block, sep = "_"),
    SSBS = paste(SSB, Site_number, sep = "_") )

#STEP 4: Sampling methods   ------------------------------------------------------------
sampling_methods_predicts <- unique(observations_predicts$Sampling_method) #all unique sampling methods in predicts dataset for europe
sampling_methods_biotime <- unique(biotime_obs$SUMMARY_METHODS) #all unique sampling methods in biotime

#now changing according to Appendix 1 of the Note
#first predicts
predicts_lookup <- tibble::tribble(
  ~Sampling_method,                ~standardized_sampling_method,
  "specimen collection",           "specimen collection",
  "sweep nets",                    "sweep netting",
  "sweep netting",                 "sweep netting",
  "sweep path",                    "sweep netting",
  "flower visitation observation", "flower visitation observation",
  "flight trap",                   "flight trap",
  "light trap",                    "light trap",
  "pit-fall traps",                "pit-fall traps",
  "baited pit-fall traps",         "baited pit-fall traps",
  "systematic searching",          "systematic searching",
  "pan traps",                     "pan traps",
  "soil core",                     "soil core",
  "aerial transect",               "aerial transect",
  "line/belt transects",           "transect",
  "transect",                      "transect",
  "window traps",                  "window traps",
  "window trap",                   "window traps",
  "suction samplers",              "suction samplers",
  "visual encounter survey",       "visual encounter survey",
  "various",                       "various",
  "fixed plots/quadrats",          "fixed plots/quadrats"
)

#then BioTIME
biotime_lookup <- tibble::tribble(
  ~SUMMARY_METHODS, ~standardized_sampling_method,
  
  "Transects", "transect",
  "Plots", "fixed plots/quadrats",
  "10 x 10 m squares", "fixed plots/quadrats",
  "Nets", "sweep netting",
  "Organisms were collected through Surber samplings (0.0506 squared metres, mesh size 500 µm) in 3 sites at different elevation", "surber sampling",
  "Using kicksampling 20 subsamples corresponding to a total of 1.25 m2 of stream bed were collected at each site and pooled into a bulk sample, with subsamples distributed among the major habitats proportionally to their coverage within the site.", "surber sampling",
  "Pitfall traps, annual sums.", "pit-fall traps",
  "pitfall traps", "pit-fall traps",
  "Pitfal traps", "pit-fall traps",
  "Pitfall traps, abundance as annual activity density", "pit-fall traps",
  "Pitfall traps installed in UTM grid cells (1 km x 1 km) all over the NPHK. Traps left for appr. 1 month, emptied & refilled for the spring-summer period (not all sites every year, but all habitats every year).", "pit-fall traps",
  "Pitfall traps established along transects and emptied and replaced fortnightly.", "pit-fall traps",
  "For 9 years, each site was monitored using a sampling grid consisting of 5 parallel transects, with 20 soil samples collected from each transect.", "soil core",
  "Tidal flats are sampled for macroinvertebrates by core sampling 4.5 cm diam. Since 2008 points are random stratified sampled in the system. Before they were on fixed transects. Taxa determined on species level (expect for Oligochaeta)", "soil core",
  "Malaise trap, weekly counts from April to October.", "malaise trap",
  "Standardized sweep netting.", "sweep netting",
  "Mist nets were used to sample bird populations. Ten mist nets were operated on 1-3 occasions each autumn (October-December), each lasting for 2-5 days.", "mist net",
  "butterfly traps", "baited trap",
  "A standard Rothamsted light trap (Williams 1948) is used to trap moths. Traps are switched on at dusk and off at dawn each day by an automatic solar dial time-switch which self-adjusts for seasonal variation in daylength. Traps are emptied daily.", "light trap",
  "Light traps - daily counts over yearly active flight period of moths. Site-year data was left out if: 1) ID was restricted to pest species, 2) there were >10 days of no records. Eupithecia spp. were excluded because difficult to ID.", "light trap",
  "Weekly light trap collections in Finnish forests deployed in April to October.", "light trap",
  "captures per minute of sampling effort", "unclear",
  "Weekly or biweekly sampling throughout water column at various sampling stations.", "unclear",
  "Adult damselfies were surveyed from multiple ponds across NE scotland, yearly between May and August. Densities were caculated as number of damselflies caught divided by the total survey time multiplied by number of surveyors.", "unclear",
  "Densities of bird populations were estimated by territory mapping and nest counts in study plots. Results were expressed in territories per km2.", "unclear",
  "Saprobic system DIN 38410 before 2000, multi-habitat method following the German macroinvertebrate monitoring protocol for the European Waterframework Directive after 2000. Data after 2000 converted to the 6 abundance classes of the saprobic system.", "surber sampling",
  "Sampling procedure followed the saprobic system DIN 38410 for all samples taken before year 2000 and Multi-habitat sampling method according to the German  macroinvertebrate monitoring standard protocol for the European Waterframework Directive for s", "surber sampling",
  "Standardised butterfly counts along fixed transects. Transects were walked weekly between 1 April and 29 September, providing the meteorological conditions were met.", "transect",
  "Communities of saproxylic beetles and spiders were sampled for 3 years on deadwood in 6 sites within Steigerwald forest, using stem emergence traps.", "emergence trap",
  "Saproxylic beetles were trapped using flight-interception traps placed at the plot centroid", "flight trap",
  "Light traps", "light trap",
  "Counts within defined area", "fixed plots/quadrats",
  "Stratified sampling of vascular plants in plots at different elevations.", "unclear",
  "Annual (April/May) point count bird assemblage monitoring along three French rivers from 1999 to 2011.", "unclear",
  "Volunteer ornithologists counted birds annually along line transects and point counts in 2 x 2 km squares across Sweden. There are 716 fixed routes. Surveys are conducted during breeding season at 04:00 but no earlier than 30 minutes before sunrise.", "unclear",
  "Sampling was done 11 times from 1921 to 2011 in 5 permanent 5-m radius circle plots. All vascular plants were determined in the 5 subplots and the cover was estimated for each species.", "fixed plots/quadrats"
)

#trim whitespace on both sides before joining
observations_predicts <- observations_predicts %>% mutate(Sampling_method = str_squish(Sampling_method))
predicts_lookup <- predicts_lookup  %>% mutate(Sampling_method = str_squish(Sampling_method))

biotime_obs <- biotime_obs  %>% mutate(SUMMARY_METHODS = str_squish(SUMMARY_METHODS))
biotime_lookup <- biotime_lookup %>% mutate(SUMMARY_METHODS = str_squish(SUMMARY_METHODS))

#join the standardized category onto each dataframe
observations_predicts <- observations_predicts %>%left_join(predicts_lookup, by = "Sampling_method")
biotime_obs <- biotime_obs %>%left_join(biotime_lookup, by = "SUMMARY_METHODS")

#check
observations_predicts %>% filter(is.na(standardized_sampling_method)) %>% distinct(Sampling_method)
biotime_obs %>% filter(is.na(standardized_sampling_method)) %>% distinct(SUMMARY_METHODS)

#STEP 5: Measurement vs Abundance and Biomass  & Coordinate ID ------------------------------------------------------------
biotime_obs <- biotime_obs %>%
  pivot_longer( cols = c(ABUNDANCE, BIOMASS),
                names_to = "Diversity_metric",
                values_to = "Measurement",
                values_drop_na = TRUE ) %>%
  mutate(Diversity_metric = recode(Diversity_metric,"ABUNDANCE" = "abundance","BIOMASS" = "biomass"))

#coordinate ID 
biotime_obs <- biotime_obs %>%
  mutate(Coordinate_ID = paste(LATITUDE, LONGITUDE, sep = ","))

#STEP 6: Site Level Database   ------------------------------------------------------------
sites_biotime <- biotime_obs %>%
  group_by(SSBS, YEAR) %>%
  summarise(Source_ID = paste(unique(Source_ID), collapse = ", "),
            Title_dataset = paste(unique(TITLE_dataset), collapse = ", "),
            Publication_title = paste(unique(Publication_title), collapse = ", "),
            Study_number = paste(unique(Study_number), collapse = ", "),
            METHODS = paste(unique(METHODS), collapse = ", "),
            SUMMARY_METHODS = paste(unique(SUMMARY_METHODS), collapse = ", "),
            standardized_sampling_method = paste(unique(standardized_sampling_method), collapse = ", "),
            Site_number = paste(unique(Site_number), collapse = ", "),
            Block = paste(unique(Block), collapse = ", "),
            SS = paste(unique(SS), collapse = ", "),
            SSB = paste(unique(SSB), collapse = ", "),
            LONGITUDE= paste(unique(LONGITUDE), collapse = ", "),
            LATITUDE= paste(unique(LATITUDE), collapse = ", "),
            Coordinate_ID = paste(unique(Coordinate_ID), collapse = ", "),
            Country = paste(unique(Country), collapse = ", "),
            resolution = paste(unique(resolution[Measurement > 0 & !is.na(Measurement) & !is.na(resolution) & resolution != ""]), collapse = ", "),
            Diversity_metric = paste(unique(Diversity_metric), collapse = ", "),
            N_class = n_distinct(class[Measurement > 0 & !is.na(Measurement) & !is.na(class) & class != ""]),
            class = paste(unique(class[Measurement > 0 & !is.na(Measurement) & !is.na(class) & class != ""]), collapse = ", "),
            N_order = n_distinct(order[Measurement > 0 & !is.na(Measurement) & !is.na(order) & order != ""]),
            order = paste(unique(order[Measurement > 0 & !is.na(Measurement) & !is.na(order) & order != ""]), collapse = ", "),
            N_family = n_distinct(family[Measurement > 0 & !is.na(Measurement) & !is.na(family) & family != ""]),
            family = paste(unique(family[Measurement > 0 & !is.na(Measurement) & !is.na(family) & family != ""]), collapse = ", "),
            N_GENUS = n_distinct(GENUS[Measurement > 0 & !is.na(Measurement) & !is.na(GENUS) & GENUS != ""]),
            GENUS = paste(unique(GENUS[Measurement > 0 & !is.na(Measurement) & !is.na(GENUS) & GENUS != ""]), collapse = ", "),
            valid_name = paste(unique(valid_name[Measurement > 0 & !is.na(Measurement) & !is.na(valid_name) & valid_name != ""]), collapse = ", "),
            TA=sum(Measurement, na.rm=TRUE),
            .groups = "drop")

taxa_richness_by_site_biotime <- biotime_obs %>%
  filter(Measurement > 0) %>%
  group_by(SSBS, YEAR) %>%
  summarise(taxa_richness = n_distinct(valid_name), .groups = "drop")

sites_biotime <- sites_biotime %>%
  left_join(taxa_richness_by_site_biotime %>% select(SSBS, YEAR, taxa_richness),by = c("SSBS", "YEAR"))


#ANALYSING BIOTIME -----------------------------------------------------------
#amount of observations
nrow(biotime_obs)
nrow(observations_predicts)

#amount of unique sites
n_distinct(sites_biotime$Coordinate_ID)
n_distinct(sites_predicts$Coordinate_ID)

#number of distinct SSBS
n_distinct(sites_biotime$SSBS)
n_distinct(sites_predicts$SSBS)

#amount of observations from particular orders
hemiptera <- biotime_obs %>%
  filter(order == "Hemiptera") 
hymenoptera <- biotime_obs %>%
  filter(order == "Hymenoptera") 

#Count the amount of different years with sampling points for each SSBS
sites_biotime_summary <- sites_biotime %>%
  group_by(SSBS) %>%
  summarise(n_years = n_distinct(YEAR), .groups = "drop")
#mean, min and max n_years
summary(sites_biotime_summary$n_years)

#count the amount of SSBS with n_year = 1
sum(sites_biotime_summary$n_years == 1)

#min YEAR and max year in the data
min(sites_biotime$YEAR, na.rm = TRUE)
max(sites_biotime$YEAR, na.rm = TRUE)

#sites in Sweden
sites_biotime_sweden <- sites_biotime%>%
  filter(Country == "Sweden")

#sites in France
sites_biotime_france <- sites_biotime%>%
  filter(Country == "France")

#country summary
unique(sites_biotime$Country)
unique(sites_predicts$Country)

#countries that are in sites_biotime but not in sites_predicts?
setdiff(unique(sites_biotime$Country), unique(sites_predicts$Country))

#looking at the taxa of observations from each country in the difference
sites_biotime_austria <- sites_biotime%>%
  filter(Country == "Austria")
sites_biotime_latvia <- sites_biotime%>%
  filter(Country == "Latvia")

#opposite direction?
setdiff(unique(sites_predicts$Country), unique(sites_biotime$Country))

#looking at the taxa of observations from each country in the difference
sites_biotime_estonia <- observations_predicts%>%
  filter(Country == "Estonia")
unique(sites_biotime_estonia$Class)
unique(sites_biotime_estonia$Order)
sites_biotime_poland <- observations_predicts%>%
  filter(Country == "Poland")
unique(sites_biotime_poland$Class)
unique(sites_biotime_poland$Order)
sites_biotime_greece <- observations_predicts%>%
  filter(Country == "Greece")
unique(sites_biotime_greece$Class)
unique(sites_biotime_greece$Order)
sites_biotime_Portugal <- observations_predicts%>%
  filter(Country == "Portugal")
unique(sites_biotime_Portugal$Class)
unique(sites_biotime_Portugal$Order)
sites_biotime_ireland <- observations_predicts%>%
  filter(Country == "Ireland")
unique(sites_biotime_ireland$Class)
unique(sites_biotime_ireland$Order)

#summaty table with amount of individual SSBS and percentage, by country in sites_biotime
summary_biotime_country <- sites_biotime %>%
  group_by(Country) %>%
  summarise(n_unique_coordinate = n_distinct(Coordinate_ID), .groups = "drop") %>%
  mutate(percentage = (n_unique_coordinate / sum(n_unique_coordinate)) * 100) %>%
  arrange(desc(n_unique_coordinate))
write.csv(summary_biotime_country, "Intermediate_dataset/summary_biotime_country.csv", row.names = FALSE)

#analysing sources
unique(biotime_obs$Source_ID)

#looking at sources of interest
vanKlink <- biotime_obs %>% #because is the author of Insect Change and want to check weather the data is a particular study of him or more data from the database
  filter(Source_ID == "2021_van Klink") 
Koivula <- biotime_obs %>% #because it is the same author as one in PREDICTS
  filter(Source_ID == "2019_Koivula") 

#looking at sources titles
Dapporto <- biotime_obs %>%
  filter(Source_ID == "2009_Dapporto") 
print(Dapporto$Publication_title)
Pilotto <- biotime_obs %>%
  filter(Source_ID == "2020_Pilotto") 
print(Pilotto$Publication_title)
Lindstrom <- biotime_obs %>%
  filter(Source_ID == "2020_Lindstrom") 
print(Lindstrom$Publication_title)
Gaget <- biotime_obs %>%
  filter(Source_ID == "2020_Gaget") 
print(Gaget$Publication_title)


#PLOTS ------------------------------------------------------------
##Spatial distribution ------------------------------------------------------------
coords_predicts <- sites_predicts %>%
  mutate(SSBS = as.character(SSBS)) %>%
  left_join(observations_predicts %>% mutate(SSBS = as.character(SSBS)) %>% 
              distinct(SSBS, Latitude, Longitude),by = "SSBS") %>%
  filter(!is.na(Latitude) & !is.na(Longitude)) %>%
  distinct(Longitude, Latitude) %>%
  mutate(Dataset = "PREDICTS")

coords_biotime <- sites_biotime %>%
  distinct(Coordinate_ID, LATITUDE, LONGITUDE) %>%
  mutate(LATITUDE = as.numeric(LATITUDE),
         LONGITUDE = as.numeric(LONGITUDE)) %>%
  rename(Latitude = LATITUDE, Longitude = LONGITUDE) %>%
  select(Longitude, Latitude) %>%
  mutate(Dataset = "BioTIME")

coords_combined <- bind_rows(coords_predicts, coords_biotime)
europe_map <- ne_countries(scale = "medium", continent = "Europe", returnclass = "sf")

ggplot() +
  geom_sf(data = europe_map, fill = "grey95", color = "grey60") +
  geom_point(data = coords_combined, aes(x = Longitude, y = Latitude, color = Dataset),
             size = 1, alpha = 0.6) +
  scale_color_manual(values = c("PREDICTS" = "#0284C7", "BioTIME" = "firebrick")) +
  coord_sf(xlim = c(-25, 45), ylim = c(34, 72)) +
  theme_minimal() +
  labs(title = paste0("PREDICTS (n=", nrow(coords_predicts), ") vs. BioTIME (n=", nrow(coords_biotime), ") sites in Europe"),
       x = "Longitude", y = "Latitude", color = "Dataset")

##Taxonomic distribution ------------------------------------------------------------
taxa_clean_data_biotime <- biotime_obs %>%
  filter(Measurement > 0) %>%
  filter(!is.na(order), !is.na(family)) %>%
  group_by(order, family) %>%
  tally(name = "obs_count") %>%
  ungroup() %>%
  mutate(percentage = (obs_count / sum(obs_count)) * 100) %>%   # percentage share of total records
  filter(percentage >= 0.5) %>%                                 # only families representing >= 0.5%, to keep it clean
  group_by(order) %>%
  mutate(order_total = sum(percentage)) %>%
  ungroup() %>%
  mutate(Family = reorder(family, percentage))

plot_taxa_grid_biotime <- ggplot(taxa_clean_data_biotime, aes(x = percentage, y = family, fill = order)) +
  geom_col(width = 0.8, alpha = 0.9, show.legend = FALSE) +
  geom_text(aes(label = paste0(round(percentage, 1), "%")), hjust = -0.1, size = 3, color = "#334155") +
  facet_wrap(~ order, scales = "free_y", ncol = 2) +
  scale_x_continuous(labels = label_percent(scale = 1), expand = expansion(mult = c(0, 0.15))) +
  scale_fill_brewer(palette = "Set2") +
  theme_bw(base_size = 11) +
  labs(title = "Taxonomic composition of BioTIME dataset",
       subtitle = "Blocks represent major Orders (>1% of total share); bars show individual Family percentage shares",
       x = "Percentage of total records", y = "") +
  theme(plot.title = element_text(face = "bold", size = 13, color = "#0F172A"),
        plot.subtitle = element_text(size = 10, color = "#475569"),
        strip.text = element_text(face = "bold", size = 11, color = "#0F172A"),
        strip.background = element_rect(fill = "#F1F5F9", color = "#CBD5E1"),
        axis.text.y = element_text(size = 9, color = "#1E293B"),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(1, "lines"))

print(plot_taxa_grid_biotime)


##Temporal distribution ------------------------------------------------------------
biotime_timeline <- biotime_obs %>%
  select(Source_ID, Coordinate_ID, YEAR) %>%
  distinct()

reference_metadata_biotime <- biotime_timeline %>%
  group_by(Source_ID) %>%
  summarise(start_yr = min(YEAR, na.rm = TRUE),
            end_yr   = max(YEAR, na.rm = TRUE), .groups = "drop")

timeline_data_biotime <- biotime_timeline %>%
  left_join(reference_metadata_biotime, by = "Source_ID") %>%
  mutate(Ref_Label = reorder(Source_ID, start_yr))

plot_event_timeline_biotime <- ggplot(timeline_data_biotime) +
  geom_hline(yintercept = seq_along(levels(timeline_data_biotime$Ref_Label)),
             color = "#E2E8F0", size = 0.4, linetype = "dashed") +
  geom_segment(data = timeline_data_biotime %>% select(Ref_Label, start_yr, end_yr) %>% distinct(),
               aes(x = start_yr, xend = end_yr, y = Ref_Label, yend = Ref_Label),
               color = "#00A3A5", size = 1.0, alpha = 0.7) +
  geom_point(aes(x = YEAR, y = Ref_Label), color = "#00A3A5", size = 2, alpha = 0.9) +
  scale_x_continuous(breaks = seq(min(timeline_data_biotime$start_yr), max(timeline_data_biotime$end_yr, na.rm = TRUE), by = 5)) +
  theme_bw(base_size = 10) +
  labs(title = "Timeline covered by reference (BioTIME)",
       subtitle = "Lines show study duration; dots indicate active sampling years",
       x = "Year", y = "") +
  theme(plot.title = element_text(face = "bold", size = 13, color = "#0F172A"),
        plot.subtitle = element_text(size = 10, color = "#475569"),
        axis.text.y = element_text(size = 8.5, color = "#334155", family = "mono"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank())

print(plot_event_timeline_biotime)


#STEP 7: Combining databases   ------------------------------------------------------------
##Observation-level ------------------------------------------------------------
#looking into the names of the columns to decide which ones I need to change
colnames(biotime_obs)
colnames(observations_predicts)

#rename BioTIME columns to match PREDICTS naming
biotime_obs_renamed <- biotime_obs %>%
  rename(Latitude = LATITUDE, 
         Longitude = LONGITUDE, 
         Genus = GENUS,
         Family = family, 
         Order = order, 
         Class = class, 
         Phylum = phylum,
         Kingdom = kingdom, 
         Species = valid_name, #this is complicated! BioTIME does not have a Species column, even though it report many observations as "Species" resolution. In value_name we find the name of the species when the resolution is species
         Rank = resolution,
         Sample_start_year = YEAR, 
         Sampling_method = SUMMARY_METHODS,
         Title = Publication_title)

#when Rank is not species, we turn Species to NA
biotime_obs_renamed <- biotime_obs_renamed %>%
  mutate(Species = ifelse(Rank != "Species", NA, Species))

#predicts observations have a _mod value for Latitude and Longitud that comes from after the curation process of script 2. 
#for simplification, we are changing the name of the original ones 
predicts_obs_renamed <- observations_predicts %>%
  rename(Longitud_pre_curation = Longitude, 
         Latitude_pre_curation = Latitude,
         Lat_Long_pre_curation = Lat_Long)

#and now we assign the modified after curation value here
predicts_obs_renamed <- predicts_obs_renamed %>%
  rename(Coordinate_ID = Lat_Long_mod,
         Longitude = Longitude_mod, 
         Latitude = Latitude_mod)

#some fields PREDICTS stores them as character so we need the same for biotime
biotime_obs_renamed <- biotime_obs_renamed %>%
  mutate(Block = as.character(Block),
         Study_number = as.character(Study_number),
         Site_number = as.character(Site_number))

#now I create the database in which I will combine the 2 datasets by first making sure the PREDICTS dataset has the same column types for those 3 columns
observations_predicts_for_combine <- predicts_obs_renamed %>%
  mutate(Block = as.character(Block),
         Study_number = as.character(Study_number),
         Site_number = as.character(Site_number))

#work out which columns exist in one dataset but not the other, after the renaming above
biotime_only_obs  <- setdiff(names(biotime_obs_renamed), names(observations_predicts_for_combine))
predicts_only_obs <- setdiff(names(observations_predicts_for_combine), names(biotime_obs_renamed))

#the unmatched columns I will prefix so nothing gets confused after we combine
biotime_obs_renamed <- biotime_obs_renamed %>%
  rename_with(~ paste0("BioTIME_", .x), all_of(biotime_only_obs)) %>% #here I add the prefix only to the ones that are only in biotime
  mutate(Database = "BioTIME") #I add a column with the name of the database

observations_predicts_for_combine <- observations_predicts_for_combine %>%
  rename_with(~ paste0("PREDICTS_", .x), all_of(predicts_only_obs)) %>% #here I add the prefix only to the ones that are only in predicts
  mutate(Database = "PREDICTS") #I add a column with the name of the database

#stack both datasets 
#shared columns (like Latitude, Measurement, SSBS) line up in the same column
#unmatched columns are filled with NA
observation_level_combined <- bind_rows(biotime_obs_renamed, observations_predicts_for_combine)

#reorder just for readability
priority_cols_obs <- c("Database", "SSBS", "SS", "SSB", "Source_ID", "Study_number",
                       "Block", "Site_number", "Country", "Latitude", "Longitude",
                       "Sample_start_year", "Rank", "Kingdom", "Phylum", "Class", "Order",
                       "Family", "Genus", "Species", 
                       "Sampling_method", "standardized_sampling_method",
                       "Diversity_metric", "Measurement", "Title")

observation_level_combined <- observation_level_combined %>%
  relocate(any_of(priority_cols_obs))

##Site-level ------------------------------------------------------------
#looking into the names of the columns to decide which ones I need to change
colnames(sites_biotime)
colnames(sites_predicts)

#rename BioTIME columns to match PREDICTS naming
sites_biotime_renamed <- sites_biotime %>% 
  rename(Sample_start_year = YEAR, 
         Sampling_method = standardized_sampling_method,
         Title = Publication_title,
         Genus = GENUS,
         Class = class,
         Order = order,
         Family = family,
         Latitude=LATITUDE,
         Longitude =LONGITUDE,
         Species= valid_name)

#if the value in rank has "species" within its text we keep the value in Species, otherwise we turn it to NA
sites_biotime_renamed <- sites_biotime_renamed %>%
  mutate(Species = ifelse(grepl("species", resolution, ignore.case = TRUE), Species, NA))

colnames(sites_biotime_renamed) #checking names

#we take out some of the columns of the sites in biotime to keep it simpler (as in sites from predicts)
sites_biotime_renamed <- sites_biotime_renamed %>%
  select(-Title_dataset, -Study_number, -Diversity_metric, -resolution, -METHODS, ,-N_class, -N_order,  -N_family,  -N_GENUS)

colnames(sites_biotime_renamed) #checking names

#we add some columns to site_predicts based on observations_predicts to make it more similar to sites_biotime
sites_predicts_renamed <- sites_predicts %>%
  left_join(sites_full_info %>% select(SSBS, Rank, Class, Order, Family, Genus, Species, taxa_richness, Source_ID, Title, Site_number, Block, Longitude_mod, Latitude_mod) %>% distinct(), by = "SSBS")
colnames(sites_predicts_renamed)

sites_predicts_renamed <- sites_predicts_renamed %>%
  rename(TA_original = TA,
         TA = TA_corrected_mod,
         Latitude = Latitude_mod,
         Longitude = Longitude_mod)

#we take out the predictors and others that are not relevant for the site level integration
sites_predicts_renamed <- sites_predicts_renamed %>%
  select(-water_pct_1000m_mod, -urban_pct_1000m_mod, -natural_pct_1000m_mod, -pasture_pct_1000m_mod, ,-crop_pct_1000m, -pesticide_1000m,  
         -fertilizer_1000m,  -crop_diversity_shannon_1000m,  -crop_diversity_machefer_alpha_1000m, -pct_very_small_and_small_field_1000m,
         -pct_very_large_and_large_field_1000m, -climate_mean_temp_c_1000m, -climate_annual_precip_mm_1000m, -Rank)

#fixing type mismatches
#Site_number and Block are ID-like labels so we keep as character on both sides
sites_biotime_renamed <- sites_biotime_renamed %>%
  mutate(Site_number = as.character(Site_number),
         Block = as.character(Block))

sites_predicts_renamed <- sites_predicts_renamed %>%
  mutate(Site_number = as.character(Site_number),
         Block = as.character(Block))

#Latitude/Longitude need to stay numeric for spatial analysis so we convert BioTIME's back to numeric
sites_biotime_renamed <- sites_biotime_renamed %>%
  mutate(Latitude = as.numeric(Latitude),
         Longitude = as.numeric(Longitude))

#rename BioTIME columns to match PREDICTS naming
biotime_only_site  <- setdiff(names(sites_biotime_renamed), names(sites_predicts_renamed))
predicts_only_site <- setdiff(names(sites_predicts_renamed), names(sites_biotime_renamed))

#adding the prefixes and the database colum
sites_biotime_renamed <- sites_biotime_renamed %>%
  rename_with(~ paste0("BioTIME_", .x), all_of(biotime_only_site)) %>%
  mutate(Database = "BioTIME")

sites_predicts_renamed <- sites_predicts_renamed %>%
  rename_with(~ paste0("PREDICTS_", .x), all_of(predicts_only_site)) %>%
  mutate(Database = "PREDICTS")

site_level_combined <- bind_rows(sites_biotime_renamed, sites_predicts_renamed)

priority_cols_site <- c("Database", "SSBS",  "TA", "taxa_richness", "Diversity_metric","SS", "SSB", "Source_ID", "Study_number",
                       "Block", "Site_number", "Country", "Latitude", "Longitude", "Coordinate_ID",
                       "Sample_start_year", "Rank", "Kingdom", "Phylum", "Class", "Order",
                       "Family", "Genus", "Species", 
                       "Sampling_method", "Title")

site_level_combined <- site_level_combined %>%
  relocate(any_of(priority_cols_site))

n_distinct(site_level_combined$SSBS)
n_distinct(site_level_combined$Coordinate_ID)

##Save ------------------------------------------------------------
saveRDS(observation_level_combined, "Intermediate_dataset/OL_EU_UK_Ins_Arc.rds")
saveRDS(site_level_combined,"Intermediate_dataset/SL_EU_UK_Ins_Arc.rds")
