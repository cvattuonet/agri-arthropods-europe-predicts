#---------------------------------------------------------------------------------------------------------------------------------------#
# CODE DETAILS                                                                                                                          #
#                                                                                                                                       #
# Author: Catalina Vattuone  (cvattuonet@gmail.com)                                                                                     #        
# Date of latest update: 19-08-2026                                                                                                     #
# Main use: Study on "Assessing how agricultural practices and landscape composition shape differences in arthropod biodiversity among  #  
#           habitat types in Europe using the PREDICTS database". Internship at INRAE                                                   #
#                                                                                                                                       #
# Content: Estimation of predictors based on the extraction of cell values from multiple agronomic databases                            #
#          within buffers of 10.000 mts around PREDICTS insects and arachnids sampling points                                           #                           
#---------------------------------------------------------------------------------------------------------------------------------------#

#LIBRARIES  ------------------------------------------------------------
pacman::p_load(terra,readr,dplyr, tibble,parallel,ecmwfr)

#LOAD ------------------------------------------------------------
sites_df <- read_delim("Intermediate_dataset/sites_mod.csv")

#FIXED EFFECTS CALCULATION ------------------------------------------------------------
#For each fixed effect the process is similar
#Step 1: We load the database
#Step 2: If necesary, we reclasify the raster
#Step 3: We define the radii in which values will be extracted and computed (we use 1000 but can be changed to other or to various)
#Step 4: We create a loop in which, by radii (if someone wants later to add more radiis) and chunk of sites, we create a buffer, and extract the values of cells directly as mean in some cases, as direct values in others to compute internally the indicator we are looking for. 
#        In the case of cropland cover and other covers from Corine Land Cover the loop has an extra step of years because the resolution is fine and takes longer
#Step 5: We save the intermediate results and then put them together
#Step 6: We add all the fixed effects to the biodiversity dataset obtained previously. 

##CROPLAND COVER ------------------------------------------------------------
#load 
#directory with all the tiles. Files are grouped within many subfolders
cropland_root_dir <- "Databases/Cropland_Liao/"
all_tile_paths <- list.files(cropland_root_dir, pattern = "\\.tif$", recursive = TRUE, full.names = TRUE)

#projection matching
sites_spatial        <- vect(sites_df, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")
first_tile           <- rast(all_tile_paths[1])
sites_projected_liao <- project(sites_spatial, crs(first_tile))

#master bounding box setup for tile filtering
europe_extent_buffered <- ext(ext(sites_projected_liao)) + 0.1
eur_xmin <- xmin(europe_extent_buffered); eur_xmax <- xmax(europe_extent_buffered)
eur_ymin <- ymin(europe_extent_buffered); eur_ymax <- ymax(europe_extent_buffered)

#tile filtering to only keep those that intersect with Europe
is_in_europe <- sapply(all_tile_paths, function(path) {
  t_ext <- ext(rast(path))
  return(t_ext$xmin <= eur_xmax && t_ext$xmax >= eur_xmin &&
           t_ext$ymin <= eur_ymax && t_ext$ymax >= eur_ymin)
})

relevant_tile_paths <- all_tile_paths[is_in_europe]

#precomputing each tile's real geographic extent for correct tile-name tracking
tile_extents <- tibble(tile_path = relevant_tile_paths,
  xmin = sapply(relevant_tile_paths, function(p) ext(rast(p))$xmin),
  xmax = sapply(relevant_tile_paths, function(p) ext(rast(p))$xmax),
  ymin = sapply(relevant_tile_paths, function(p) ext(rast(p))$ymin),
  ymax = sapply(relevant_tile_paths, function(p) ext(rast(p))$ymax))

#parameters
num_cores       <- detectCores() - 1. #detecting the number of cores available for parallel process given the heavy of the process
chunk_size_crop <- 50 #creating chunks of this size to avoid memory issues when extracting from the tiles
radii           <- c(1000) #if other radii want to be tried, if should be changed here

#creating a clean tibble to store the results
cropland_results <- tibble(.rows = nrow(sites_df))

#loop creating the columns to store the results. Important if there are more than one radii
for (r in radii) {
  col_name <- paste0("crop_pct_", r, "m")
  cropland_results[[col_name]] <- NA_real_
}

cropland_results$assigned_tile_name <- NA_character_ 

#creating a directory to save the intermediate results
dir.create("Intermediate_dataset/fixed_effects_intermediate/", showWarnings = FALSE)

#selecting all the years to run the extraction loop. Loop will go year by year 
unique_years <- sort(unique(sites_df$Sample_start_year[sites_df$Sample_start_year >= 2000 & !is.na(sites_df$Sample_start_year)]))

#Year by year, and chunk by chunk, it selects sites, creates a a 1000 mts buffer, looks for the correct band in the Liao layer to extract the values
for (yr in unique_years) {
  cat("PROCESSING SITES FOR SAMPLING YEAR:", yr, "\n") #year being process
  year_indices <- which(sites_df$Sample_start_year == yr)  #the index of the year
  num_sites    <- length(year_indices) #estimates the amount of sites in that year
  cat("Total sites for this year:", num_sites) #reports the amoung of sites of the year being process
  
  if (num_sites == 0) next #moves to the next year rapidly if there is not sites in that year
  
  target_band <- yr - 1999 #The Liao raster has bands for each year so we have to define what is the right band to look at to extract the values 
  current_points <- sites_projected_liao[year_indices, ]   #selects the sites of the year being processed and projects them to the Liao raster projection
  
  tile_collection <- vrt(relevant_tile_paths, options = c("-b", as.character(target_band)))  #selecting the right band for the year being processed
  chunks <- split(1:num_sites, ceiling(seq_along(1:num_sites) / chunk_size_crop)) #defining the chunks to avoid memory issues when extracting from the tiles
  
  for (r in radii) { #in case many radiis are used
    use_exact <- (r <= 1000) #only for radii 1000 and less than 10000 it keeps the exact = TRUE category in the extraction package. For larger it will extract only the points with the center inside the buffer
    year_polygons <- buffer(current_points, width = r) #creates the buffer for the sites of the year being processed

    chunk_results_list <- mclapply(seq_along(chunks), function(ch_idx) {   #process each chunk in parallel
      sub_indices    <- chunks[[ch_idx]]
      chunk_polygons <- year_polygons[sub_indices, ]
      
      if (use_exact) {
        extracted <- terra::extract(tile_collection, chunk_polygons, exact = TRUE) #extract the values from the correct bank, for each site in the chunk, with the exact = TRUE option to consider the weights of the pixels that are partially inside the buffer 
        if (nrow(extracted) == 0) return(rep(0, length(sub_indices))) #if there is no values extracted, it returns 0 for all the sites in the chunk
        
        colnames(extracted) <- c("ID", "value", "weight") #defines the column names of the extracted values
        id_factor           <- factor(extracted$ID, levels = 1:length(sub_indices)) #creates a factor to group the extracted values by site ID
        
        binary_vals   <- ifelse(extracted$value == 10, 1, 0) #creates a binary vector to identify the pixels that are cropland (value 10 in Liao)
        weighted_sums <- tapply(binary_vals * extracted$weight, id_factor, sum, na.rm = TRUE) #calculates the weighted sum of cropland pixels for each site in the chunk
        weight_sums   <- tapply(extracted$weight, id_factor, sum, na.rm = TRUE) #calculates the total weight of pixels for each site in the chunk
        
        pcts <- (weighted_sums / weight_sums) * 100 #calculates the percentage of cropland for each site in the chunk
        pcts[is.na(pcts) | is.nan(pcts)] <- 0 #replaces NA and NaN values with 0
        return(as.numeric(pcts)) #returns the percentage of cropland for each site in the chunk
        
      } else { #this part repeats the proces but for radii that over 1000mts where the exact = FALSE option is used in the extraction. The weights are not considered and only the pixels that have their center inside the buffer are used
        extracted <- terra::extract(tile_collection, chunk_polygons, exact = FALSE)
        if (nrow(extracted) == 0) return(rep(0, length(sub_indices)))
        
        colnames(extracted) <- c("ID", "value")
        id_factor           <- factor(extracted$ID, levels = 1:length(sub_indices))
        
        crop_counts  <- tapply(extracted$value == 10, id_factor, sum, na.rm = TRUE)
        valid_counts <- tapply(!is.na(extracted$value), id_factor, sum, na.rm = TRUE)
        
        pcts <- (crop_counts / valid_counts) * 100
        pcts[is.na(pcts) | is.nan(pcts)] <- 0
        return(as.numeric(pcts))
      }
    }, mc.cores = num_cores) #process each chunk in parallel
    
    col_name <- paste0("crop_pct_", r, "m") #save results into the external tibble
    cropland_results[year_indices, col_name] <- unlist(chunk_results_list) 
  }
  
  #here we define the coordinates of the sites to check which tile they belong to later
  site_coords <- crds(current_points)
  
  #in this function we check for each site which tile contains its coordinates and return the tile name. If no tile contains the site, we return NA.
  assigned_tile <- sapply(seq_len(nrow(site_coords)), function(i) {
    x <- site_coords[i, 1]; y <- site_coords[i, 2]
    match_idx <- which(x >= tile_extents$xmin & x <= tile_extents$xmax &
                         y >= tile_extents$ymin & y <= tile_extents$ymax)
    if (length(match_idx) == 0) return(NA_character_)
    basename(tile_extents$tile_path[match_idx[1]])
  })
  
  #we save the tile name for each site in the results tibble
  cropland_results$assigned_tile_name[year_indices] <- assigned_tile
  
  #intermediate checkpoint save
  temp_checkpoint <- bind_cols(sites_df, cropland_results)
  saveRDS(temp_checkpoint, paste0("Intermediate_dataset/fixed_effects_intermediate/sites_df_cropland_checkpoint_", yr, ".rds")) 
  rm(tile_collection, year_polygons, temp_checkpoint) #remove large objects to free memory
  gc() #clean memory 
}

cropland_results <- cropland_results %>% #we add teh SSBS name to the results so that later we can join them all having the value.
  mutate(SSBS = sites_df$SSBS)  

#save 
saveRDS(cropland_results, "Intermediate_dataset/fixed_effects_intermediate/fixed_effect_cropland.rds")

##WATER COVER ------------------------------------------------------------
#load
waw_100m <- rast("Databases/Water_Copernicus/WAW_2018_100m_eu_03035_V1_0.tif") #we used the 100mts resolution layer for speed but there is also a layer with more resolution available

#reclassification of the water mask to keep only pemanent and temporary water and wetness
water_mask <- classify(waw_100m, rcl = matrix(c(
  0,   0,   #dry -> 0
  1,   1,   #permanent water -> 1 
  2,   1,   #temporary water -> 1
  3,   1,   #permanent wet -> 1
  4,   1,   #temporary wet -> 1
  253, NA,  #sea water -> NA
  254, NA,  #unclassifiable -> NA
  255, NA   #outside area -> NA
), ncol = 2, byrow = TRUE))

#Projection matching 
waw_crs             <- crs(water_mask)
sites_spatial       <- vect(sites_df, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")
sites_projected_waw <- project(sites_spatial, waw_crs)

#parameters setup
total_sites_wow <- nrow(sites_projected_waw)  #amount of sites
chunk_size_wow  <- 500 #size of the chunks
radii           <- c(1000)

#creating a clean tibble to store results
wow_results <- tibble(SSBS = sites_df$SSBS,.rows = total_sites_wow)

#run the loop following the same logic as case 1 except we do not go by year because the resolution is coarser
for (r in radii) {
  site_polygons <- buffer(sites_projected_waw, width = r)
  col_pct  <- paste0("water_pct_", r, "m")
  col_area <- paste0("water_area_m2_", r, "m")
  
  if (r <= 1000) { 
    pct_results <- numeric(total_sites_wow)
    chunks      <- split(1:total_sites_wow, ceiling(seq_along(1:total_sites_wow) / chunk_size_wow))
    
    for (i in seq_along(chunks)) {
      indices <- chunks[[i]]
      cat(sprintf("  -> Processing chunk %d/%d (Sites %d to %d)... ",  i, length(chunks), indices[1], indices[length(indices)]))
      
      extracted_chunk <- terra::extract(water_mask, site_polygons[indices, ], fun = mean, na.rm = TRUE,exact = TRUE) #we use mean because given the layer is 1 and 0 s, mean will give us the percentage of water in the buffer
      pct_results[indices] <- extracted_chunk[, 2] #no need to use the weights for the calculation as the layer is already 1 and 0 s, so mean will give us the percentage of water in the buffer
    }
    
    wow_results[[col_pct]] <- pct_results * 100 #we only need to multiply by 100 to obtain the percentage
    
  } else {
    extracted <- terra::extract(water_mask, site_polygons, fun = mean, na.rm = TRUE, exact = FALSE) #if the radii is more than 1000 we only use the cell that have their center inside the buffer.
    wow_results[[col_pct]] <- extracted[, 2] * 100
  }
  
  #we also want to save the area of water (not only the percentage so total area of the buffer is estimated and with the percentage the area of water
  total_buffer_area_m2    <- pi * (r^2)
  wow_results[[col_area]] <- total_buffer_area_m2 * (wow_results[[col_pct]] / 100)
}

#save results
saveRDS(wow_results, "Intermediate_dataset/fixed_effects_intermediate/fixed_effect_wow.rds")

##SEA WATER COVER  ------------------------------------------------------------
#This is not a predictor that will be used, but an indicator we need to rescale later the land covers of the landscape.

#we load the same later
waw_100m <- rast("Databases/Water_Copernicus/WAW_2018_100m_eu_03035_V1_0.tif")

#reclassification to have only the sea water
sea_mask <- classify(waw_100m, rcl = matrix(c(
  0,   0,   #dry -> 0
  1,   0,   #permanent water -> 0 (already counted in water_pct)
  2,   0,   #temporary water -> 0
  3,   0,   #permanent wet -> 0
  4,   0,   #temporary wet -> 0
  253, 1,   #sea water -> 1
  254, NA,  #unclassifiable -> NA
  255, NA   #outside area -> NA
), ncol = 2, byrow = TRUE))


#Projection matching 
waw_crs             <- crs(sea_mask)
sites_spatial       <- vect(sites_df, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")
sites_projected_waw <- project(sites_spatial, waw_crs)

#parameters setup
total_sites_wow <- nrow(sites_projected_waw)
chunk_size_wow  <- 500
radii           <- c(1000)

#creating the clean table to store results
sea_results <- tibble(SSBS = sites_df$SSBS, .rows = total_sites_wow)

for (r in radii) { #same loop as for water
  cat("STARTING SEA WATER EXTRACTION FOR RADIUS:", r, "METERS\n")
  site_polygons <- buffer(sites_projected_waw, width = r)
  col_pct <- paste0("sea_pct_", r, "m")
  
  if (r <= 1000) {
    pct_results <- numeric(total_sites_wow)
    chunks <- split(1:total_sites_wow, ceiling(seq_along(1:total_sites_wow) / chunk_size_wow))
    
    for (i in seq_along(chunks)) {
      indices <- chunks[[i]]
      cat(sprintf("  -> Processing chunk %d/%d (Sites %d to %d)... ",  i, length(chunks), indices[1], indices[length(indices)]))
      extracted_chunk <- terra::extract(sea_mask, site_polygons[indices, ],  fun = mean, na.rm = TRUE,exact = TRUE) #layer is 1 and 0s so mean will give us directly the proportion
      pct_results[indices] <- extracted_chunk[, 2]
    }
    
    sea_results[[col_pct]] <- pct_results * 100 #we transform the proportion to percentage
    
  } else {
    extracted <- terra::extract(sea_mask, site_polygons, fun = mean, na.rm = TRUE, exact = FALSE) #again false if radii is over 1000 to save speed
    sea_results[[col_pct]] <- extracted[, 2] * 100
  }
}

#save results
saveRDS(sea_results, "Intermediate_dataset/fixed_effects_intermediate/fixed_effect_sea.rds")

##CROP DIVERSITY SHANNON  ------------------------------------------------------------

#load 
shannon_europe <- crop(rast("Intermediate_dataset/Crop_Diversity_CROPGRIDS/shannon_map.tif"),  
                       ext(vect(sites_df, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")) + 0.15)

#projection matching
sites_spatial_shannon <- vect(sites_df, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")
sites_proj_shannon    <- project(sites_spatial_shannon, crs(shannon_europe))

#parameters
radii <- c(1000)

#creating clean tibble to store results
shannon_results <- tibble(SSBS = sites_df$SSBS, .rows = nrow(sites_df))

#run for every radii. We don't separate by years as in case 1 of croplands because the resolution is coarser (faster)
for (r in radii) {
  cat("STARTING SHANNON INDEX EXTRACTION FOR RADIUS:", r, "METERS\n")
  poly <- buffer(sites_proj_shannon, width = r)
  col_name <- paste0("crop_diversity_shannon_", r, "m") #create a dynamic destination column name
  
  #adapt extraction performance strategy based on scale (matches water loop logic)
  if (r <= 1000) { 
    ext_data <- terra::extract(shannon_europe, poly, fun = mean, na.rm = TRUE, exact = TRUE)  #fun = mean to obtain the average value in the buffer
  } else {
    ext_data <- terra::extract(shannon_europe, poly, fun = mean, na.rm = TRUE, exact = FALSE)
  }
  
  #save the results directly into the tibble
  shannon_results[[col_name]] <- ext_data[, 2] #here we dont need to multiply by 100 because we keet the average directly 
}

#save results
saveRDS(shannon_results, "Intermediate_dataset/fixed_effects_intermediate/fixed_effect_shannon.rds")

##URBAN, NATURAL VEGETATION AND PASTURE COVER  ------------------------------------------------------------

#directory and list of files
clc_root_dir <- "Databases/CLC_all_years/"
sample_tile    <- list.files(clc_root_dir, pattern = "\\.tif$", full.names = TRUE)[1]

#Corine Cover codes per category
urban_codes       <- 1:11     
veg_codes         <- 23:29    
pasture_codes     <- 18       

#projection matching
clc_crs        <- crs(rast(sample_tile))
sites_spatial  <- vect(sites_df, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")
sites_metric   <- project(sites_spatial, clc_crs)

#parameters
radii        <- c(1000)
chunk_size   <- 50  #only 50 because the layer takes longer
num_cores    <- detectCores() - 1 #this predictor takes long so we use parallel processing

#creating clean tibbles for the results of each category
urban_results       <- tibble(.rows = nrow(sites_df))
veg_results         <- tibble(.rows = nrow(sites_df))
pasture_results     <- tibble(.rows = nrow(sites_df))

#we create the columns to store the results of each category
for (r in radii) { 
  urban_results[[paste0("urban_pct_", r, "m")]]         <- NA_real_
  veg_results[[paste0("natural_pct_", r, "m")]]         <- NA_real_
  pasture_results[[paste0("pasture_pct_", r, "m")]]     <- NA_real_
}

#metadata tracking rows
urban_results$assigned_clc_file         <- NA_character_
veg_results$assigned_clc_file           <- NA_character_
pasture_results$assigned_clc_file       <- NA_character_

#given it takes more time, we follow the cropland procedure of adding a loop per year
unique_years <- sort(unique(sites_df$Sample_start_year[sites_df$Sample_start_year >= 2000 & !is.na(sites_df$Sample_start_year)]))

for (yr in unique_years) {
  cat("PROCESSING CORINE LAND COVER FOR SAMPLING YEAR:", yr, "\n")
  year_indices <- which(sites_df$Sample_start_year == yr)
  num_sites    <- length(year_indices)
  cat("Total sites for this year:", num_sites, "\n")
  if (num_sites == 0) next
 
  #depending on the year that is being analysed, we select the correct file. We have Corine Land Cover files of yeats 2000, 2006, 2012, an 2018. Given our data we only use the first 2 but in other situation the others could be needded
  target_file <- case_when(
    yr >= 2000 & yr <= 2005 ~ "U2006_CLC2000_V2020_20u1.tif",
    yr >= 2006 & yr <= 2011 ~ "U2012_CLC2006_V2020_20u1.tif",
    yr >= 2012 & yr <= 2017 ~ "U2018_CLC2012_V2020_20u1.tif",
    yr >= 2018              ~ "U2018_CLC2018_V2020_20u1.tif",
    TRUE                    ~ NA_character_)
  
  if (is.na(target_file)) { #if for some reason we have a database that has other year it is skip here.
    warning(paste("Year", yr, "falls outside specified CLC assignment ranges! Skipping."))
    next
  }
  
  current_raster <- rast(paste0(clc_root_dir, target_file))  #here we load the correct file given what we decided on target_file
  current_points <- sites_metric[year_indices, ] 
  
  chunks <- split(1:num_sites, ceiling(seq_along(1:num_sites) / chunk_size)) 
  
  for (r in radii) { #here we follow the same procedure as always
    
    use_exact <- (r <= 1000)
    year_polygons <- buffer(current_points, width = r)
    
    chunk_results_list <- mclapply(seq_along(chunks), function(ch_idx) {
      sub_indices    <- chunks[[ch_idx]]
      chunk_polygons <- year_polygons[sub_indices, ]
      
      extracted <- terra::extract(current_raster, chunk_polygons, exact = use_exact) #the extraction step
      
      if (nrow(extracted) == 0) {
        return(tibble(urban = rep(0, length(sub_indices)),  natural = rep(0, length(sub_indices)),  pasture = rep(0, length(sub_indices))))
      }
      
      if (use_exact) { #different situations 
        colnames(extracted) <- c("ID", "value", "weight")
      } else {
        colnames(extracted) <- c("ID", "value")
        extracted$weight    <- 1 
      }
      
      id_factor <- factor(extracted$ID, levels = 1:length(sub_indices))
      
      #we have have the weather of each pixel in the buffer, and we can calculate the weighted sum of each category and the total weight to calculate the percentage of each category in the buffer
      urb_weights     <- tapply(ifelse(extracted$value %in% urban_codes, extracted$weight, 0), id_factor, sum, na.rm = TRUE)
      nat_weights     <- tapply(ifelse(extracted$value %in% veg_codes, extracted$weight, 0), id_factor, sum, na.rm = TRUE)
      pas_weights     <- tapply(ifelse(extracted$value %in% pasture_codes, extracted$weight, 0), id_factor, sum, na.rm = TRUE)
      total_weights   <- tapply(extracted$weight, id_factor, sum, na.rm = TRUE)
      
      #percentage estimation
      urb_pct <- (urb_weights / total_weights) * 100
      nat_pct <- (nat_weights / total_weights) * 100
      pas_pct <- (pas_weights / total_weights) * 100
      
      #we give them 0 if there is no value
      urb_pct[is.na(urb_pct)] <- 0; nat_pct[is.na(nat_pct)] <- 0; pas_pct[is.na(pas_pct)] <- 0
      
      #we get the results in a tibble to be able to bind them later
      return(tibble(urban = as.numeric(urb_pct), natural = as.numeric(nat_pct), pasture = as.numeric(pas_pct)))
    }, mc.cores = num_cores)
    
    #combine the results of all the chunks into a single tibble
    combined_chunks <- bind_rows(chunk_results_list)
    
    #adding the results to the main results tibbles
    urban_results[year_indices, paste0("urban_pct_", r, "m")] <- combined_chunks$urban
    veg_results[year_indices, paste0("natural_pct_", r, "m")] <- combined_chunks$natural
    pasture_results[year_indices, paste0("pasture_pct_", r, "m")] <- combined_chunks$pasture
  }
  
  #we save the name of the file used for each site in the results tibbles
  urban_results$assigned_clc_file[year_indices] <- target_file
  veg_results$assigned_clc_file[year_indices]  <- target_file
  pasture_results$assigned_clc_file[year_indices] <- target_file
  
  #clean up memory
  rm(current_raster, year_polygons) 
  gc()
}

#adding SSBS for later joining all predictors
urban_results         <- urban_results %>% mutate(SSBS = sites_df$SSBS)
veg_results           <- veg_results %>% mutate(SSBS = sites_df$SSBS)
pasture_results       <- pasture_results %>% mutate(SSBS = sites_df$SSBS)

#save outputs
saveRDS(urban_results,        "Intermediate_dataset/fixed_effects_intermediate/fixed_effect_urban.rds")
saveRDS(veg_results,          "Intermediate_dataset/fixed_effects_intermediate/fixed_effect_natural_veg.rds")
saveRDS(pasture_results,      "Intermediate_dataset/fixed_effects_intermediate/fixed_effect_pastures.rds")


##ALPHA CROP DIVERSITY------------------------------------------------------------
#IMPORTANT NOTE: This was the code run for the predictor, however, I considers all pixels including zeros (where no cropland pixels). This creates a problem as the predict is signaling this way not only diversity but cropland presence
#The analysis that is presented later is trying to capture the correlation between the results obtained with this method and the results obtained if zeros were excluded to check how much the compounding effect could be playing a role

#load
machefer_alpha <- rast("Databases/Crop_Diversity_Machefer/crop_div_alpha_1km.tif")

#projection matching
sites_spatial_jrc <- vect(sites_df, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")
sites_proj_jrc    <- project(sites_spatial_jrc, crs(machefer_alpha))

#parameters
radii <- c(1000)

#creating a clean tibble
jrc_results <- tibble(SSBS = sites_df$SSBS, .rows = nrow(sites_df))

#run the loop as previously
for (r in radii) {
  poly <- buffer(sites_proj_jrc, width = r)   #generate polygon boundaries for this distance
  col_name <- paste0("crop_diversity_machefer_alpha_", r, "m")#create a destination column
  if (r <= 1000) {  #adapt extraction performance strategy based on scale (if radii is less than 1000, we weight the pixels that are partially inside the buffer, if not we only consider the pixels that have their center inside the buffer)
    ext_data <- terra::extract(machefer_alpha, poly, fun = mean, na.rm = TRUE, exact = TRUE)
  } else {
    ext_data <- terra::extract(machefer_alpha, poly, fun = mean, na.rm = TRUE, exact = FALSE)
  }
  jrc_results[[col_name]] <- ext_data[, 2]#save the results directly into our clean results tibble
}

#save results
saveRDS(jrc_results, "Intermediate_dataset/fixed_effects_intermediate/fixed_effect_jrc.rds")

#analysis when only cropland is present, see IMPORTANT NOTE above
sites_proj_jrc <- project(vect(sites_df, geom = c("Longitude","Latitude"), crs = "EPSG:4326"),crs(machefer_alpha))
poly <- buffer(sites_proj_jrc, width = 1000)
ex <- terra::extract(machefer_alpha, poly, exact = TRUE)
names(ex) <- c("ID", "alpha", "w")

jrc_new <- ex %>%
  filter(!is.na(alpha)) %>%
  group_by(ID) %>%
  summarise(crop_div_cropcells = if (any(alpha > 0))  #diversity where cropland is present (zeros excluded)
      sum(alpha[alpha > 0] * w[alpha > 0]) / sum(w[alpha > 0]) else NA_real_,
    crop_cell_share    = sum(w[alpha > 0]) / sum(w), #share of buffer with any cropland 
    crop_div_buffermean = sum(alpha * w) / sum(w),  #original composite to compare
    .groups = "drop") %>%
  mutate(SSBS = sites_df$SSBS[ID])

correlation1 <- cor(jrc_new$crop_div_buffermean, jrc_new$crop_div_cropcells, use = "complete.obs")
correlation2 <- cor(jrc_new$crop_div_buffermean, jrc_new$crop_cell_share, use = "complete.obs")
print(correlation1)
print(correlation2)

##PESTICIDES   ------------------------------------------------------------

#load
total_pesticide_M <- rast("Intermediate_dataset/Pesticides/total_pesticide_medium_scenario.tif")
total_pesticide_M[is.na(total_pesticide_M)] <- 0 #where we have NA we assign 0 because it means that there is no pesticide application in that pixel

#projection matching
sites_spatial   <- vect(sites_df, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")
raster_crs      <- crs(total_pesticide_M) 
sites_projected <- project(sites_spatial, raster_crs)

#parameters
total_sites <- nrow(sites_projected)
chunk_size  <- 500
radii       <- c(1000)

#creating a clean tibble for storing the results
pesticide_results <- tibble(SSBS = sites_df$SSBS, .rows = nrow(sites_df))

#run
for (r in radii) {
  site_polygons <- buffer(sites_projected, width = r) #generate circular polygon 
  col_name <- paste0("pesticide_", r, "m") #create a column for the results

  if (r <= 1000) {   #adapt extraction on scale
    radius_results <- numeric(total_sites)
    chunks <- split(1:total_sites, ceiling(seq_along(1:total_sites) / chunk_size))
    
    for (i in seq_along(chunks)) {
      indices <- chunks[[i]]
      extracted_chunk <- terra::extract(total_pesticide_M,site_polygons[indices, ],fun = mean,na.rm = TRUE,exact = TRUE)
      radius_results[indices] <- extracted_chunk[, 2]
    }
    pesticide_results[[col_name]] <- radius_results #save to tibble
    
  } else {
    extracted <- terra::extract(total_pesticide_M,site_polygons,fun = mean,na.rm = TRUE, exact = FALSE)
    pesticide_results[[col_name]] <- extracted[, 2]  #save to tibble
  }
}

#save results
saveRDS(pesticide_results, "Intermediate_dataset/fixed_effects_intermediate/fixed_effect_pesticides.rds")

##NITROGEN FERTILIZATION ------------------------------------------------------------
#load
npk_grids <- rast("Intermediate_dataset/Fertilizers_NPKGRIDS/n_areal_load_kg_per_ha_land.tif")
npk_grids[is.na(npk_grids)] <- 0

#projection matching
sites_spatial   <- vect(sites_df, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")
raster_crs      <- crs(npk_grids) 
sites_projected <- project(sites_spatial, raster_crs)

#parameters
total_sites <- nrow(sites_projected)
chunk_size  <- 500
radii       <- c(1000)

#creating clean tibble for storing results
fertilization_results <- tibble(SSBS = sites_df$SSBS, .rows = nrow(sites_df))

#run
for (r in radii) {
  site_polygons <- buffer(sites_projected, width = r) #generate buffer
  col_name <- paste0("fertilizer_", r, "m")   #create column to put results
  if (r <= 1000) {   #adapt extraction based on scale as in previous cases
    radius_results <- numeric(total_sites)
    chunks <- split(1:total_sites, ceiling(seq_along(1:total_sites) / chunk_size))
    
    for (i in seq_along(chunks)) {
      indices <- chunks[[i]]
      extracted_chunk <- terra::extract(npk_grids,site_polygons[indices, ],fun = mean,na.rm = TRUE,exact = TRUE)
      radius_results[indices] <- extracted_chunk[, 2]
    }
    
    fertilization_results[[col_name]] <- radius_results #save chunked local results to tibble
    
  } else {
    extracted <- terra::extract(npk_grids,site_polygons,fun = mean,na.rm = TRUE,exact = FALSE)
    fertilization_results[[col_name]] <- extracted[, 2] #save chunked local results to our independent tibble
  }
}

#save results
saveRDS(fertilization_results, "Intermediate_dataset/fixed_effects_intermediate/fixed_effect_field_fertilization.rds")

##FIELD SIZE SMALL AND LARGE FIELDS PREVALENCE  ------------------------------------------------------------
#load
field_raster <- rast("Databases/FieldSize_Lesiv/dominant_field_size_categories.tif")

#projection matching
field_crs             <- crs(field_raster)
sites_spatial         <- vect(sites_df, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")
sites_projected_field <- project(sites_spatial, field_crs)

#parameters
radii       <- c(1000)
total_sites <- nrow(sites_projected_field)

#creating clean tibble for storing results to store results
field_size_results <- tibble(SSBS = sites_df$SSBS, .rows = nrow(sites_df))

#run extraction  with an additional step to calculate the metric we are interested on
for (r in radii) {
  site_polygons <- buffer(sites_projected_field, width = r)  #generate circular buffers
  
  #extraction of values 
  if (r <= 1000) {   #extract pixel configurations based on scale
    extracted_raw <- terra::extract(field_raster, site_polygons, ID = TRUE, exact = TRUE)
    colnames(extracted_raw) <- c("ID", "clc_value", "weight")
  } else {
    extracted_raw <- terra::extract(field_raster, site_polygons, ID = TRUE, exact = FALSE)
    colnames(extracted_raw) <- c("ID", "clc_value")
    extracted_raw$weight <- 1 
  }
  
  #we estimate the pct for each category in the raster
  pct_summary <- extracted_raw %>%
    filter(!is.na(clc_value)) %>% 
    group_by(ID) %>%
    summarise(total_weight = sum(weight),
      pct_no_field         = if_else(total_weight > 0, (sum(weight[clc_value == 0]) / total_weight) * 100, 0),
      pct_very_large_field = if_else(total_weight > 0, (sum(weight[clc_value == 3502]) / total_weight) * 100, 0),
      pct_large_field      = if_else(total_weight > 0, (sum(weight[clc_value == 3503]) / total_weight) * 100, 0),
      pct_medium_field     = if_else(total_weight > 0, (sum(weight[clc_value == 3504]) / total_weight) * 100, 0), #this is not acyually needed but I leave it as someone might be interested in extracting it
      pct_small_field      = if_else(total_weight > 0, (sum(weight[clc_value == 3505]) / total_weight) * 100, 0),
      pct_very_small_field = if_else(total_weight > 0, (sum(weight[clc_value == 3506]) / total_weight) * 100, 0),
      .groups = "drop" ) %>%
    mutate( pct_very_small_and_small_field = pct_very_small_field + pct_small_field, #then we estimate the indicators we are going to use
      pct_very_large_and_large_field = pct_very_large_field + pct_large_field)

  
  #metrics we actually want to keep
  target_metrics <- c("pct_very_small_and_small_field", "pct_very_large_and_large_field")
  
  for (metric in target_metrics) {
    col_name <- paste0(metric, "_", r, "m")
    vals <- rep(0, total_sites)  #default to 0 for sites with no valid pixels in the buffer
    vals[pct_summary$ID] <- pct_summary[[metric]]
    field_size_results[[col_name]] <- vals
  }
}

#save results
saveRDS(field_size_results, "Intermediate_dataset/fixed_effects_intermediate/fixed_effect_field_size.rds")


##CLIMATE ------------------------------------------------------------
#Data was first downloaded from the website and then loaded for use. 
#Here, I leave the code for downloading the data but only commented. User and keys are left as "XXX" for anyone to complete with their own credentials.

# wf_set_key( user = "XXX", key  = "XXX")
# target_years <- as.character(2000:2016)
# 
# cds_request <- list(
#   dataset_short_name = "reanalysis-era5-land-monthly-means",
#   product_type  = "monthly_averaged_reanalysis",
#   variable = c("2m_temperature", "total_precipitation"),
#   year  = target_years,
#   month  = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"),
#   time = "00:00",
#   area  = c(71.5, -12.5, 34.0, 35.0),  
#   format = "netcdf",
#   target= "era5_land_europe_monthly_2000_2016.nc")
# 
# #directory creation
# if(!dir.exists("Databases/Climate_Copernicus")) {
#   dir.create("Databases/Climate_Copernicus")
# }
# 
# #request
# wf_request(request  = cds_request, user = "XXX",transfer = TRUE,path = "Databases/Climate_Copernicus")
# 
# #load 
# zip_file <- "Databases/Climate_Copernicus/era5_land_europe_monthly_2000_2016.zip"
# unzip(zip_file, exdir = "Databases/Climate_Copernicus")
# extracted_files <- list.files("Databases/Climate_Copernicus", pattern = "\\.nc$", full.names = TRUE)
# if(length(extracted_files) > 0) {
#   climate_raw <- rast(extracted_files[1])
#   print(climate_raw)
#   message("Success! The climate raster stack is loaded.")
# } else {
#   stop("Extraction failed or no .nc file was found inside the zip archive.")
# }


#load
climate_raw <- rast("Databases/Climate_Copernicus/data_stream-moda.nc")

#subset from the loaded data
temp_indices   <- grep("^t2m", names(climate_raw))
precip_indices <- grep("^tp", names(climate_raw))
temp_kelvin <- subset(climate_raw, temp_indices)
precip_m    <- subset(climate_raw, precip_indices)

#units conversion
temp_celsius <- temp_kelvin - 273.15
precip_mm    <- precip_m * 1000

#time attributes back to the layers
time_steps <- time(climate_raw)[1:204]
time(temp_celsius) <- time_steps
time(precip_mm)    <- time_steps
names(temp_celsius) <- as.character(time(temp_celsius))

#parameters
radii <- c(1000)
unique_years <- sort(unique(sites_df$Sample_start_year[sites_df$Sample_start_year >= 2000]))
raster_years <- rep(2000:2016, each = 12)
raster_vars  <- c(rep("t2m", 204), rep("tp", 204))

#creating clean tibble to store the data
climate_results <- tibble(SSBS = sites_df$SSBS)

#creating columns to store the results
for (r in radii) {
  climate_results[[paste0("climate_mean_temp_c_", r, "m")]]   <- NA_real_
  climate_results[[paste0("climate_annual_precip_mm_", r, "m")]] <- NA_real_
}

#points kept in a metric CRS so buffer() widths are accurate in meters
sites_spatial <- vect(sites_df, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")
sites_spatial_metric <- project(sites_spatial, "EPSG:3035")

#we will need to know the days in month to convert monthly precipitation to annual precipitation 
days_in_month <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

#run
for (yr in unique_years) { #given the amoount of data we go by years
  
  cat("PROCESSING CLIMATE DATA FOR SAMPLING YEAR:", yr, "\n")
  site_indices <- which(sites_df$Sample_start_year == yr)
  if (length(site_indices) == 0) next #if there is no site for this year we skip it
  
  #isolate the 12 monthly layers belonging to the target year
  temp_layers   <- climate_raw[[which(raster_years == yr & raster_vars == "t2m")]]
  precip_layers <- climate_raw[[which(raster_years == yr & raster_vars == "tp")]]
  
  #estimate the indicators 
  annual_mean_temp    <- mean(temp_layers) - 273.15
  precip_monthly_total <- precip_layers * days_in_month
  annual_total_precip <- sum(precip_monthly_total) * 1000
  
  #clean non-finite values on the raster raster, before any extraction
  annual_total_precip[is.nan(annual_total_precip)] <- NA
  annual_total_precip[is.infinite(annual_total_precip)] <- NA
  
  current_points <- sites_spatial_metric[site_indices, ]
  
  #iterative loop over your radii spatial scales
  for (r in radii) {
    site_polygons_metric <- buffer(current_points, width = r)  #buffers built in metric CRS, then projected into the raster's native CRS
    site_polygons_native <- project(site_polygons_metric, crs(annual_mean_temp))
    
    extracted_temp   <- terra::extract(annual_mean_temp,    site_polygons_native, fun = mean, na.rm = TRUE, exact = TRUE)[, 2] #extraction
    extracted_precip <- terra::extract(annual_total_precip, site_polygons_native, fun = mean, na.rm = TRUE, exact = TRUE)[, 2] #extraction
    
    #select column
    col_temp   <- paste0("climate_mean_temp_c_", r, "m")
    col_precip <- paste0("climate_annual_precip_mm_", r, "m")
    
    climate_results[[col_temp]][site_indices] <- extracted_temp #save values in columns
    climate_results[[col_precip]][site_indices] <- extracted_precip #save values in columns
  }
}

#save results
saveRDS(climate_results, "Intermediate_dataset/fixed_effects_intermediate/fixed_effect_climate.rds")


#There are Nan values that need to be fix. These correspond to coastal points
#We run them again with the land cell values only

#list of failed sites
nan_indices_global <- which(is.na(climate_results$climate_mean_temp_c_1000m))
sites_patch <- sites_df[nan_indices_global, ]

#projections
sites_spatial_patch  <- vect(sites_patch, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")
sites_spatial_metric <- project(sites_spatial_patch, "EPSG:3035")

#run the loop again for the patch sites

unique_patch_years <- sort(unique(sites_patch$Sample_start_year[sites_patch$Sample_start_year >= 2000])) #defining the years to run

for (yr in unique_patch_years) {
  #row positions relative to this patch subset
  sub_indices <- which(sites_patch$Sample_start_year == yr)
  if (length(sub_indices) == 0) next #we skip the year if there are no values
  
  #annual aggregate layers
  annual_mean_temp    <- mean(temp_celsius[[which(raster_years == yr)]])
  precip_year_stack    <- subset(climate_raw, precip_indices)[[which(raster_years == yr)]]  # raw, meters
  precip_monthly_total <- precip_year_stack * days_in_month
  annual_total_precip  <- sum(precip_monthly_total) * 1000
  
  #clean non-finite values on original raster before any extraction
  annual_total_precip[is.nan(annual_total_precip)] <- NA
  annual_total_precip[is.infinite(annual_total_precip)] <- NA
  
  current_points_metric <- sites_spatial_metric[sub_indices, ]
  
  for (r in radii) {
    #given the coastal points are falling into grids of approax 9x9 km that are "water" and therefore NaN
    #I need a safe buffer from which to extract the values that is big enough to be able to gather info from the nearest pixel
    #I will use a 12 km buffer for that
    rescue_buffer_metric <- buffer(current_points_metric, width = 12000)
    rescue_buffer_native  <- project(rescue_buffer_metric, crs(annual_mean_temp))
    
    extracted_temp   <- terra::extract(annual_mean_temp,    rescue_buffer_native, fun = mean, na.rm = TRUE, exact = TRUE)[, 2]
    extracted_precip <- terra::extract(annual_total_precip, rescue_buffer_native, fun = mean, na.rm = TRUE, exact = TRUE)[, 2]
    
    #sub-indices back to the true absolute row indices in climate_results
    global_target_indices <- nan_indices_global[sub_indices]
    
    #selecting the column to modify
    col_temp   <- paste0("climate_mean_temp_c_", r, "m")
    col_precip <- paste0("climate_annual_precip_mm_", r, "m")
    
    #overwrite the NaNs with the safe land values
    climate_results[[col_temp]][global_target_indices]   <- extracted_temp
    climate_results[[col_precip]][global_target_indices] <- extracted_precip
  }
}

#save final patched results
saveRDS(climate_results, "Intermediate_dataset/fixed_effects_intermediate/fixed_effect_climate.rds")
  

#JOIN FIXED EFFECTS ------------------------------------------------------------
cropland_data <- readRDS("Intermediate_dataset/fixed_effects_intermediate/fixed_effect_cropland.rds") 
wow_data <- readRDS("Intermediate_dataset/fixed_effects_intermediate/fixed_effect_wow.rds") 
urban_data <- readRDS("Intermediate_dataset/fixed_effects_intermediate/fixed_effect_urban.rds")
veg_data <- readRDS("Intermediate_dataset/fixed_effects_intermediate/fixed_effect_natural_veg.rds")
pasture_data <- readRDS("Intermediate_dataset/fixed_effects_intermediate/fixed_effect_pastures.rds")
sea_data <- readRDS("Intermediate_dataset/fixed_effects_intermediate/fixed_effect_sea.rds")
shannon_data <- readRDS("Intermediate_dataset/fixed_effects_intermediate/fixed_effect_shannon.rds")
alpha_data <- readRDS("Intermediate_dataset/fixed_effects_intermediate/fixed_effect_jrc.rds") 
pesticide_data <- readRDS("Intermediate_dataset/fixed_effects_intermediate/fixed_effect_pesticides.rds")
field_size_data <- readRDS("Intermediate_dataset/fixed_effects_intermediate/fixed_effect_field_size.rds")
fertilization_data <- readRDS("Intermediate_dataset/fixed_effects_intermediate/fixed_effect_field_fertilization.rds") 
climate_data <- readRDS("Intermediate_dataset/fixed_effects_intermediate/fixed_effect_climate.rds")

sites_df_final <- sites_df %>%
  left_join(cropland_data, by = "SSBS")  %>%
  left_join(urban_data,    by = "SSBS") %>%
  left_join(wow_data,      by = "SSBS")  %>%
  left_join(sea_data,      by = "SSBS")  %>%
  left_join(veg_data,      by = "SSBS")  %>%
  left_join(pasture_data,  by = "SSBS")  %>%
  left_join(shannon_data,  by = "SSBS")  %>%
  left_join(alpha_data,    by = "SSBS")  %>%
  left_join(pesticide_data, by = "SSBS")  %>%
  left_join(field_size_data, by = "SSBS")  %>%
  left_join(fertilization_data, by = "SSBS")  %>%
  left_join(climate_data, by = "SSBS")

#WATER DENOMINATOR CORRECTION ------------------------------------------------------------
#the water indicator needs to be corrected 
#in the water_mask reclassification, sea water is mapped NA (not 0), whily dry land is 0 and the others 1. 
#when I did fun=mean and na.rm=TRUE that drops the NA (sea) cells for the average  so the raw water_pct_1000m is not the percentage of the wholl 1000m buffer but the percerntage of the land portion
#so later, when we analyse the total land use cover, we would get wrong results when doing the sum of water + sea water percentages. 
sites_df_final <- sites_df_final %>%
  mutate(water_pct_1000m_landfrac = water_pct_1000m, #keep the original land-fraction values for traceability
    water_pct_1000m = water_pct_1000m_landfrac * (1 - sea_pct_1000m / 100)) #overwrite with the whole-buffer fraction


#LAND USE REESCALING ------------------------------------------------------------
#we used various sources to obtain land cover in the landscape. Which means that in some cases we can have over 100% coverage. So we need to correct that

sites_df_final<- sites_df_final %>%  mutate(total_landuse_1000m = urban_pct_1000m + crop_pct_1000m + water_pct_1000m + sea_pct_1000m + natural_pct_1000m + pasture_pct_1000m) #total land use
sites_df_final <- sites_df_final %>% mutate(available_land_1000m = pmax(0, 100 - crop_pct_1000m - sea_pct_1000m)) #if the total_land_use is more than 100%, we estimate how much land is left after both crop_pct (Liao) and sea_pct are held fixed as anchors and not rescaled (we don´t expect sea to have change or present many errors)
sites_df_final <- sites_df_final %>%  mutate(crop_pct_1000m_mod = if_else(crop_pct_1000m + sea_pct_1000m > 100, pmax(0, 100 - sea_pct_1000m), crop_pct_1000m)) #where crop_pct + sea_pct alone exceed 100%, I trust sea_pct and cap crop_pct to whatever room remains
sites_df_final <- sites_df_final %>% mutate(available_land_1000m = pmax(0, 100 - crop_pct_1000m_mod - sea_pct_1000m)) #recompute available_land using the capped crop value instead of raw crop_pct

#now, when total_landuse_ is more than 100% we rescale urban/water/natural/pasture (no crop, no sea) to fit within available_land
sites_df_final <- sites_df_final %>%
  mutate(urban_pct_1000m_mod = if_else(total_landuse_1000m > 100 & (total_landuse_1000m - crop_pct_1000m_mod - sea_pct_1000m) > 0,
                                  (urban_pct_1000m / (total_landuse_1000m - crop_pct_1000m_mod - sea_pct_1000m)) * available_land_1000m, urban_pct_1000m),
    water_pct_1000m_mod = if_else(total_landuse_1000m > 100 & (total_landuse_1000m - crop_pct_1000m_mod - sea_pct_1000m) > 0,
                                  (water_pct_1000m / (total_landuse_1000m - crop_pct_1000m_mod - sea_pct_1000m)) * available_land_1000m, water_pct_1000m),
    natural_pct_1000m_mod = if_else(total_landuse_1000m > 100 & (total_landuse_1000m - crop_pct_1000m_mod - sea_pct_1000m) > 0,
                                    (natural_pct_1000m / (total_landuse_1000m - crop_pct_1000m_mod - sea_pct_1000m)) * available_land_1000m, natural_pct_1000m),
    pasture_pct_1000m_mod = if_else(total_landuse_1000m > 100 & (total_landuse_1000m - crop_pct_1000m_mod - sea_pct_1000m) > 0,
                                    (pasture_pct_1000m / (total_landuse_1000m - crop_pct_1000m_mod - sea_pct_1000m)) * available_land_1000m, pasture_pct_1000m) )

sites_df_final<- sites_df_final %>% mutate(total_landuse_1000m_mod = urban_pct_1000m_mod + crop_pct_1000m_mod + water_pct_1000m_mod + sea_pct_1000m + natural_pct_1000m_mod + pasture_pct_1000m_mod)
sum(sites_df_final$total_landuse_1000m_mod > 100.001, na.rm = TRUE) #should be zero


#SAVE ------------------------------------------------------------
write.csv(sites_df_final, "Intermediate_dataset/sites_df_predictors.csv", row.names = FALSE)




