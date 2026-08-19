#---------------------------------------------------------------------------------------------------------------------------------------#
# CODE DETAILS                                                                                                                          #
#                                                                                                                                       #
# Author: Catalina Vattuone  (cvattuonet@gmail.com)                                                                                     #        
# Date of latest update: 19-08-2026                                                                                                     #
# Main use: Study on "Assessing how agricultural practices and landscape composition shape differences in arthropod biodiversity among  #  
#           habitat types in Europe using the PREDICTS database". Internship at INRAE                                                   #
#                                                                                                                                       #
# Content: Preparation of agronomic databases for estimating landscape predictors                                                       #                           
#---------------------------------------------------------------------------------------------------------------------------------------#


#LIBRARIES  ------------------------------------------------------------
pacman::p_load(terra, dplyr)

#PESTICIDES  ------------------------------------------------------------
#The source comes with a folder with individual files for individual active ingredients. 
#Here we prepare a single raster layer with the total pesticide application for the medium scenario

#Directory
raster_folder <- "Databases/Pesticides_Porta/maps_geotiff"

#Create the list of names of all teh files of medium scenario in the folder
pesticide_files <- list.files(path = raster_folder, pattern = "_M\\.tif$",  full.names = TRUE)

#Load all the files to R as a SpatRaster collection 
pesticide_stack <- rast(pesticide_files) 

#Transform the -1 values to NA values so that when the usm is performed there is no problems
pesticide_stack[pesticide_stack == -1] <- NA

#Sum all layers into a single raster layer
total_pesticide_M <- sum(pesticide_stack, na.rm = TRUE)

#Save
dir.create("Intermediate_dataset/Pesticides", showWarnings = FALSE)
writeRaster(total_pesticide_M, "Intermediate_dataset/Pesticides/total_pesticide_medium_scenario.tif", overwrite = TRUE) 

#N FERTILIZER  ------------------------------------------------------------
#NPKGRIDS is used as source, but each crop's Nrate is expressed per hectare the crop
#Therefore, CROPGRIDS is also used following the formulas below
# total_n_kg= sum_i (Nrate_i * harvarea_i) , this gives the kg N in the cell
# n_areal_load= total_n_kg / cell_area_ha, this gives the kg N per ha of land
# n_mean_rate= total_n_kg / sum_i(harvarea_i), this gives the kg N per ha of cropped land

#Directories
npk_dir  <- "Databases/Fertilizers_NPKGRIDS/"
crop_dir <- "Databases/CROPGRID/"

#Variables of interest
npk_var  <- "Nrate"
crop_var <- "harvarea"   #not the croparea because NPKGRIDS rates are per ha of harvested area

#limits to cut the rasters to Europe
europe_bounding_box <- ext(-25, 45, 34, 72)

#we create a list of the files for each case, considering other files could be in the folders
npk_files  <- list.files(npk_dir,  pattern = "^NPKGRIDS.*\\.nc$",  full.names = TRUE)
crop_files <- list.files(crop_dir, pattern = "^CROPGRIDS.*\\.nc$", full.names = TRUE)

#we clean the files names to obtain in each case the name of the crop that is being represented. 
get_crop_name <- function(paths) tolower(sub("^[^_]*_", "", sub("\\.nc$", "", basename(paths))))

npk_crops  <- get_crop_name(npk_files)
crop_crops <- get_crop_name(crop_files)

#creating a crop matching inventory to check every crop is in both datasets
crop_inventory <- data.frame(crop = sort(union(npk_crops, crop_crops))) %>%
  mutate(in_npkgrids  = crop %in% npk_crops,
         in_cropgrids = crop %in% crop_crops,
         matched      = in_npkgrids & in_cropgrids)

matched_crops <- crop_inventory$crop[crop_inventory$matched]

#creating the template raster to accumulate the results of the process
template <- crop(rast(npk_files[match(matched_crops[1], npk_crops)], subds = npk_var), europe_bounding_box)
total_n_kg    <- template; total_n_kg[]    <- 0
total_area_ha <- template; total_area_ha[] <- 0

#creating a log to keep track of the processing status of each crop
processing_log <- data.frame(crop = character(), status = character())

#a loop that goes by each of the crops and updates the total_n_kg and total_area_ha rasters. 
for (i in seq_along(matched_crops)) {
  cname     <- matched_crops[i] #the name of the crop being processed
  npk_path  <- npk_files[match(cname, npk_crops)] #selects the NPKGRIDS file for the crop
  crop_path <- crop_files[match(cname, crop_crops)] #selects the CROPGRIDS file for the crop
  
  status <- tryCatch({ 
    r_rate <- crop(rast(npk_path,  subds = npk_var),  europe_bounding_box) #extracts the N rate for the crop
    r_area <- crop(rast(crop_path, subds = crop_var), europe_bounding_box) #extracts the harvested area for the crop
    
    r_rate[r_rate == -1 | is.na(r_rate)] <- 0  #clean no data flags (NPKGRIDS marks ocean as -1)
    r_area[r_area < 0   | is.na(r_area)] <- 0  #clean no data flags (CROPGRIDS marks ocean as -1)
    
    if (!compareGeom(r_rate, r_area, stopOnError = FALSE)) { #if the rasters are not aligned, resample the area raster to match the rate raster
      r_area <- resample(r_area, r_rate, method = "bilinear") #resample the area raster to match the rate raster
      r_area[r_area < 0 | is.na(r_area)] <- 0 
    }
    
    total_n_kg    <- total_n_kg + (r_rate * r_area) #computes the total N in kg for the crop and adds it to the total N raster
    total_area_ha <- total_area_ha + r_area #computes the total harvested area in ha for the crop and adds it to the total area raster
    "OK"
  }, error = function(e) paste("FAILED:", conditionMessage(e))) #shows the error message if the process fails for a crop
  
  cat(status, "\n")
  processing_log <- rbind(processing_log, data.frame(crop = cname, status = status)) #saves the status of the processing for the crop in the log
  if (i %% 25 == 0) gc() #calls garbage collection (gc) every 25 iterations to free up memory
}

#defines the area of a cell in hectares, which is used to estimate the areal loading of N on land
cell_area_ha <- cellSize(total_n_kg, unit = "ha")  

#estimates the areal loading of N on land, which is the total N in kg divided by the area of the cell in hectares
n_areal_load <- total_n_kg / cell_area_ha

#estimates the mean rate of N application on cropped land
n_mean_rate  <- total_n_kg / total_area_ha

#We set the areal loading and mean rate to NA where the total area of cropped land is zero to indicate that there is no cropped land in those cells.
n_mean_rate[total_area_ha == 0] <- NA             

#Saving all rasters
dir.create("Intermediate_dataset/Fertilizers_NPKGRIDS", showWarnings = FALSE)
writeRaster(n_areal_load, "Intermediate_dataset/Fertilizers_NPKGRIDS/n_areal_load_kg_per_ha_land.tif",  overwrite = TRUE)
writeRaster(n_mean_rate, "Intermediate_dataset/Fertilizers_NPKGRIDS/n_mean_rate_kg_per_ha_crop.tif",   overwrite = TRUE)
writeRaster(total_area_ha, "Intermediate_dataset/Fertilizers_NPKGRIDS/total_crop_harvested_area_ha.tif", overwrite = TRUE)


#SHANNON INDEX ------------------------------------------------------------
#load
crop_dir <- "Databases/CROPGRID/"
crop_var <- "harvarea"
crop_files <- list.files(crop_dir, pattern = "^CROPGRIDS.*\\.nc$", full.names = TRUE)

crop_rasters <- lapply(crop_files, function(x){rast(x, subds = crop_var)})  #opens each file and pulls just the harvested area set
all_crops_raster <- rast(crop_rasters) #stak all separate single layer rasers into one multi-layer SpatRaster

shannon_index <- function(p) {
  p <- p[!is.na(p) & p > 0]  
  if (length(p) == 0) {
    return(NA) 
  }
  p <- p / sum(p) 
  return(-sum(p * log(p)))
}

shannon_map <- terra::app(all_crops_raster, shannon_index) #function to each pixel

dir.create("Intermediate_dataset/Crop_Diversity_CROPGRIDS", showWarnings = FALSE)
writeRaster(shannon_map, "Intermediate_dataset/Crop_Diversity_CROPGRIDS/shannon_map.tif", overwrite = TRUE)
