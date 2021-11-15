
# Monthy NDVI Phenology for PlanetScope Imagery
#   Generates NDVI Phenology Images for PlanetScope Tiles
#   Assumes tiles downloaded and formatted with naming convention: EEEEEEE_NNNNNNN_****_SR_clip.tif (for surface reflectance)
#      where EEEEEEE is easting and NNNNNNN is northing for each tile
#   Filters NDVI by Unusable Data Mask from PlanetScope (file naming convention EEEEEEE_NNNNNNN_****_SR_clip.tif)
#   Expands upon outputs from 'planet_phenology_basic.R' by providing not just monthly stats (e.g. mean, stdev)
#      but a resampled timeseries of values using actual dates within the month to estimate the phenology curve
#   Currently, just takes nearest neighbor values of each temporal point in the timeseries and assigns a new value
#      as an inverse-distance weighted average 

library(tidyverse)
library(here)
library(raster)

# ******************************************************
# Load all information on number of scenes, metadata, etc.
# ******************************************************

# Directory containing all planetscope download subfolders
target_directory <- "D:/SERDP/Mission_Creek/Planet/"

# Load Metadata
source(here::here("planet_metadata_review.R"))

# for now, filter out rasters which erroneously lost their location in the name...
all_metadata <- metadata
metadata <- metadata %>% 
  filter((substr(metadata$file, 1,1) != "_"))

# Get coordinates for each tile (from filenames)
coords <- str_extract_all(metadata$file, "[\\d]+")
nth_element <- function(list, n)
{
  return(list[[n]])
}
metadata$northing <- as.numeric(lapply(coords, nth_element, n=1))
metadata$easting <- as.numeric(lapply(coords, nth_element, n=2))
# Get filenames for data raster and for unusable data mask 
metadata$raster_file <- paste(sub("metadata_clip\\.xml", "", metadata$file), "SR_clip.tif", sep="")
metadata$quality_file <- paste(sub("metadata_clip\\.xml", "", metadata$file), "DN_udm_clip.tif", sep="")

# Statistics across each tile, by month
scene_data <- metadata %>% 
  group_by(easting, northing, month) %>%
  summarize(number_of_images = n())
# Statistics across each tile, including all months
tile_data <- metadata %>% 
  group_by(easting, northing, month) %>%
  summarize(number_of_images = n())
number_total_tiles <- length(unique(paste(scene_data$easting, scene_data$northing)))


# ******************************************************
# Check image quality - mask clouds, pixels with malfunctioned sensors, etc. 
# ******************************************************
# Function to check whether NDVI values will be bad, based on Usable Data Mask from the Planetscope standard. 
#   Specification: https://www.planet.com/products/satellite-imagery/files/1610.06_Spec%20Sheet_Combined_Imagery_Product_Letter_ENGv1.pdf
#   8-bit image, with each bit representing a different kind of problem
#   Bit 0 --> backfill (if 1)
#   Bit 1 --> clouds (if 1) - NOTE this will miss small clouds, miss haze, and identify snow as clouds
#   Bit 4 --> errors in Red band (if 1)
#   Bit 6 --> errors in the Red Edge band (if 1)
#     Other bits are for other bands (not relevant for NDVI)
isPixelBad <- function(quality_value)
{
  quality_bits <- as.numeric(intToBits(quality_value)[1:8])
  return((quality_bits[[1]]+quality_bits[[2]]+quality_bits[[5]]+quality_bits[[7]]))
}

# ******************************************************
# Functions related to inverse distance resampling on raster
# ******************************************************

# Get nearest temporal neighbors to a target date-of-year among input scene dates
get_point_neighbors <- function(target_day, data, num_neighbors)
{
  data$distance <- abs(target_day - data$day_of_year)
  data$index <- 1:nrow(data)
  data <- data %>% arrange(distance)
  if(num_neighbors > nrow(data))
  {
    print(paste("WARNING - requesting more neighbors than total number of points... Requested ", num_neighbors, " but only ", nrow(data), " rows in original timeseries. Keeping all points.", sep=""))
    return(data)
  }
  return(data[1:num_neighbors,])
}

# Get inverse time-distance between target date and all neighbor dates
inverse_distance_coefficients <- function(target_day, neighbors)
{
  return(list(target_day, 
              1/(7+neighbors$distance), 
              neighbors$index))
}

# Apply both above functions (get neighbors and then get inverse distance coefficients)
get_resampling_matrix <- function(target_day, data, num_neighbors)
{
  neighbors <- get_point_neighbors(target_day, data, num_neighbors)
  return(inverse_distance_coefficients(target_day, neighbors))
}

# Resample original rasters onto new date
temporally_resample_raster <- function(resampling_matrix, phenoseries_in)
{
  # Temporal neighbor scenes for analysis
  neighbor_scenes <- phenoseries_in[[resampling_matrix[[3]]]]
  # Resample based on temporal inverse distance
  phenoseries_out <- sum(neighbor_scenes*resampling_matrix[[2]], na.rm=TRUE) / sum(resampling_matrix[[2]])
}


# ******************************************************
# Actually reprocess all input scenes to new temporal space
# ******************************************************

# Iterate across all tiles, generate resampled NDVI phenoseries for each individual scene, and write outputs to disk 
tiles_processed <- 0
for(current_northing in unique(metadata$northing))
{
  for(current_easting in unique(metadata$easting))
  {
    for(current_year in unique(metadata$year))
    {
    # ********* Load and generate masked NDVI scenes ***********
    
    # Get all scenes within the current tile (unique Northing/Easting values)
    all_current_scenes <- metadata %>% 
      filter(northing == current_northing) %>%
      filter(easting == current_easting) %>%
      filter(year == current_year) %>%
      arrange(doy)
    # If this combination of northing and easting doesn't exist, continue...
    #   This occurs if the whole rectangular extent isn't occupied by tiles
    if(nrow(all_current_scenes) < 1)
      next
    
    print(paste("For tile at ", current_northing, "x", current_easting, " in year ", current_year, ", we found ", nrow(all_current_scenes), " scenes.", sep=""))
    
    # Assemble all rasters for this scene, and generate NDVI (masking out bad pixels) 
    all_rasters <- raster()
    all_rasters_evi <- raster()
    all_metadata <- all_current_scenes[0,]
    all_quality_images <- raster()
    all_bad_data_masks <- raster()
    for(ind in 1:nrow(all_current_scenes))
    {
      # Verify that image exists and is not empty
      file_size <- file.info(paste(target_directory, all_current_scenes[ind,]$raster_file, sep=""))$size
      if(is.na(file_size)) 
        next
      if(file_size < 1000000)
        next
      
      # Generate NDVI
      image <- stack(paste(target_directory, all_current_scenes[ind,]$raster_file, sep=""))
      ndvi <- (image[[4]] - image[[3]]) / (image[[4]] + image[[3]])
      # NOTE EVI formula below is modified because PlanetScope imagery comes with reflectance scaled up by 10000
      evi <- (2.5*(image[[4]]-image[[3]])) / (10000 + image[[4]] + 6*image[[3]] - 7.5*image[[1]])
      
      # For some reason, some of these quality rasters can't be loaded in R... skip scenes when that's the case
      if(!tryCatch(
        expr = {
          raster(paste(target_directory, all_current_scenes[ind,]$quality_file, sep="")); TRUE
        }, 
        error = function(e) {
          print(paste("   There's an issue with the quality image ", all_current_scenes[ind,]$quality_file, sep=""))
          return(FALSE) 
        }))
        next 
      
      # Generate Quality Mask
      quality <- raster(paste(target_directory, all_current_scenes[ind,]$quality_file, sep=""))
      bad_data_vector <- as.numeric(lapply(getValues(quality), isPixelBad))
      bad_data_mask <- quality
      bad_data_mask[,] <- bad_data_vector 
      
      # Output Raster 
      all_rasters <- addLayer(all_rasters, mask(ndvi, bad_data_mask, maskvalue=1))
      all_rasters_evi <- addLayer(all_rasters_evi, mask(evi, bad_data_mask, maskvalue=1))
      all_metadata <- rbind(all_metadata, all_current_scenes[ind,])
      all_quality_images <- addLayer(all_quality_images, quality)
      all_bad_data_masks <- addLayer(all_bad_data_masks, bad_data_mask)
    }
    print(paste("   Generated ", nrow(all_metadata), " NDVI scenes.", sep=""))
    if(nrow(all_metadata) == 0)
    {
      print(paste("For the current scene, NO output NDVI scenes were generated, probably because all the quality images are bad. "))
      next
    }
    
    # Find inverse distance coefficients for all target dates 
    source_dates <- data.frame(day_of_year = all_metadata$doy)
    target_dates <- (1:26)*14
    resampling_matrix_list <- lapply(target_dates, get_resampling_matrix, data=source_dates, num_neighbors=nrow(source_dates)/4)
    
    # Resample all rasters for current target scene
    final_phenoseries <- lapply(resampling_matrix_list, temporally_resample_raster, phenoseries_in=all_rasters)
    final_phenoseries_evi <- lapply(resampling_matrix_list, temporally_resample_raster, phenoseries_in=all_rasters_evi)
    # Collapse list of raster layers to a raster stack
    final_phenoseries <- stack(final_phenoseries)
    final_phenoseries_evi <- stack(final_phenoseries_evi)
    
    writeRaster(final_phenoseries, 
                paste(target_directory, "phenoseries/", current_easting, "_", current_northing, "_", current_year, "_inverse_distance.tif", sep=""), 
                overwrite=TRUE)
    writeRaster(final_phenoseries_evi, 
                paste(target_directory, "phenoseries/", current_easting, "_", current_northing, "_", current_year, "_inverse_distance_evi.tif", sep=""), 
                overwrite=TRUE)
    
    tiles_processed <- tiles_processed + 1
    print(paste("   Saved resampled rasters! Have finished ", tiles_processed, " tiles out of ", number_total_tiles*length(unique(metadata$year)), ", or ", round(tiles_processed/number_total_tiles*100,2), "%.", sep=""))
    }
  }
}

