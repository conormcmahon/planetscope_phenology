
library(tidyverse)
library(here)
library(raster)

target_directory <- "D:/SERDP/Huachuca/Planet/output_imagery/"

source(here::here("planet_metadata_review.R"))

# for now, filter out rasters which erroneously lost their location in the name...
all_metadata <- metadata
metadata <- metadata %>% 
  filter((substr(metadata$file, 1,1) != "_"))

coords <- str_extract_all(metadata$file, "[\\d]+")
nth_element <- function(list, n)
{
  return(list[[n]])
}
metadata$northing <- as.numeric(lapply(coords, nth_element, n=1))
metadata$easting <- as.numeric(lapply(coords, nth_element, n=2))
metadata$raster_file <- paste(sub("metadata_clip\\.xml", "", metadata$file), "SR_clip.tif", sep="")
metadata$quality_file <- paste(sub("metadata_clip\\.xml", "", metadata$file), "DN_udm_clip.tif", sep="")

scene_data <- metadata %>% 
  group_by(easting, northing, month) %>%
  summarize(number_of_images = n())
tile_data <- metadata %>% 
  group_by(easting, northing, month) %>%
  summarize(number_of_images = n())
number_total_tiles <- length(unique(paste(scene_data$easting, scene_data$northing)))


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


# Iterate across all tiles 
tiles_processed <- 0
for(current_northing in unique(metadata$northing))
{
  for(current_easting in unique(metadata$easting))
  {
      
    all_current_scenes <- metadata %>% 
      filter(northing == current_northing) %>%
      filter(easting == current_easting) %>%
      arrange(doy)
    # If this combination of northing and easting doesn't exist, continue...
    if(nrow(all_current_scenes) < 1)
      next
    
    print(paste("For tile at ", current_northing, "x", current_easting, ", we found ", nrow(all_current_scenes), " scenes.", sep=""))
    
    # Assemble all rasters for this scene, and generate NDVI (masking out bad pixels) 
    all_rasters <- raster()
    all_metadata <- all_current_scenes[0,]
    for(ind in 1:length(all_current_scenes))
    {
      # Verify that image exists and is not empty
      file_size <- file.info(paste(target_directory, all_current_scenes[ind,]$raster_file, sep=""))$size
      if(is.na(file_size)) 
        next
      if(file_size < 1000000)
        next
      
      image <- stack(paste(target_directory, all_current_scenes[ind,]$raster_file, sep=""))
      ndvi <- (image[[4]] - image[[3]]) / (image[[4]] + image[[3]])

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
      
      quality <- raster(paste(target_directory, all_current_scenes[ind,]$quality_file, sep=""))
      bad_data_vector <- as.numeric(lapply(getValues(quality), isPixelBad))
      bad_data_mask <- quality
      bad_data_mask[,] <- bad_data_vector
      
      all_rasters <- addLayer(all_rasters, mask(ndvi, bad_data_mask, maskvalue=1))
      all_metadata <- rbind(all_metadata, all_current_scenes[ind,])
    }
    print(paste("   Generated ", nrow(all_metadata), " NDVI scenes.", sep=""))
    if(nrow(all_metadata) == 0)
    {
      print(paste("For the current scene, NO output NDVI scenes were generated, probably because all the quality images are bad. "))
      next
    }
    
    # Generate monthly statistics for scene
    months <- unique(all_metadata$month)
    monthly_rasters_max <- raster()                 # Maximum NDVI in each pixel
    monthly_rasters_mean <- raster()                # Mean NDVI in each pixel
    monthly_rasters_min <- raster()                 # Minimum NDVI in each pixel
    monthly_rasters_sd <- raster()                  # Standard deviation of NDVI in each pixel
    monthly_rasters_n_images <- all_rasters[[1]]*0  # How many good data points are there for each pixel?
    for(month in months)
    {
      current_month <- all_rasters[[which(all_metadata$month == month)]]
      monthly_rasters_max <- addLayer(monthly_rasters_max, max(current_month, na.rm=TRUE))
      monthly_rasters_mean <- addLayer(monthly_rasters_mean, mean(current_month, na.rm=TRUE))
      monthly_rasters_min <- addLayer(monthly_rasters_min, min(current_month, na.rm=TRUE))
      monthly_rasters_sd <- addLayer(monthly_rasters_sd, calc(current_month, fun=function(x){sd(x, na.rm=TRUE)}))
      monthly_rasters_n_images <- monthly_rasters_n_images + !is.na(current_month)
    }
    print(paste("   Generated all max, min, mean, sd, and n values."))
    
    # Output stats rasters 
    writeRaster(monthly_rasters_max, 
                paste(target_directory, "../phenoseries/", all_metadata[1,]$easting, "_", all_metadata[1,]$northing, "_max.tif", sep=""), 
                overwrite=TRUE)
    writeRaster(monthly_rasters_mean, 
                paste(target_directory, "../phenoseries/", all_metadata[1,]$easting, "_", all_metadata[1,]$northing, "_mean.tif", sep=""), 
                overwrite=TRUE)
    writeRaster(monthly_rasters_min, 
                paste(target_directory, "../phenoseries/", all_metadata[1,]$easting, "_", all_metadata[1,]$northing, "_min.tif", sep=""), 
                overwrite=TRUE)
    writeRaster(monthly_rasters_sd, 
                paste(target_directory, "../phenoseries/", all_metadata[1,]$easting, "_", all_metadata[1,]$northing, "_sd.tif", sep=""), 
                overwrite=TRUE)
    writeRaster(monthly_rasters_n_images, 
                paste(target_directory, "../phenoseries/", all_metadata[1,]$easting, "_", all_metadata[1,]$northing, "_n_images.tif", sep=""), 
                overwrite=TRUE)
    
    tiles_processed <- tiles_processed + 1
    print(paste("   Saved all output rasters! Have finished ", tiles_processed, " tiles out of ", number_total_tiles, ", or ", round(tiles_processed/number_total_tiles*100,2), "%.", sep=""))
  }
}