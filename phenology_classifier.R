
library(tidyverse)
library(raster)
library(rgdal)
library(here)
library(caret)
library(snow)
library(randomForest)
library(rlist)
library(rgeos)
library(exactextractr)


set.seed(7)

# ***** Load Data for Model ****

# Load Canopy Height Information
chm <- raster("D:/SERDP/Walnut_Gulch/LiDAR/tiles/chm/chm_entire_watershed.tif")
# Load phenology rasters
#phenoseries <- stack("D:/SERDP/Walnut_Gulch/planet/phenoseries/phenoseries_mosaic.tif")
phenoseries <- stack("C:/Users/grad/Downloads/phenoseries_mosaic.tif")
phenoseries <- addLayer(phenoseries, chm)

# **** Load training polygons ****
segments <- readOGR("D:/SERDP/Walnut_Gulch/field_data/sampling_poly_wgs84_30.shp")
segments <- spTransform(segments, crs(phenoseries))
names(segments) <- c("id","class")
segments_new <- readOGR("D:/SERDP/Walnut_Gulch/field_data/sampling_poly_wgs84_30_new.shp")
segments_new <- spTransform(segments_new, crs(phenoseries))
names(segments_new) <- names(segments)
segments <- rbind(segments, segments_new)
segments$area <- area(segments)/900
as.data.frame(segments) %>% group_by(class) %>% summarize(total_area = sum(area)) %>% arrange(-total_area)
segments <- crop(segments, extent(phenoseries))

cropAndExtract <- function(segment, image, buffer)
{
  image_subset <- crop(image, extent(buffer(segment, width=buffer)))
  exact_extract(image_subset, segment, include_xy=TRUE)
}
segmentAllPolygons <- function(segments, image)
{
  data <- data.frame()
  for(ind in 1:nrow(segments))
  {
    new_data <- cropAndExtract(segments[ind,], image, 3)[[1]]
    new_data$segment_ind <- ind
    if(nrow(data) == 0)
      data <- new_data
    else
      data <- rbind(data, new_data)
    print(paste("Segmented ",ind," segments out of ", nrow(segments)," or about ",round(ind/nrow(segments),4)*100,"%. So far, extracted ",nrow(data)," cells.",sep=""))    
  }
  return(data)
}
  
# **** Set Up Training Data ****
# --- Cottonwood ---
print("")
print("Extracting cottonwood training data...")
cottonwood_segments <- segments[segments$class=="cottonwood",]
cottonwood_df <- segmentAllPolygons(cottonwood_segments, phenoseries)
cottonwood_df$class <- "Cottonwood"
# --- Mesquite ---
print("")
print("Extracting mesquite training data...")
mesquite_segments <- segments[segments$class=="mesquite",]
mesquite_df <- segmentAllPolygons(mesquite_segments, phenoseries)
mesquite_df$class <- "Mesquite"
# ---  ---
print("")
print("Extracting shrub training data...")
shrub_segments <- segments[segments$class=="shrub",]
shrub_df <- segmentAllPolygons(shrub_segments, phenoseries)
shrub_df$class <- "Shrub"
# --- Grass ---
print("")
print("Extracting grass training data...")
grass_segments <- segments[segments$class=="grass",]
grass_df <- segmentAllPolygons(grass_segments, phenoseries)
grass_df$class <- "Grass"
# --- All Label Data ---
label_df <- rbind(cottonwood_df, mesquite_df, shrub_df, grass_df)
label_df <- label_df %>% 
  drop_na()
label_df$class <- factor(label_df$class, 
                         levels = c("Cottonwood", "Mesquite", "Shrub", "Grass"))
names(label_df) <- c(as.character(1:26),"chm","x","y","coverage_fraction","segment_ind","class")
# Remove any shapefiles for which ALL data is 0 (these aren't inside the image extent)
label_df$ndvi_total <- rowSums(as.matrix(label_df[,1:10]))
label_df <- label_df %>% 
  filter(ndvi_total != 0) %>%
  filter(chm != 0)

# **** Visualize Greenness Phenology **** 
pheno_df_long <- label_df %>% pivot_longer((1:26), names_to="month", values_to="ndvi")
pheno_df_long$week <- as.numeric(pheno_df_long$month)*2
pheno_stats <- pheno_df_long %>% 
  group_by(class, week) %>%
  summarise(p_05=quantile(ndvi, probs=c(0.05)), 
            p_25=quantile(ndvi, probs=c(0.25)),
            p_50=quantile(ndvi, probs=c(0.50)),
            p_75=quantile(ndvi, probs=c(0.75)),
            p_95=quantile(ndvi, probs=c(0.95)))
# Natural Vegetation Phenology
veg_pheno_plot <- ggplot(data=pheno_stats[pheno_stats$class %in% factor(c("Cottonwood","Mesquite","Shrub","Grass")),], aes(x=week)) + 
  geom_line(aes(y=p_05), linetype="dashed", color="red") + 
  geom_line(aes(y=p_25), color = "cyan3") + 
  geom_line(aes(y=p_50)) + 
  geom_line(aes(y=p_75), color = "cyan3") + 
  geom_line(aes(y=p_95), linetype="dashed", color="red") + 
  geom_hline(yintercept=0, color="gray") +
  facet_wrap(~class, ncol=1) + 
  labs(x = "Week of Year",
       y = "NDVI") + 
  ggtitle("Greenness Phenology by Vegetation Type")
veg_pheno_plot_together <- ggplot(data=pheno_stats[pheno_stats$class %in% factor(c("Cottonwood","Mesquite","Shrub","Grass")),], aes(x=week, group=class)) + 
  geom_line(aes(y=p_50)) + 
  geom_hline(yintercept=0, color="gray") +
  labs(x = "Week of Year",
       y = "NDVI") + 
  ggtitle("Greenness Phenology by Vegetation Type")
# Plot Height Distributions
veg_height_plot <- ggplot(data=label_df) +
  geom_boxplot(aes(y=chm, group=class, col=class))
height_summary <- label_df %>%
  group_by(class) %>%
  summarize(height_mean = mean(chm),
            height_min = min(chm),
            height_max = max(chm),
            height_sd = sd(chm))

# Rename month columns to have text names (not numbers) 
names(label_df) <- c(paste("NDVI_week_",(1:26)*2,sep=""),"chm",
                     "x","y","coverage_fraction","segment_ind","class","ndvi_total")

# ------ Generate Training Data ------
min_valid_samples_by_class <- label_df %>%
  group_by(class) %>%
  summarize(num_samples = n())
num_samples <- 500 #min(min_valid_samples_by_class$num_samples) / 2
# --- Cottonwood ---
cottonwood_train_df <- label_df[label_df$class=="Cottonwood",]
cottonwood_train_indices <- sample(1:nrow(cottonwood_train_df),num_samples)
train <- cottonwood_train_df[cottonwood_train_indices,]
leftover <- cottonwood_train_df[-cottonwood_train_indices,]
# --- Mesquite ---
mesquite_train_df <- label_df[label_df$class=="Mesquite",]
mesquite_train_indices <- sample(1:nrow(mesquite_train_df),num_samples)
train <- rbind(train, mesquite_train_df[mesquite_train_indices,])
leftover <- rbind(leftover, mesquite_train_df[-mesquite_train_indices,])
# --- Shrub ---
shrub_train_df <- label_df[label_df$class=="Shrub",]
shrub_train_indices <- sample(1:nrow(shrub_train_df),num_samples)
train <- rbind(train, shrub_train_df[shrub_train_indices,])
leftover <- rbind(leftover, shrub_train_df[-shrub_train_indices,])
# --- Grass ---
grass_train_df <- label_df[label_df$class=="Grass",]
grass_train_indices <- sample(1:nrow(grass_train_df),num_samples)
train <- rbind(train, grass_train_df[grass_train_indices,])
leftover <- rbind(leftover, grass_train_df[-grass_train_indices,])


# Generate some validation data
#   to keep it balanced, we'll figure out the minimum number of remaining samples across all classes:
remaining_per_class <- leftover %>% group_by(class) %>% summarise(num=n())
min_remaining <- min(remaining_per_class$num)
# --- Cottonwood ---
cottonwood_leftover_df <- leftover[leftover$class=="Cottonwood",]
cottonwood_test_indices <- sample(1:nrow(cottonwood_leftover_df),min_remaining)
valid <- cottonwood_leftover_df[cottonwood_test_indices,]
leftover_final <- cottonwood_leftover_df[-cottonwood_test_indices,]
# --- Mesquite ---
mesquite_leftover_df <- leftover[leftover$class=="Mesquite",]
mesquite_test_indices <- sample(1:nrow(mesquite_leftover_df),min_remaining)
valid <- rbind(valid, mesquite_leftover_df[mesquite_test_indices,])
leftover_final <- rbind(leftover_final, mesquite_leftover_df[-mesquite_test_indices,])
# --- Shrub ---
shrub_leftover_df <- leftover[leftover$class=="Shrub",]
shrub_test_indices <- sample(1:nrow(shrub_leftover_df),min_remaining)
valid <- rbind(valid, shrub_leftover_df[shrub_test_indices,])
leftover_final <- rbind(leftover_final, shrub_leftover_df[-shrub_test_indices,])
# --- Grass ---
grass_leftover_df <- leftover[leftover$class=="Grass",]
grass_test_indices <- sample(1:nrow(grass_leftover_df),min_remaining)
valid <- rbind(valid, grass_leftover_df[grass_test_indices,])
leftover_final <- rbind(leftover_final, grass_leftover_df[-grass_test_indices,])



# Fit Random Forest Model
# Model based on ALL weeks - run initially to choose which weeks are most important, then run again with fewer weeks + LiDAR
#  modFit_rf <- train(as.factor(class) ~ NDVI_week_2 + NDVI_week_4 + NDVI_week_6 + NDVI_week_8 + NDVI_week_10 + NDVI_week_12 + NDVI_week_14 + NDVI_week_16 + NDVI_week_18 + NDVI_week_20 + NDVI_week_22 + NDVI_week_24 + NDVI_week_26 + NDVI_week_28 + NDVI_week_30 + NDVI_week_32 + NDVI_week_34 + NDVI_week_36 + NDVI_week_38 + NDVI_week_40 + NDVI_week_42 + NDVI_week_44 + NDVI_week_46 + NDVI_week_48 + NDVI_week_50 + NDVI_week_52, 
#                     method="rf", data = train)
# New model based on a subset of the 10 best discriminating weeks:
modFit_rf <- train(as.factor(class) ~ chm + NDVI_week_14 + NDVI_week_16 + NDVI_week_18 + NDVI_week_12 + NDVI_week_28 + NDVI_week_4 + NDVI_week_20 + NDVI_week_42 + NDVI_week_2 + NDVI_week_10,
                   method="rf", data = train)
# --- Validation  Model (Balanced) ---
validation_results <- raster::predict(modFit_rf, newdata=valid)
# Confusion Matrix with Validation Data
confusionMatrix(validation_results, as.factor(valid$class))
# --- Validation Model (All) ---
validation_results_all <- raster::predict(modFit_rf, newdata=leftover)
confusionMatrix(validation_results_all, as.factor(leftover$class))
# --- Variable Importance ---
varImp(modFit_rf)

# Get entire input image as a data frame for training
phenoseries_df <- as.data.frame(phenoseries, xy=TRUE)
names(phenoseries_df) <- c("x","y",paste("NDVI_week_",(1:26)*2,sep=""),"chm")
# Filter out locations with missing Planet or LiDAR data
phenoseries_df$ndvi_total <- rowSums(as.matrix(phenoseries_df[,3:28]))
phenoseries_df <- phenoseries_df %>% 
  drop_na() %>%
  filter(chm != 0) %>%
  filter(ndvi_total != 0) %>%
  select(1:(ncol(phenoseries_df)-1))
#predictions_all <- raster::predict(modFit_rf, newdata=phenoseries_df)
fit_size <- 100000
steps <- (1:floor(nrow(phenoseries_df)/fit_size))
predictions_all <- raster::predict(modFit_rf, newdata=phenoseries_df[1:(nrow(phenoseries_df) %% fit_size),])
for(ind in steps)
{
  new_predictions <- raster::predict(modFit_rf, newdata=phenoseries_df[((ind-1)*fit_size+1):(ind*fit_size) + (nrow(phenoseries_df) %% fit_size),])
  predictions_all <- c(predictions_all, new_predictions)
  print(paste("Finished a set of predictions... so far we've done ",length(predictions_all)," out of ",nrow(phenoseries_df),", or ",round(length(predictions_all)/nrow(phenoseries_df)*100,2),"%.",sep=""))
}
#predictions_all_num <- 1:length(predictions_all)
#predictions_all_num[predictions_all=="Cottonwood"] <- 1 
#predictions_all_num[predictions_all=="Mesquite"] <- 2
#predictions_all_num[predictions_all=="Shrub"] <- 3
#predictions_all_num[predictions_all=="Grass"] <- 4
predictions_df <- phenoseries_df %>% dplyr::select(1,2)
predictions_df$class <- predictions_all
predictions_raster <- rasterFromXYZ(predictions_df)
crs(predictions_raster) <- crs(phenoseries)
writeRaster(predictions_raster,
            paste("D:/SERDP/Walnut_Gulch/classes/classes_lidar_planet_10_weeks.tif",sep="")
            ,overwrite=TRUE)
