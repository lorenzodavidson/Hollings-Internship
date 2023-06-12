# Libraries ####

suppressPackageStartupMessages(suppressWarnings({
  #library(devtools)
  #devtools::install_github("rspatial/dismo")
  library(sp)
  library(rgdal)
  library(raster)
  library(dismo)
  library(adehabitatLT) # help(package='adehabitat') # help.search('angle',package='adehabitat')
  library(maps)       # for map.where
  library(mapdata)    # for worldHires
  library(sf)
  library(maptools)
  library(mgcv)
  library(ape)
  library(ncf)
  library(ncdf4)
  library(spdep)
  library(ROCR)
  library(gbm)
  library(tidyverse)
  library(viridis)
  #library(rJava)
  #library(ggmap)
  #library(RgoogleMaps)
  library(ggplot2)
  library(geosphere)
}))

# Open and rename file  ---------------------------------------------------

bwfw <- read.csv("~/Desktop/Hollingsinternship/Code/Hollings-Internship/bwfw_data.csv")
bwfw <- bwfw %>%
  rename("lat" = "latitude") %>%
  rename("lon" = "longitude") %>%
  rename("tag" = "DeploymentID") %>%
  rename("dtime" = "timestamp_gmt")
bwfw$dtime <- mdy_hms(bwfw$dtime)

# Land removal (from Hazen tutorial) --------------------------------------

# Function
removeland<-function(dataset,sp.obj){
  pt = SpatialPoints(matrix(c(dataset$lon,dataset$lat), nrow=length(dataset$lon)), proj4string=CRS("+proj=longlat +datum=WGS84"))
  place = over(pt, sp.obj)
  return(dataset[is.na(place),])
}

# Land map
map_background = maps::map('worldHires', fill=T, col='transparent')
back.IDs = sapply(strsplit(map_background$names, ":"), function(x) x[1])
back.sp = map2SpatialPolygons(map_background, IDs=back.IDs, proj4string=CRS("+proj=longlat +datum=WGS84"))

# Remove Data
bwfw_clean<-removeland(bwfw,back.sp)

# Plot data
data(wrld_simpl)
plot(wrld_simpl, xlim=c(-120,-80), ylim=c(0,60), axes=TRUE, col="light yellow")
coordinates(bwfw_clean) = ~lon+lat
points(bwfw_clean, pch=16, col=rgb(1, 0, 0, alpha=0.5), cex=0.3)

# Calculate timestep, distance and residual ---------------------------------------------------

bwfw_test <- bwfw_clean
unique_tags <- unique(bwfw_clean$tag) # List of unique tags

# Function to find timestep and distance between subsequent points of unique tag
find_timedist <- function(df) {
  timestep <- rep(NA, length(df$tag))
  distance <- rep(NA, length(df$tag))
  dist_error <- rep(NA, length(df$tag))
  for (i in seq_along(df$tag)) {
    if (i == 1) {
      timestep[[i]] <- NA
      distance[[i]] <- NA
    } else {
      timestep[[i]] <- difftime(df$dtime[[i]], df$dtime[[i-1]], units="hours") 
      distance[[i]] <- distm(c(df$lon[[i-1]], df$lat[[i-1]]), c(df$lon[[i]], df$lat[[i]]), fun = distHaversine)/1000
      # distance with error_radius of 2 data points subtracted
      dist_error[[i]] <- distance[[i]] - df$error_radius[[i]]/1000 - df$error_radius[[i-1]]/1000
    }
  }
  output <- cbind(df,timestep) %>%
    cbind(distance) %>%
    cbind(dist_error)
}

max_speed <- 100 # km/hour (twice the average speed of a blue whale cruising speed)

# Create the 5 new columns
apply_calcs <- function(df) {
  for (i in seq_along(unique_tags)) {
    if (i == 1) {
      bwfw_unique <- filter(df, tag == unique_tags[[i]])
      bwfw_unique <- find_timedist(bwfw_unique)
      bwfw_timedist <- bwfw_unique
    } else {
      bwfw_unique <- filter(df, tag == unique_tags[[i]])
      bwfw_unique <- find_timedist(bwfw_unique)
      bwfw_timedist <- rbind(bwfw_timedist,bwfw_unique)
    }
  }
  residual <- bwfw_timedist$distance - bwfw_timedist$timestep*max_speed
  residual_error <- bwfw_timedist$dist_error - bwfw_timedist$timestep*max_speed
  output <- cbind(bwfw_timedist,residual) %>%
    cbind(residual_error)
}

bwfw_calculations <- data.frame(apply_calcs(bwfw_test))
write.csv(bwfw_calculations, "~/Desktop/Hollingsinternship/Code/Hollings-Internship/bwfw_calculations.csv", row.names = FALSE)

# Find Anomalous Points ---------------------------------------------------

bwfw_anom1 <- filter(bwfw_calculations, residual > 0)
bwfw_anom2 <- filter(bwfw_calculations, residual_error > 0)

# Find unique tags with anomalous points
bwfw_anomtags1 <- unique(bwfw_anom1$tag)
bwfw_anomtags2 <- unique(bwfw_anom2$tag)

ggplot(bwfw_anom1,aes(x = tag)) +
  stat_count()

ggplot(bwfw_anom2,aes(x = tag)) +
  stat_count()

# Visualize Top 4 Anomalous Whales ----------------------------------------

weirdo_whales <- bwfw_anom2 %>%
  group_by(tag) %>%
  summarize(
    n = n()
  ) %>%
  filter(n>20)

weirdo_points <- semi_join(bwfw_anom2, weirdo_whales)
weirdo_whaletracks <- semi_join(bwfw_calculations,weirdo_whales)
weirdo_whaletracks$anom <- 0

# Data points for weirdo whales
weirdo_whaletracks <- full_join(bwfw_calculations,weirdo_whales)

# Setting anomalous points to 1
for (i in seq_along(weirdo_points$feature_id)) {
  index <- which(weirdo_whaletracks$feature_id == weirdo_points$feature_id[[i]])
  weirdo_whaletracks$anom[[index]] <- 1
}

non_anomalous <- filter(weirdo_whaletracks,"anom"==0)
anomalous <- filter(weirdo_whaletracks,"anom"==1)
# Plot anomalous tracks with anomalous points in different color
plot(wrld_simpl, xlim=c(-120,-80), ylim=c(0,60), axes=TRUE, col="light yellow")
coordinates(weirdo_whaletracks) = ~lon+lat
points(filter(weirdo_whaletracks,"anom"==0), col='red', pch=16, cex=0.5)

######







find_anomalies2 <- function(df,speed) {
  for (i in seq_along(unique_tags)) {
    if (i == 1) {
      bwfw_unique <- filter(df, tag == unique_tags[[i]])
      bwfw_unique <- find_timedist(bwfw_unique)
      bwfw_timedist <- bwfw_unique
    } else {
      bwfw_unique <- filter(df, tag == unique_tags[[i]])
      bwfw_unique <- find_timedist(bwfw_unique)
      bwfw_timedist <- rbind(bwfw_timedist,bwfw_unique)
    }
  }
  residual <- bwfw_timedist$dist_error - bwfw_timedist$timestep*max_speed
  bwfw_residual <- cbind(bwfw_timedist,residual)
  
# Find anomalous points
  anomalies <- filter(bwfw_residual, residual > 0)
  anomaly_points <- semi_join(bwfw_residual,anomalies)
  output <- anomaly_points
}
      
bwfw_anomalous1 <- find_anomalies1(bwfw_test,max_speed)
bwfw_anomalous2 <- find_anomalies2(bwfw_test,max_speed)
unique_tagsanom <- unique(bwfw_anomalous$tag)



weirdo_whale <- filter(bwfw_clean, DeploymentID == "2014CA-SPOT5-05784")
weird_whales <-
  



