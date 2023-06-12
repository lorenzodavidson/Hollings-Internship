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

# Calculate timestep, distance and residual ---------------------------------------------------

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

max_speed <- 100 # km/hour (much more than max speed)

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
  # Accounting for gps error --> this is used for spotting anomalous points
  residual_error <- bwfw_timedist$dist_error - bwfw_timedist$timestep*max_speed 
  output <- cbind(bwfw_timedist,residual) %>%
    cbind(residual_error)
}

bwfw_calculations50 <- data.frame(apply_calcs(bwfw_clean))
write.csv(bwfw_calculations, "~/Desktop/Hollingsinternship/Code/Hollings-Internship/bwfw_calculations.csv", row.names = FALSE)

# Find and Visualize Anomalous Points ---------------------------------------------------

bwfw_anom <- filter(bwfw_calculations, residual_error > 0)
bwfw_anomlarge <- filter(bwfw_calculations, residual_error > 100)
bwfw_anomsmall <- filter(bwfw_anom, residual_error < 100)

# Find unique tags with anomalous points
bwfw_anomtags <- unique(bwfw_anom$tag) # 101 have anomalies
bwfw_anomtagslarge <- unique(bwfw_anomlarge$tag) # 43 (~half has 1 large anom)
bwfw_anomtagssmall <- unique(bwfw_anomsmall$tag) # 99 (all but 2 have small timestep anom)

ggplot(bwfw_anom,aes(x = tag)) +
  stat_count()
ggplot(bwfw_anomlarge,aes(x = tag)) +
  stat_count()
ggplot(bwfw_anomsmall,aes(x = tag)) +
  stat_count()

bwfw_calculations50$anom <- 0
# Setting large anomaly points to 1 and small anomaly points to 0.5
for (i in seq_along(bwfw_anomlarge$feature_id)) {
  index <- which(bwfw_calculations$feature_id == bwfw_anomlarge$feature_id[[i]])
  bwfw_calculations$anom[[index]] <- 1
}
for (i in seq_along(bwfw_anomsmall$feature_id)) {
  index <- which(bwfw_calculations$feature_id == bwfw_anomsmall$feature_id[[i]])
  bwfw_calculations$anom[[index]] <- 0.5
}

non_anomalous <- filter(bwfw_calculations,anom==0)
small_anomalous <- filter(bwfw_calculations,anom==0.5)
large_anomalous <- filter(bwfw_calculations,anom==1)
coordinates(non_anomalous) = ~lon+lat
coordinates(small_anomalous) = ~lon+lat
coordinates(large_anomalous) = ~lon+lat

data(wrld_simpl)
plot(wrld_simpl, xlim=c(-120,-80), ylim=c(0,60), axes=TRUE, col="light yellow")
points(non_anomalous, pch=16, col=rgb(0, 0, 1, alpha=0.5), cex=0.3)
points(small_anomalous, pch=16, col=rgb(0, 1, 0, alpha=0.5), cex=0.3)
points(large_anomalous, pch=16, col=rgb(1, 0, 0, alpha=0.5), cex=0.3)


# Remove 1st Anomaly of Pairs and Run Detention Again ---------------------

for (i in seq_along(bwfw_calculations50$feature_id)) {
  x <- i
  y < i-1
  if (bwfw_calculations50$anom[[x]] == 1 & bwfw_calculations50$anom[[y]] == 1) {
    bwfw_calculations50[-y]
  } else if (bwfw_calculations50$anom[[x]] == 1 & is.na(bwfw_calculations50$anom[[y]]) == TRUE) {
    bwfw_calculations50[-y]
  } else {
    bwfw_calculations50[[i]] <- bwfw_calculations50[[i]]
  }
}





# Visualize Top 4 Anomalous Whales ----------------------------------------

weirdo_whales <- bwfw_anom2 %>%
  group_by(tag) %>%
  summarize(
    n = n()
  ) %>%
  filter(n>20)

weirdo_points <- semi_join(bwfw_anom2, weirdo_whales)
weirdo_whaletracks <- semi_join(bwfw_calculations,weirdo_whales) # Data points for weirdo whales
weirdo_whaletracks$anom <- 0

# Setting anomalous points to 1
for (i in seq_along(weirdo_points$feature_id)) {
  index <- which(weirdo_whaletracks$feature_id == weirdo_points$feature_id[[i]])
  weirdo_whaletracks$anom[[index]] <- 1
}

non_anomalous <- filter(weirdo_whaletracks,anom==0)
anomalous <- filter(weirdo_whaletracks,anom==1)
coordinates(non_anomalous) = ~lon+lat
coordinates(anomalous) = ~lon+lat

data(wrld_simpl)
plot(wrld_simpl, xlim=c(-120,-80), ylim=c(0,60), axes=TRUE, col="light yellow")
points(non_anomalous, pch=16, col=rgb(0, 0, 1, alpha=0.5), cex=0.3)
points(anomalous, pch=16, col=rgb(1, 0, 0, alpha=0.5), cex=0.3)



