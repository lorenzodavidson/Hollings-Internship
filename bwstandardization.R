library(ggplot2)
library(maps)
install.packages("ggmap")
library(ggmap)
library(tidyverse)

# Open and tidy file
bwfw <- read.csv("~/Desktop/Hollingsinternship/Code/myrepo/bwfw_data.csv")

bwfw_date_time <- separate(bwfw,timestamp_gmt, into = c("date", "time"), sep = " ")

# Frequency plot
bwfw_date_time %>%
  ggplot(aes(date)) + 
  stat_count()

# Standardize to 1 data per day per tagged whale
bwfw_standardized <- group_by(bwfw_date_time,DeploymentID,date) %>%
  summarise(mean_lat = mean(latitude,na.rm = TRUE), mean_lon = mean(longitude,na.rm = TRUE))

bwfw_standardized %>%
  ggplot(aes(date)) + 
  stat_count()

# Map
map <- get_stamenmap(
  bbox = c(left=-135,bottom=0,right=-75,top=50),
  maptype = "terrain",
  zoom=4)
map <- ggmap(map, legend = "topright")

# Plotting telemetry data
map +
  geom_point(bwfw_standardized,mapping=aes(x=mean_lon,y=mean_lat,color=DeploymentID),
             show.legend = FALSE)

map +
  stat_density2d(bwfw_standardized,
                 mapping=aes(x = mean_lon, y = mean_lat),
                 linewidth = 2, bins = 100, alpha = 0.3, na.rm = TRUE,
                 geom = "polygon"
  )