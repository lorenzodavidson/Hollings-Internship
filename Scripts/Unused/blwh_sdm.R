# Libraries ---------------------------------------------------------------
library(aniMotum)
require(tidyverse)
require(patchwork)
require(sf)

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
# Clean up data for Animotum ----------------------------------------------

bwfw <- read.csv("~/Desktop/Hollingsinternship/Code/Hollings-Internship/bwfw_data.csv")
bwfw <- bwfw %>%
  rename("lat" = "latitude", "lon" = "longitude", "date" = "timestamp_gmt", "id" = "DeploymentID",
        "lc" = "loc_class", "smaj" = "semi_major", "smin" = "semi_minor","eor" = "ellipse")
bwfw$date <- mdy_hms(bwfw$date)
bwfw <- subset(bwfw, select = c("id","lon","lat","date","lc","smaj","smin","eor"))

library("readxl")
bwfw_meta <- data.frame(read_excel("~/Desktop/Hollingsinternship/Code/Hollings-Internship/bwfw_metadata.xlsx"))
bwfw_meta <- bwfw_meta[,11:17]
bwfw_meta <- bwfw_meta[,-2:-6] %>%
  rename("id" = "DeploymentID") %>%
  filter(SpeciesCommonName == "Blue whale") 

blwh_prep <- format_data(semi_join(bwfw,bwfw_meta))
