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
  library(anytime)
  library(RColorBrewer)
  library(rasterVis)
}))