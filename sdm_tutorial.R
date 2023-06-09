#Load libraries 
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
}))


# Data Preparation --------------------------------------------------------
library(dismo)

file <- paste0(system.file(package="dismo"), "/ex/bradypus.csv")
bradypus <- read.table(file,header=TRUE,sep=",")

bwfw <- read.csv("~/Desktop/Hollingsinternship/Code/myrepo/bwfw_data.csv")
bwfw <- bwfw %>%
  rename("lat" = "latitude") %>%
  rename("lon" = "longitude")

# Land removal (from Hazen tutorial) --------------------------------------

# Function + Remove
removeland<-function(dataset,sp.obj){
  pt = SpatialPoints(matrix(c(dataset$lon,dataset$lat), nrow=length(dataset$lon)), proj4string=CRS("+proj=longlat +datum=WGS84"))
  place = over(pt, sp.obj)
  return(dataset[is.na(place),])
}

# Land map
map_background = maps::map('worldHires', fill=T, col='transparent')
back.IDs = sapply(strsplit(map_background$names, ":"), function(x) x[1])
back.sp = map2SpatialPolygons(map_background, IDs=back.IDs, proj4string=CRS("+proj=longlat +datum=WGS84"))

# Clean data
bwfw_clean<-removeland(bwfw,back.sp)
# bwfw_filtered <- bwfw_clean %>%
  # filter(0 < lat | lat < 60 | -140 < lon | lon < -75)

# Plot
library(maptools)
data(wrld_simpl)
plot(wrld_simpl, xlim=c(-120,-80), ylim=c(0,60), axes=TRUE, col="light yellow")
coordinates(bwfw_clean) = ~lon+lat
points(bwfw_clean, pch=16, col=rgb(1, 0, 0, alpha=0.5), cex=0.3)

# Sampling Bias -----------------------------------------------------------
library(sp)
coordinates(bwfw) <- ~lon+lat
crs(bwfw) <- crs(wrld_simpl)
class(bwfw)

r <- raster(bwfw) # create a RasterLayer with the extent of acgeo
res(r) <- 1 # set the resolution of the cells to (for example) 1 degree
r <- extend(r, extent(r)+1)
bwfw_sel <- gridSample(bwfw, r, n=1) # sample
p <- rasterToPolygons(r)
plot(p, border='gray')
points(bwfw)
points(bwfw_sel, cex=1, col='red', pch='x') # selected points in red

# Background/Absence Data -------------------------------------------------

# Unable to mask out land
files <- list.files(path=paste(system.file(package="dismo"), '/ex',
                               sep=''),  pattern='grd',  full.names=TRUE )
land_mask <- raster(files[1])
# set.seed(1963) # select 500 random points with same seed
bg <- randomPoints(!mask, 500)

par(mfrow=c(1,2))
plot(!is.na(mask), legend=FALSE)
points(bg, cex=0.5)

e <- extent(-140, -75, 0, 60)
bg2 <- randomPoints(mask, 50, ext=e)
plot(!is.na(mask), legend=FALSE)
plot(e, add=TRUE, col='black')
points(bg2, cex=0.5)

