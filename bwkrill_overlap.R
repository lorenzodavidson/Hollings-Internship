# Libraries ---------------------------------------------------------------
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
}))

# Creating monthly blwh, krill, and bwkr files --------------------------------

# Open whale grd files for months of interest into xyz and raster format
fpath <- "~/Dropbox/blwh_sst_monthly/blwh/"
dfnames <- list.files(fpath)
dfnames <- list.files(fpath)[stringr::str_detect(dfnames,".grd")]

months <- c("03","04","05","06","07")
keywords1 <- c("_03_","_04_","_05_","_06_","_07_")
years <- c(1990:2020)
keywords2 <- as.character(years)

dfnames <- dfnames[Reduce("|", lapply(keywords1, function(x) grepl(x, dfnames)))] # Filter for years and months
dfnames <- dfnames[Reduce("|", lapply(keywords2, function(x) grepl(x, dfnames)))]

blwh_monthly_xyz <- list()
blwh_monthly_raster <- list()
for (i in dfnames) {
  raster <- raster(paste0(fpath,"/",i))
  blwh_monthly_raster[[i]] <- raster
  xyz <- data.frame(rasterToPoints(blwh_monthly_raster[[i]])) %>%
    rename("blwh" = "layer")
  blwh_monthly_xyz[[i]] <- xyz
}

# Open krill files and turn into xyz and raster format
totalkrill <- nc_open('TotalKril_CPUE.nc')
lon <- round(ncvar_get(totalkrill, "Longitude"),digits = 2)
lat <- round(ncvar_get(totalkrill, "Latitude"),digits = 2)
time <- anytime(ncvar_get(totalkrill, "Time"))
TotalKrill_CPUE <- ncvar_get(totalkrill, "TotalKrill_CPUE")

krill_monthly_xyz <- list()
krill_monthly_raster <- list()
for (t in seq_along(time)) {
  for (i in seq_along(lon)) {
    for (j in seq_along(lat)) {
      if (i == 1 & j == 1) {
        krill_rasterprep <- c(lon[[i]],lat[[j]],TotalKrill_CPUE[[i,j,t]])
      } else {
        new_row <- c(lon[[i]],lat[[j]],TotalKrill_CPUE[[i,j,t]])
        krill_rasterprep <- rbind(krill_rasterprep, new_row)
      }
    }
  }
  krill_monthly_xyz[[t]] <- data.frame(na.omit(krill_rasterprep)) %>%
    rename("x" = "X1") %>%
    rename("y" = "X2") %>%
    rename("krill" = "X3")
  #krill_raster <- rasterFromXYZ(krill_rasterprep,crs=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  #krill_monthly_raster[[t]] <- krill_raster
}

# Creating bwkr dataset
bwkr_monthly <- list()
for (t in seq_along(time)) {
  bwkr_monthly[[t]] <- inner_join(blwh_monthly_xyz[[t]],krill_monthly_xyz[[t]]) %>%
    rename("blwh" = "layer")
}

# Save all 5 datasets as csv or raster files
fpath_bwkr <- "~/Dropbox/blwh_sst_monthly/blwh_krill/bwkr_monthly"
fpath_blwh <- "~/Dropbox/blwh_sst_monthly/blwh_krill/blwh_monthly"
fpath_krill <- "~/Dropbox/blwh_sst_monthly/blwh_krill/krill_monthly"
for (i in seq_along(years)) {
  for (j in seq_along(months)) {
    x <- i*5 - 5 + j
    bwkr_name <- paste0("bwkr_", years[[i]], "_", months[[j]],".csv") # Need to figure out
    write.csv(bwkr_monthly[[x]],paste0(fpath_bwkr,"/",bwkr_name),row.names = FALSE)
    
    blwh_name <- paste0("blwh_", years[[i]], "_", months[[j]])
    write.csv(blwh_monthly_xyz[[x]],paste0(fpath_blwh,"/",blwh_name,".csv"),row.names = FALSE)
    writeRaster(blwh_monthly_raster[[x]],paste0(fpath_blwh,"/",blwh_name,".grd"),format = "raster") # Figure out how to save rasters
    
    krill_name <- paste0("krill_", years[[i]], "_", months[[j]])
    write.csv(krill_monthly_xyz[[x]],paste0(fpath_krill,"/",krill_name,".csv"),row.names = FALSE)
    writeRaster(krill_monthly_raster[[x]],paste0(fpath_krill,"/",krill_name,".grd"),format = "raster")
  }
}





# Trying out with 1 df at a time
blwh_t1 <- data.frame(rasterToPoints(blwh_monthly_long[[1]]))
krill_t1 <- data.frame(krill_monthly_long[[1]])
krill_t1na <- data.frame(na.omit(krill_monthly_long[[1]])) %>%
  rename("x" = "X1") %>%
  rename("y" = "X2") %>%
  rename("krill" = "X3")

blkr_t1test <- inner_join(blwh_t1,krill_t1na)


blwh_t1 <- blwh_monthly_xyz[[1]]
krill_t1 <- krill_monthly_xyz[[1]]

write.csv(krill_t1na, "~/Desktop/Hollingsinternship/Code/Hollings-Internship/krill_t1.csv", row.names = FALSE)
write.csv(blwh_t1, "~/Desktop/Hollingsinternship/Code/Hollings-Internship/blwh_t1.csv", row.names = FALSE)

krill_t1$x %>%
  round(digits = 2)

round(krill_t1$y, digits = 2)

write.csv(krill_t1, "~/Desktop/Hollingsinternship/Code/Hollings-Internship/krill_t1.csv", row.names = FALSE)

filter(krill_t1,y==39.95)

for (i in seq_along(lon)) {
  for (j in seq_along(lat)) {
    if (i == 1 & j == 1) {
      krill_rasterprep <- c(lon[[i]],lat[[j]],TotalKrill_CPUE[[i,j,1]])
    } else {
      new_row <- c(lon[[i]],lat[[j]],TotalKrill_CPUE[[i,j,1]])
      krill_rasterprep <- rbind(krill_rasterprep, new_row)
    }
  }
}
krill_t1na <- data.frame(na.omit(krill_rasterprep))
  
krill_raster <- rasterFromXYZ(krill_rasterprep)
krill_t1test <- data.frame(rasterToPoints(krill_raster))

head(filter(krill_t1test,y==39.95))











blwh_t1 <- data.frame(rasterToPoints(blwh_monthly[[1]]))
krill_t1 <- data.frame(rasterToPoints(krill_monthly[[1]]))



krill_wtf <- krill_t1test %>%
  filter(x == -127.05 | y == 47.95)

blwh_t1[[1,1]]

# Opening whale and krill data --------------------------------------------
fpath_blwh <- "~/Dropbox/blwh_sst_monthly/blwh_krill/blwh_monthly"
fpath_krill <- "~/Dropbox/blwh_sst_monthly/blwh_krill/krill_monthly"
fpath_bwkr <- "~/Dropbox/blwh_sst_monthly/blwh_krill/bwkr_monthly"

dfnames_bwkr <- list.files(fpath_bwkr)
dfnames_blwh <- list.files(fpath_blwh)
dfnames_krill <- list.files(fpath_krill)

blwh_monthly_raster <- list()
for(i in dfnames_blwh[stringr::str_detect(dfnames_blwh, ".grd")]){
  tt <- raster(paste0(fpath_blwh,"/",i))
  blwh_monthly_raster[[i]] <- tt
}

krill_monthly_raster <- list()
for(i in dfnames_krill[stringr::str_detect(dfnames_krill, ".grd")]){
  tt <- raster(paste0(fpath_krill,"/",i))
  krill_monthly_raster[[i]] <- tt
}

blwh_monthly_xyz <- list()
for(i in dfnames_blwh[stringr::str_detect(dfnames_blwh, ".csv")]){
  tt <- read.csv(paste0(fpath_blwh,"/",i))
  blwh_monthly_xyz[[i]] <- tt
}

krill_monthly_xyz <- list()
for(i in dfnames_krill[stringr::str_detect(dfnames_krill, ".csv")]){
  tt <- read.csv(paste0(fpath_krill,"/",i))
  krill_monthly_xyz[[i]] <- tt
}

bwkr_monthly <- list()
for(i in dfnames_bwkr[stringr::str_detect(dfnames_bwkr, ".csv")]){
  tt <- read.csv(paste0(fpath_bwkr,"/",i))
  bwkr_monthly[[i]] <- tt
}


# Overlap Metrics Code ----------------------------------------------------







writeRaster(krill_monthly, "~/Desktop/Hollingsinternship/Code/Hollings-Internship/krill_monthly.nc", format="CDF")

write.csv(krill_monthly_long, "~/Desktop/Hollingsinternship/Code/Hollings-Internship/krill_monthly_long.csv", row.names = FALSE)


krill_t1_rasterprep1 <- c(lon[[1]],lat[[1]],krill_t1[[1,1]])
new_row <- c(lon[[1]],lat[[2]],krill_t1[[1,2]])
krill_t1_rasterprep1 <- rbind(krill_t1_rasterprep1, new_row)
krill_t1_raster <- rasterFromXYZ(krill_t1_rasterprep1,crs=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

plot(krill_t1_raster)

krill_monthly<- list()
for (i in seq_along(time)) {
  x1 <- TotalKrill_CPUE[,,i]
  x2 <- t(x1)
  krill_monthly[[i]] <- raster(x2,
                               xmn=range(lon[1]), xmx=range(lon[2])
                               ymn=range(lat[1]), ymx=range(lat[2]),
                               crs=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
                               )
}

dat1=list()
dat1$x=lon
dat1$y=lat
dat1$z=t(tmp.array.day)

r <-raster(
  dat1$z,
  xmn=range(dat1$x)[1], xmx=range(dat1$x)[2],
  ymn=range(dat1$y)[1], ymx=range(dat1$y)[2], 
  crs=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
)


krill_t1 <- TotalKrill_CPUE[,,1]
krill_t1_transpose <- t(krill_t1)
plot(krill_monthly[[2]])

# Flipping Krill, Lat and Lon Arrays for Visualization
krill_flipped <- array(dim = c(185,180,155))
for (t in 1:155) {
  for (i in 1:185) {
    for (j in 1:180) {
      krill_flipped[[i,181-j,t]] <- krill_monthly[[i,j,t]]
    }
  }
}

lon_flipped <- array(NA,185)
for (i in 1:185) {
  lon_flipped[[186-i]] <- -lon[[i]]
}

lat_flipped <- array(NA,180)
for (j in 1:180) {
  lat_flipped[[181-j]] <- lat[[j]]
}

# Rasterizing Krill Data
krill_t1 <- krill_monthly[,,1]
krill_t1_flipped <- krill_flipped[,,1]

rotate <- function(x) t(apply(x, 2, rev))

krill_t1_transpose <- t(krill_t1)

plot(raster(krill_t1))
plot(raster(krill_t1_flipped))
plot(raster(krill_t1_rotated))
plot(raster(krill_t1_transpose))



krill_t1_raster <- data_frame()
raster()
for (i in seq_along(lon)) {
  for (j in seq_along(lat)) {
    
  }
}



for (i in seq_along(t)) {
  
}




plot(krill_flipped[,,6])

image(lon_flipped,lat_flipped,krill_flipped[,,6])






lon_cal <- matrix(NA,185,180)
lat_cal <- matrix(NA,185,180)
for (i in 1:185) {
  for (j in 1:180) {
    lon_cal[[i,j]] <- lon[[i]]
    lat_cal[[i,j]] <- lat[[j]]
  }
}

write.csv(krill_t1, "~/Desktop/Hollingsinternship/Code/Hollings-Internship/krill_t1.csv", row.names = FALSE)
write.csv(lon_cal, "~/Desktop/Hollingsinternship/Code/Hollings-Internship/lon_cal.csv", row.names = FALSE)
write.csv(lat_cal, "~/Desktop/Hollingsinternship/Code/Hollings-Internship/lat_cal.csv", row.names = FALSE)


himage(lon,lat,krill_t1)

plot(krill_t1)

data(wrld_simpl)
plot(wrld_simpl, xlim=c(-120,-80), ylim=c(0,60), axes=TRUE, col="light yellow")
points(non_anomalous, pch=16, col=rgb(0, 0, 1, alpha=0.5), cex=0.3)
points(anomalous, pch=16, col=rgb(1, 0, 0, alpha=0.5), cex=0.3)