#Libraries ---------------------------------------------------------------
suppressPackageStartupMessages(suppressWarnings({
  library(raster)
  library(anytime)
  library(ncf)
  library(ncdf4)
  library(tidyverse)
  library(rgdal)
  library(mondate)
}))

# Open whale grd files for months of interest into xyz and raster --------

fpath <- "~/Dropbox/blwhsst_HS_monthly_raw/blwh/"
dfnames <- list.files(fpath)
dfnames <- list.files(fpath)[stringr::str_detect(dfnames,".grd")]

keywords1 <- c("_04_","_05_","_06_","_07_","_08_") # filtering for months
years <- c(1990:2020)
keywords2 <- as.character(years)

# Filter for years and months (March-July of 1990-2020)
dfnames <- dfnames[Reduce("|", lapply(keywords1, function(x) grepl(x, dfnames)))] 
dfnames <- dfnames[Reduce("|", lapply(keywords2, function(x) grepl(x, dfnames)))]

blwh_monthly_xyz <- list()
blwh_monthly_raster <- list()
for (i in dfnames) {
  raster <- raster(paste0(fpath,"/",i))
  blwh_monthly_raster[[i]] <- raster
  xyz <- data.frame(rasterToPoints(blwh_monthly_raster[[i]])) %>%
    rename("lon" = "x", "lat" = "y", "blwh" = "layer")
  blwh_monthly_xyz[[i]] <- xyz
}

# Open krill netCDF file into xyz and raster format --------

totalkrill <- nc_open("~/Dropbox/krill_raw/TotalKril_CPUE.nc", write = TRUE)
lon <- ncvar_get(totalkrill, "Longitude")
lat <- ncvar_get(totalkrill, "Latitude", verbose = F)
tt <- ncvar_get(totalkrill, "Time")
dates <- as.Date(as.POSIXct(tt, tz = "GMT", origin = "1970-01-01"))
years <- year(dates)
months <- month.abb[as.numeric(month(dates))]
monthnums <- as.numeric(month(dates))
dates <- format(dates,"%Y-%m")
TotalKrill_CPUE <- ncvar_get(totalkrill, "TotalKrill_CPUE")

# Function to create raster from netCDF
create_ROMS_RASTER=function(nc,dname,month_nc,template){ #nc file, name of variable, month of interest, template raster
  nc.data=nc_open(nc)
  lat <- ncvar_get(nc.data,'Latitude')
  lon <- ncvar_get(nc.data,'Longitude')
  nrows <- length(lat); ncols <- length(lon)
  tt <- ncvar_get(totalkrill, "Time")
  dates <- format(as.Date(as.POSIXct(tt, tz = "GMT", origin = "1970-01-01")),"%Y-%m")
  tmp.array <- ncvar_get(nc.data, dname)
  fillvalue <- ncatt_get(nc.data, dname, "_FillValue")
  names=unlist(lapply(dates,function(x)gsub(" UTC","",x)))
  index=grep(month_nc,names)
  
  # Create rasters
  tmp.array.day=tmp.array[,,index]
  tmp.array.day[tmp.array.day==fillvalue$value]=NA #setting fill value
  
  dat1=list()
  dat1$x=lon
  dat1$y=lat
  dat1$z=t(tmp.array.day)
  
  r <-raster(
    dat1$z,
    xmn=range(dat1$x)[1], xmx=range(dat1$x)[2],
    ymn=range(dat1$y)[1], ymx=range(dat1$y)[2], 
  )
  r <- raster::resample(r, template, method="bilinear")
}

# krill monthly rasters and xyz datasets
nc <- "~/Dropbox/krill_raw/TotalKril_CPUE.nc"
dname = 'TotalKrill_CPUE'
temp <- raster(ncol=185, nrow=180, xmn=-134, xmx=-115.5, ymn=30, ymx=48)

krill_monthly_raster <- list()
krill_monthly_xyz <- list()
for (t in seq_along(dates)) {
  month_nc <- dates[[t]]
  krill_monthly_raster[[t]] <- create_ROMS_RASTER(nc,dname,month_nc,temp)
  
  krill_monthly_xyz[[t]] <- data.frame(na.omit(rasterToPoints(krill_monthly_raster[[t]]))) %>%
    rename("lon"="x","lat"="y","krill"="layer")
}

# Creating bwkr_monthly_xyz and bwkr_all datasets --------------------------

bwkr_monthly_xyz <- list()
for (t in seq_along(dates)) {
  bwkr_monthly_xyz[[t]] <- inner_join(blwh_monthly_xyz[[t]],krill_monthly_xyz[[t]])
}

bwkr_all_xyz <- data.frame(bwkr_monthly_xyz[[1]])
bwkr_all_xyz$date <- dates[[1]]
bwkr_all_xyz$year <- years[[1]]
bwkr_all_xyz$month <- months[[1]]
bwkr_all_xyz$monthnum <- monthnums[[1]]
for (i in 2:155) {
  bwkr_add <- bwkr_monthly_xyz[[i]]
  bwkr_add$date <- dates[[i]]
  bwkr_add$year <- years[[i]]
  bwkr_add$month <- months[[i]]
  bwkr_add$monthnum <- monthnums[[i]]
  bwkr_all_xyz <- rbind(bwkr_all_xyz, bwkr_add)
}

# Creating final blwh and krill raster datasets  --------------------------
ext <- extent(-127.5, -115.5, 30, 48) # Extent of bwkr raster with edges

blwh_all_raster <- stack() # Stack of 155 rasters containing only bwkr shared datapoints
krill_all_raster <- stack()
for (t in seq_along(dates)) {
  bwkr_rasterprep <- bwkr_monthly_xyz[[t]][,-5:-7]
  if (t ==1) {
    blwh_all_raster <- extend(rasterFromXYZ(bwkr_rasterprep[,-4]),ext)
    krill_all_raster <- extend(rasterFromXYZ(bwkr_rasterprep[,-3]),ext)
  } else {
    blwh_raster <- extend(rasterFromXYZ(bwkr_rasterprep[,-4]),ext)
    krill_raster <- extend(rasterFromXYZ(bwkr_rasterprep[,-3]),ext)
    blwh_all_raster <- stack(blwh_all_raster, blwh_raster)
    krill_all_raster <- stack(krill_all_raster, krill_raster)  
  }
}
blwh_all_raster <- setZ(blwh_all_raster, dates)
krill_all_raster <- setZ(krill_all_raster, dates)

# Monthly raster stacks (31 layers each)
blwh_apr_raster <- stack(); krill_apr_raster <- stack()
blwh_may_raster <- stack(); krill_may_raster <- stack()
blwh_jun_raster <- stack(); krill_jun_raster <- stack()
blwh_jul_raster <- stack(); krill_jul_raster <- stack()
blwh_aug_raster <- stack(); krill_aug_raster <- stack()

apr <- seq(1,151,by=5); may <- seq(2,152,by=5); jun <- seq(3,153,by=5)
jul <- seq(4,154,by=5); aug <- seq(5,155,by=5)

blwh_apr_raster <- subset(blwh_all_raster, apr)
krill_apr_raster <- subset(krill_all_raster, apr)
blwh_may_raster <- subset(blwh_all_raster, may)
krill_may_raster <- subset(krill_all_raster, may)
blwh_jun_raster <- subset(blwh_all_raster, jun)
krill_jun_raster <- subset(krill_all_raster, jun)
blwh_jul_raster <- subset(blwh_all_raster, jul)
krill_jul_raster <- subset(krill_all_raster, jul)
blwh_aug_raster <- subset(blwh_all_raster, aug)
krill_aug_raster <- subset(krill_all_raster, aug)


# Save bwkr_all, blwh_all_raster, krill_all_raster, and monthly  --------

processed_path <- "~/Desktop/Hollingsinternship/Hollings-Internship/Processed/"

write.csv(bwkr_all_xyz,paste0(processed_path,"bwkr_all.csv"),row.names = FALSE) # bwkr_all
writeRaster(blwh_all_raster,paste0(processed_path,"blwh_all_raster.grd"),format = "raster")
writeRaster(krill_all_raster,paste0(processed_path,"krill_all_raster.grd"),format = "raster")

# Monthly
writeRaster(blwh_apr_raster,paste0(processed_path,"blwh_apr_raster.grd"),format = "raster")
writeRaster(krill_apr_raster,paste0(processed_path,"krill_apr_raster.grd"),format = "raster")
writeRaster(blwh_may_raster,paste0(processed_path,"blwh_may_raster.grd"),format = "raster")
writeRaster(krill_may_raster,paste0(processed_path,"krill_may_raster.grd"),format = "raster")
writeRaster(blwh_jun_raster,paste0(processed_path,"blwh_jun_raster.grd"),format = "raster")
writeRaster(krill_jun_raster,paste0(processed_path,"krill_jun_raster.grd"),format = "raster")
writeRaster(blwh_jul_raster,paste0(processed_path,"blwh_jul_raster.grd"),format = "raster")
writeRaster(krill_jul_raster,paste0(processed_path,"krill_jul_raster.grd"),format = "raster")
writeRaster(blwh_aug_raster,paste0(processed_path,"blwh_aug_raster.grd"),format = "raster")
writeRaster(krill_aug_raster,paste0(processed_path,"krill_aug_raster.grd"),format = "raster")
