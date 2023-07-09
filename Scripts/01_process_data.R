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

# Function to create raster from netCDF
create_ROMS_RASTER=function(nc,dname,month_nc,template){ #nc file, name of variable, month of interest, template raster
  nc.data=nc_open(nc)
  lat <- ncvar_get(nc.data,'Latitude')
  lon <- ncvar_get(nc.data,'Longitude')
  nrows <- length(lat); ncols <- length(lon)
  tt <- ncvar_get(nc.data, "Time")
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

# Total krill
totalkrill <- nc_open("~/Dropbox/krill_raw/TotalKril_CPUE.nc", write = TRUE)
lon <- ncvar_get(totalkrill, "Longitude"); lat <- ncvar_get(totalkrill, "Latitude", verbose = F)
tt <- ncvar_get(totalkrill, "Time")
dates <- as.Date(as.POSIXct(tt, tz = "GMT", origin = "1970-01-01"))
years <- year(dates); months <- as.numeric(month(dates))
monthabs <- month.abb[as.numeric(month(dates))]
dates <- format(dates,"%Y-%m")

nc <- "~/Dropbox/krill_raw/TotalKril_CPUE.nc"
dname = 'TotalKrill_CPUE'
var = "totalkrill"
temp <- raster(ncol=185, nrow=180, xmn=-134, xmx=-115.5, ymn=30, ymx=48)
krill_monthly_raster <- list()
krill_monthly_xyz <- list()
for (t in seq_along(dates)) {
  month_nc <- dates[[t]]
  krill_monthly_raster[[t]] <- create_ROMS_RASTER(nc,dname,month_nc,temp)
  
  krill_monthly_xyz[[t]] <- data.frame(na.omit(rasterToPoints(krill_monthly_raster[[t]]))) %>%
    rename("lon"="x","lat"="y","krill"="layer")
}

# EPAC
epackrill <- nc_open("~/Dropbox/krill_raw/EPAC_CPUE.nc", write = TRUE)
tt_epac <- ncvar_get(epackrill, "Time")
dates_epac <- as.Date(as.POSIXct(tt_epac, tz = "GMT", origin = "1970-01-01"))
years_epac <- year(dates_epac); months_epac <- as.numeric(month(dates_epac))
monthabs_epac <- month.abb[as.numeric(month(dates_epac))]
dates_epac <- format(dates_epac,"%Y-%m")

nc <- "~/Dropbox/krill_raw/EPAC_CPUE.nc"
dname = 'EPAC_CPUE'
temp <- raster(ncol=185, nrow=180, xmn=-134, xmx=-115.5, ymn=30, ymx=48)
epac_monthly_raster <- list()
epac_monthly_xyz <- list()
for (t in seq_along(dates_epac)) {
  month_nc <- dates_epac[[t]]
  epac_monthly_raster[[t]] <- create_ROMS_RASTER(nc,dname,month_nc,temp)
  
  epac_monthly_xyz[[t]] <- data.frame(na.omit(rasterToPoints(epac_monthly_raster[[t]]))) %>%
    rename("lon"="x","lat"="y","epac"="layer")
}

# TSPIN
tspinrill <- nc_open("~/Dropbox/krill_raw/TSPIN_CPUE.nc", write = TRUE)
tt_tspin <- ncvar_get(tspinrill, "Time")
dates_tspin <- as.Date(as.POSIXct(tt_tspin, tz = "GMT", origin = "1970-01-01"))
years_tspin <- year(dates_tspin); months_tspin <- as.numeric(month(dates_tspin))
monthabs_tspin <- month.abb[as.numeric(month(dates_tspin))]
dates_tspin <- format(dates_tspin,"%Y-%m")

nc <- "~/Dropbox/krill_raw/TSPIN_CPUE.nc"
dname = 'TSPIN_CPUE'
temp <- raster(ncol=185, nrow=180, xmn=-134, xmx=-115.5, ymn=30, ymx=48)
tspin_monthly_raster <- list()
tspin_monthly_xyz <- list()
for (t in seq_along(dates_tspin)) {
  month_nc <- dates_tspin[[t]]
  tspin_monthly_raster[[t]] <- create_ROMS_RASTER(nc,dname,month_nc,temp)
  
  tspin_monthly_xyz[[t]] <- data.frame(na.omit(rasterToPoints(tspin_monthly_raster[[t]]))) %>%
    rename("lon"="x","lat"="y","tspin"="layer")
}

# Creating bwkr_monthly_xyz and bwkr_all datasets --------------------------

bwkr_monthly_xyz <- list()
for (t in seq_along(dates)) {
  if (t <= 60) {
    bwkr_monthly_xyz[[t]] <- inner_join(blwh_monthly_xyz[[t]],krill_monthly_xyz[[t]])
  } else {
    x <- t-60
    bwkr_monthly_xyz[[t]] <- inner_join(blwh_monthly_xyz[[t]],krill_monthly_xyz[[t]])
    bwkr_monthly_xyz[[t]] <- left_join(bwkr_monthly_xyz[[t]],epac_monthly_xyz[[x]])
    bwkr_monthly_xyz[[t]] <- left_join(bwkr_monthly_xyz[[t]],tspin_monthly_xyz[[x]])
  }
}

# Includes epac and tspin
bwkr_all_xyz <- data.frame(bwkr_monthly_xyz[[1]])
bwkr_all_xyz$epac <- NA
bwkr_all_xyz$tspin <- NA
bwkr_all_xyz$date <- dates[[1]]
bwkr_all_xyz$year <- years[[1]]
bwkr_all_xyz$month <- months[[1]]
bwkr_all_xyz$monthabs <- monthabs[[1]]
for (i in 2:155) {
  if (i <= 60) {
    bwkr_add <- bwkr_monthly_xyz[[i]]
    bwkr_add$epac <- NA
    bwkr_add$tspin <- NA
    bwkr_add$date <- dates[[i]]
    bwkr_add$year <- years[[i]]
    bwkr_add$month <- months[[i]]
    bwkr_add$monthabs <- monthabs[[i]]
    bwkr_all_xyz <- rbind(bwkr_all_xyz, bwkr_add)
  } else {
    bwkr_add <- bwkr_monthly_xyz[[i]]
    bwkr_add$date <- dates[[i]]
    bwkr_add$year <- years[[i]]
    bwkr_add$month <- months[[i]]
    bwkr_add$monthabs <- monthabs[[i]]
    bwkr_all_xyz <- rbind(bwkr_all_xyz, bwkr_add)
  }
}

# Creating final blwh and krill raster datasets  --------------------------
ext <- extent(-127.5, -115.5, 30, 48) # Extent of bwkr raster with edges

blwh_all_raster <- stack() # Stack of 155 rasters containing only bwkr shared datapoints
krill_all_raster <- stack()
epac_all_raster <- stack() # Stack of 95 rasters
tspin_all_raster <- stack()
for (t in seq_along(dates)) {
  bwkr_rasterprep <- bwkr_monthly_xyz[[t]][,-7:-10]
  if (t ==1) {
    blwh_all_raster <- extend(rasterFromXYZ(bwkr_rasterprep[c(-4:-6)]),ext)
    krill_all_raster <- extend(rasterFromXYZ(bwkr_rasterprep[c(-3,-5,-6)]),ext)
  } else {
    if (t <= 60) {
      blwh_raster <- extend(rasterFromXYZ(bwkr_rasterprep[,-4]),ext)
      krill_raster <- extend(rasterFromXYZ(bwkr_rasterprep[,-3]),ext)
      blwh_all_raster <- stack(blwh_all_raster, blwh_raster)
      krill_all_raster <- stack(krill_all_raster, krill_raster) 
    } else {
      if (t == 61) {
        blwh_raster <- extend(rasterFromXYZ(bwkr_rasterprep[c(-4:-6)]),ext)
        krill_raster <- extend(rasterFromXYZ(bwkr_rasterprep[c(-3,-5:-6)]),ext)
        blwh_all_raster <- stack(blwh_all_raster, blwh_raster)
        krill_all_raster <- stack(krill_all_raster, krill_raster)
        
        epac_all_raster <- extend(rasterFromXYZ(bwkr_rasterprep[c(-3:-4,-6)]),ext)
        tspin_all_raster <- extend(rasterFromXYZ(bwkr_rasterprep[c(-3:-5)]),ext)
      } else {
        blwh_raster <- extend(rasterFromXYZ(bwkr_rasterprep[c(-4:-6)]),ext)
        krill_raster <- extend(rasterFromXYZ(bwkr_rasterprep[c(-3,-5:-6)]),ext)
        epac_raster <- extend(rasterFromXYZ(bwkr_rasterprep[c(-3:-4,-6)]),ext)
        tspin_raster <- extend(rasterFromXYZ(bwkr_rasterprep[c(-3:-5)]),ext)
        
        blwh_all_raster <- stack(blwh_all_raster, blwh_raster)
        krill_all_raster <- stack(krill_all_raster, krill_raster)
        epac_all_raster <- stack(epac_all_raster, epac_raster)
        tspin_all_raster <- stack(tspin_all_raster, tspin_raster)
      }
    }
  }
}

# Save bwkr_all, blwh_all_raster, krill_all_raster, and monthly  --------

processed_path <- "~/Desktop/Hollingsinternship/Hollings-Internship/Processed/"

write.csv(dates,paste0(processed_path,"dates.csv"),row.names = FALSE)
write.csv(bwkr_all_xyz,paste0(processed_path,"bwkr_all.csv"),row.names = FALSE) # bwkr_all
writeRaster(blwh_all_raster,paste0(processed_path,"blwh_all/","blwh.grd"),
            format = "raster",bylayer=TRUE,suffix=dates)
writeRaster(krill_all_raster,paste0(processed_path,"krill_all/","krill.grd"),
            format = "raster",bylayer=TRUE,suffix=dates)
writeRaster(epac_all_raster,paste0(processed_path,"epac_all/","epac.grd"),
            format = "raster",bylayer=TRUE,suffix=dates_epac)
writeRaster(tspin_all_raster,paste0(processed_path,"tspin_all/","tspin.grd"),
            format = "raster",bylayer=TRUE,suffix=dates_tspin)
