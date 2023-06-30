# Libraries ---------------------------------------------------------------

suppressPackageStartupMessages(suppressWarnings({
  #library(devtools)
  #devtools::install_github("rspatial/dismo")
  library(sp)
  library(rgdal)
  library(raster)
  # library(dismo)
  # library(adehabitatLT) # help(package='adehabitat') # help.search('angle',package='adehabitat')
  # library(maps)       # for map.where
  # library(mapdata)    # for worldHires
  library(sf)
  # library(maptools)
  # library(mgcv)
  # library(ape)
  # library(ncf)
  # library(ncdf4)
  # library(spdep)
  # library(ROCR)
  # library(gbm)
  library(tidyverse)
  # library(viridis)
  #library(rJava)
  #library(ggmap)
  #library(RgoogleMaps)
  library(ggplot2)
  library(geosphere)
  library(anytime)
  library(RColorBrewer)
  library(rasterVis)
  library(terra)
}))


# Open files --------------------------------------------------------------

processed_path <- "~/Desktop/Hollingsinternship/Hollings-Internship/Processed/"
dates <- t(read.csv(paste0(processed_path,"/dates.csv")))

blwh_ROMS <- readRDS("~/Dropbox/blwh_telemetry_raw/blwh_ROMSvar_list.rds")

krill_path <- paste0(processed_path,"krill_all/")
krill_names <- list.files(krill_path)
krill_names <- list.files(krill_path)[stringr::str_detect(krill_names,".grd")]

krill_all_raster <- stack()
for (i in krill_names) {
  raster <- raster(paste0(krill_path,"/",i))
  krill_all_raster <- stack(krill_all_raster,raster)
}
krill_all_raster <- setZ(krill_all_raster, dates)
names(krill_all_raster) <- dates

blwh_path <- paste0(processed_path,"blwh_all/")
blwh_names <- list.files(blwh_path)
blwh_names <- list.files(blwh_path)[stringr::str_detect(blwh_names,".grd")]

blwh_all_raster <- stack() # Stack of 155 layers (31 years)
for (i in blwh_names) {
  raster <- raster(paste0(blwh_path,"/",i))
  blwh_all_raster <- stack(blwh_all_raster,raster)
}
blwh_all_raster <- setZ(blwh_all_raster, dates)
names(blwh_all_raster) <- dates

# Function to extract krill values from nc file ---------------------------

extract_fn <- function(df) {
  
  # Create new column for krill
  df$krill <- NA
  # Krill nc file
  totalkrill <- nc_open("~/Dropbox/krill_raw/TotalKril_CPUE.nc", write = TRUE)
  TotalKrill_CPUE <- ncvar_get(totalkrill, "TotalKrill_CPUE")
  fillvalue <- ncatt_get(totalkrill, "TotalKrill_CPUE", "_FillValue")
  TotalKrill_CPUE[TotalKrill_CPUE == fillvalue$value] <- NA
  dates <- format(as.Date(as.POSIXct(ncvar_get(totalkrill, "Time"),tz = "GMT", origin = "1970-01-01")),"%Y-%m")
  lon <- ncvar_get(totalkrill, "Longitude")
  lat <- ncvar_get(totalkrill, "Latitude", verbose = F)
  
  for (i in 1:nrow(df)) {
    
    # Finding matching date
    time <- paste0(df$year[i],"-0",df$month[i])
    layer_idx <- grep(time, dates, value = FALSE)
    
    if (length(layer_idx)==0) { # If date is outside April-August
      
      df$krill[i] <- NA 
    
      } else {
      
      # Lon and lat
      c <- which.min(abs(lon-df$lon[i]))
      r <- which.min(abs(lat-df$lat[i]))
      df$krill[i] <- ncvar_get(totalkrill,'TotalKrill_CPUE',start=c(c,r,layer_idx), count=c(1,1,1),verbose=FALSE)
    
      }
  }
  
  output <- df # Output is df + krill column added
  
}


# Extracting for all 40 presabs datasets -----------------------------------

blwh_ROMS_krill <- blwh_ROMS
row_num <- matrix(NA,1,40)
for (i in 1:40) {
  blwh_ROMS_krill[[i]] <- extract_fn(blwh_ROMS[[i]])
  row_num[i] <-   nrow(blwh_ROMS_krill[[i]][complete.cases(blwh_ROMS_krill[[i]]$krill),])
}

hist(row_num,xlab="Number of data points with krill")
abline(v = mean(row_num),col = "red",lwd = 3)
