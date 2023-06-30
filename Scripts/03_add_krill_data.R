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
blwh_ROMS_krillchl_filtered <- list() # Only data points with krill values
for (i in 1:40) {
  blwh_ROMS_krill[[i]] <- extract_fn(blwh_ROMS[[i]])
  blwh_ROMS_krillchl_filtered[[i]] <- blwh_ROMS_krill[[i]][rowSums(is.na(blwh_ROMS_krill[[i]][c(1,54)])) == 0,][c(1:2,5:8,50,54)]
}

# General statistics on krill data points ----------------------------------

row_num <- matrix(NA,1,40)
row_num_krill <- matrix(NA,1,40)
row_num_spatial <- matrix(NA,1,40) # Data points Apr-Aug including those without krill
for (i in 1:40) {
  row_num[i] <- nrow(blwh_ROMS_krill[[i]])
  row_num_krill[i] <-   nrow(blwh_ROMS_krillchl_filtered[[i]])
  row_num_spatial[i] <- nrow(filter(blwh_ROMS_krill[[i]],month==4|month==5|month==6|month==7|month==8))
}

# Blwh data points with krill (all 40 presabs dataset)
avg_total_points <- round(mean(row_num),0)

hist(row_num_krill,xlab="Data points",main="Total Number of blwh data points with krill")
abline(v = mean(row_num_krill),col = "red",lwd = 3)
text(x=3325,y=12,labels=paste0("Out of total ",avg_total_points," points"),col = "red")

# Fraction of blwh data points in Apr-Aug that have krill (all 40 presabs dataset)
fraction_num_spatial <- row_num_krill/row_num_spatial
percent_missed_spatial <- round(100*(1-mean(fraction_num_spatial)),2)

hist(fraction_num_spatial,xlab="Data points",main="Fraction of Apr-Aug blwh data points that have krill")
abline(v = mean(fraction_num_spatial),col = "red",lwd = 3)
text(x=0.715,y=8,labels=paste0("Avg points missed = ",percent_missed_spatial,"%"),col = "red")

# Fraction of blwh data points that have krill (all 40 presabs dataset)
fraction_num_temporal <- row_num_krill/row_num
percent_missed_temporal <- round(100*(1-mean(fraction_num_temporal)),2)

hist(fraction_num_temporal,xlab="Data points",main="Fraction of blwh data points that have krill")
abline(v = mean(fraction_num_temporal),col = "red",lwd = 3)
text(x=0.145,y=8,labels=paste0("Avg points missed = ",percent_missed_temporal,"%"),col = "red")


# Comparing pres and abs krill and chl statistics -------------------------

blwh_ROMS_1 <- blwh_ROMS_krill[[1]]
blwh_krillchl_1 <- blwh_ROMS_1[complete.cases(blwh_ROMS_krill[[i]]$krill),][c(1:2,5:8,50,54)]
blwh_pres_1 <- filter(blwh_krillchl_1,presabs==1)
blwh_abs_1 <- filter(blwh_krillchl_1,presabs!=1)
blwh_pres_num <- blwh_pres_1[complete.cases(blwh_pres_1$krill),]
blwh_abs_num <- blwh_abs_1[complete.cases(blwh_abs_1$krill),]

# Krill
ggplot() +
  geom_area(data=blwh_pres_1, aes(x=krill), stat = "bin",binwidth=.2,
            color= "blue",fill="lightblue",alpha=0.7) +
  geom_area(data=blwh_abs_1, aes(x=krill), stat = "bin",binwidth=.2,
            color= "red",fill="pink",alpha=0.7) +
  legend
  
  


  

# Plot
layer(sp.points(xy, pch=ifelse(xy$z1 < 0.5, 2, 3), cex=2, col=1), columns=1)
