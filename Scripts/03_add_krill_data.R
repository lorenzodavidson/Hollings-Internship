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
 
# Function to extract daily krill values with month centered around 16th ----

extract16th_fn <- function(df) {
  
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
    
    # Finding matching year-month
    year_month <- paste0(df$year[i],"-0",df$month[i])
    layer_idx <- grep(year_month, dates, value = FALSE)
    
    if (length(layer_idx)==0) { # If date is outside April-August
      
      df$krill[i] <- NA 
      
    } else {
      
      # Lon and lat
      c <- which.min(abs(lon-df$lon[i]))
      r <- which.min(abs(lat-df$lat[i]))
      month <- paste0(df$month[i])
      day <- as.numeric(df$day[i])
      
      if (month == "4") { # April (no previous month)
        
        if (day <= 15) {
          
          df$krill[i] <- ncvar_get(totalkrill,'TotalKrill_CPUE',start=c(c,r,layer_idx), count=c(1,1,1),verbose=FALSE)
        
        } else {
          
          # Weighted mean
          weight1 <- 1-(day-16)/30
          weight2 <- 1-weight1
          krill1 <- weight1*ncvar_get(totalkrill,'TotalKrill_CPUE',start=c(c,r,layer_idx), count=c(1,1,1),verbose=FALSE)
          krill2 <- weight2*ncvar_get(totalkrill,'TotalKrill_CPUE',start=c(c,r,layer_idx+1), count=c(1,1,1),verbose=FALSE)
          df$krill[i] <- krill1+krill2
        
        }
        
      } else if (month == "5" | month == "7") { # May and July have 31 days
        
        if (day <= 15) {
          
          # Weighted mean with previous month
          weight1 <- 1-abs(day-16)/30
          weight2 <- 1-weight1
          krill1 <- weight1*ncvar_get(totalkrill,'TotalKrill_CPUE',start=c(c,r,layer_idx), count=c(1,1,1),verbose=FALSE)
          krill2 <- weight2*ncvar_get(totalkrill,'TotalKrill_CPUE',start=c(c,r,layer_idx-1), count=c(1,1,1),verbose=FALSE)
          df$krill[i] <- krill1+krill2
          
        } else {
          
          # Weighted mean
          weight1 <- 1-abs(day-16)/31
          weight2 <- 1-weight1
          krill1 <- weight1*ncvar_get(totalkrill,'TotalKrill_CPUE',start=c(c,r,layer_idx), count=c(1,1,1),verbose=FALSE)
          krill2 <- weight2*ncvar_get(totalkrill,'TotalKrill_CPUE',start=c(c,r,layer_idx+1), count=c(1,1,1),verbose=FALSE)
          df$krill[i] <- krill1+krill2
          
        } 
      
      } else if (month == "6") {
          
          if (day <= 15) {
            
            # Weighted mean with previous month
            weight1 <- 1-abs(day-16)/31
            weight2 <- 1-weight1
            krill1 <- weight1*ncvar_get(totalkrill,'TotalKrill_CPUE',start=c(c,r,layer_idx), count=c(1,1,1),verbose=FALSE)
            krill2 <- weight2*ncvar_get(totalkrill,'TotalKrill_CPUE',start=c(c,r,layer_idx-1), count=c(1,1,1),verbose=FALSE)
            df$krill[i] <- krill1+krill2
            
          } else {
            
            # Weighted mean
            weight1 <- 1-abs(day-16)/30
            weight2 <- 1-weight1
            krill1 <- weight1*ncvar_get(totalkrill,'TotalKrill_CPUE',start=c(c,r,layer_idx), count=c(1,1,1),verbose=FALSE)
            krill2 <- weight2*ncvar_get(totalkrill,'TotalKrill_CPUE',start=c(c,r,layer_idx+1), count=c(1,1,1),verbose=FALSE)
            df$krill[i] <- krill1+krill2
            
          }
        
      } else { # August has no month after
        
        if (day > 15) {
          
          df$krill[i] <- ncvar_get(totalkrill,'TotalKrill_CPUE',start=c(c,r,layer_idx), count=c(1,1,1),verbose=FALSE)
          
        } else {
          
          # Weighted mean
          weight1 <- 1-abs(day-16)/30
          weight2 <- 1-weight1
          krill1 <- weight1*ncvar_get(totalkrill,'TotalKrill_CPUE',start=c(c,r,layer_idx), count=c(1,1,1),verbose=FALSE)
          krill2 <- weight2*ncvar_get(totalkrill,'TotalKrill_CPUE',start=c(c,r,layer_idx-1), count=c(1,1,1),verbose=FALSE)
          df$krill[i] <- krill1+krill2
          
        }
      }
    }
  }
 
  output <- df # Output is df + krill column added
  
}

# Function to extract daily krill values with month centered around 1st ----
extract1st_fn <- function(df) {
  
  # Create new column for krill
  df$krill <- NA
  # Krill nc file
  totalkrill <- nc_open("~/Dropbox/krill_raw/TotalKril_CPUE.nc", write = TRUE)
  TotalKrill_CPUE <- ncvar_get(totalkrill, "TotalKrill_CPUE")
  fillvalue <- ncatt_get(totalkrill, "TotalKrill_CPUE", "_FillValue")
  TotalKrill_CPUE[TotalKrill_CPUE == fillvalue$value] <- NA
  dates <- format(as.Date(as.POSIXct(ncvar_get(totalkrill, "Time"),tz = "GMT", origin = "1970-01-01")),"%Y-%m-%d")
  lon <- ncvar_get(totalkrill, "Longitude")
  lat <- ncvar_get(totalkrill, "Latitude", verbose = F)
  
  for (i in 1:nrow(df)) {
    
    # Finding matching year-month
    year_month <- paste0(df$year[i],"-0",df$month[i])
    layer_idx <- grep(year_month, dates, value = FALSE)
    
    if (length(layer_idx)==0) { # If date is outside April-August
      
      df$krill[i] <- NA 
      
    } else {
      
      # Lon and lat
      c <- which.min(abs(lon-df$lon[i]))
      r <- which.min(abs(lat-df$lat[i]))
      month <- paste0(df$month[i])
      
      if (month == "4" | month == "6") { # April and June have 30 days
        
        # Weighted mean
        weight1 <- 1-(df$day-1)/30
        weight2 <- 1-weight1
        krill1 <- weight1*ncvar_get(totalkrill,'TotalKrill_CPUE',start=c(c,r,layer_idx), count=c(1,1,1),verbose=FALSE)
        krill2 <- weight2*ncvar_get(totalkrill,'TotalKrill_CPUE',start=c(c,r,layer_idx+1), count=c(1,1,1),verbose=FALSE)
        df$krill[i] <- krill1+krill2
        
      } else if (month == "5" | month == "7") { # May and July have 31 days
        
        # Weighted mean
        weight1 <- 1-(df$day-1)/31
        weight2 <- 1-weight1
        krill1 <- weight1*ncvar_get(totalkrill,'TotalKrill_CPUE',start=c(c,r,layer_idx), count=c(1,1,1),verbose=FALSE)
        krill2 <- weight2*ncvar_get(totalkrill,'TotalKrill_CPUE',start=c(c,r,layer_idx+1), count=c(1,1,1),verbose=FALSE)
        df$krill[i] <- krill1+krill2
        
      } else { # August has no month after
        
        df$krill[i] <- ncvar_get(totalkrill,'TotalKrill_CPUE',start=c(c,r,layer_idx), count=c(1,1,1),verbose=FALSE)
          
      }
    }
  }
  
  output <- df # Output is df + krill column added
  
}


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

# Turn month column into character
for (i in 1:40) {
  blwh_ROMS[[i]]$month <- as.character(blwh_ROMS[[i]]$month)
}

blwh_ROMS_krill <- blwh_ROMS
blwh_ROMS_filtered <- list() # Only data points with krill values
for (i in 1:40) {
  blwh_ROMS_krill[[i]] <- extract_fn(blwh_ROMS[[i]])
  blwh_ROMS_filtered[[i]] <- blwh_ROMS_krill[[i]][rowSums(is.na(blwh_ROMS_krill[[i]][c(1,54)])) == 0,][c(1:2,5:9,16,50,54)]
}

# Extract daily krill values with 2 functions
blwh_ROMS_krill_1st <- blwh_ROMS
blwh_ROMS_filtered_1st <- list() # Only data points with krill values
blwh_ROMS_krill_16th <- blwh_ROMS
blwh_ROMS_filtered_16th <- list() # Only data points with krill values
for (i in 1:40) {
  blwh_ROMS_krill_1st[[i]] <- extract1st_fn(blwh_ROMS[[i]])
  blwh_ROMS_filtered_1st[[i]] <- blwh_ROMS_krill_1st[[i]][rowSums(is.na(blwh_ROMS_krill_1st[[i]][c(1,54)])) == 0,][c(1:2,5:9,16,50,54)]

  blwh_ROMS_krill_16th[[i]] <- extract16th_fn(blwh_ROMS[[i]])
  blwh_ROMS_filtered_16th[[i]] <- blwh_ROMS_krill_16th[[i]][rowSums(is.na(blwh_ROMS_krill_16th[[i]][c(1,54)])) == 0,][c(1:2,5:9,16,50,54)]
}

check_1st <- blwh_ROMS_filtered_1st[[1]]
check_16th <- blwh_ROMS_filtered_16th[[1]]

# Compare 3 different ways of creating daily krill ------------------------

for (i in 1:40) {
  if (i ==1) {
    blwh_krill1_all <- blwh_ROMS_filtered[[i]]
    blwh_krill2_all <- blwh_ROMS_filtered_1st[[i]]
    blwh_krill3_all <- blwh_ROMS_filtered_16th[[i]]
  } else {
    blwh_krill1_all <- rbind(blwh_krill1_all,blwh_ROMS_filtered[[i]])
    blwh_krill2_all <- rbind(blwh_krill2_all,blwh_ROMS_filtered_1st[[i]])
    blwh_krill3_all <- rbind(blwh_krill3_all,blwh_ROMS_filtered_16th[[i]])
  }
}

blwh_krill_all <- cbind(blwh_krill1_all,blwh_krill2_all$krill,blwh_krill3_all$krill)
colnames(blwh_krill_all)[11] <- "krill2"
colnames(blwh_krill_all)[12] <- "krill3"

mean_krill <- mean(blwh_krill_all$krill,na.rm=TRUE) # 7.14
mean_krill2 <- mean(blwh_krill_all$krill2,na.rm=TRUE) # 7.13
mean_krill3 <- mean(blwh_krill_all$krill3,na.rm=TRUE) # 7.14

ggplot() + 
  geom_density(data=blwh_krill_all, aes(x=krill),alpha=0.4,fill="lightblue") +
  geom_density(data=blwh_krill_all, aes(x=krill2),alpha=0.4,fill="pink") +
  geom_density(data=blwh_krill_all, aes(x=krill3),alpha=0.4,fill="lightgreen") +
  ggtitle("Krill Distribution for 3 Ways of Creating Daily Krill Data")

# General statistics on krill data points for the 40 sets -------------------

row_num <- matrix(NA,1,40)
row_num_krill <- matrix(NA,1,40)
row_num_spatial <- matrix(NA,1,40) # Data points Apr-Aug including those without krill
for (i in 1:40) {
  row_num[i] <- nrow(blwh_ROMS_krill[[i]])
  row_num_krill[i] <-   nrow(blwh_ROMS_krillchl_filtered[[i]])
  row_num_spatial[i] <- nrow(filter(blwh_ROMS_krill[[i]],month==4|month==5|month==6|month==7|month==8))
}

blwh_krill_1 <- blwh_ROMS_krillchl_filtered[[1]]

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

# Comparing pres and abs krill, chl and sst statistics ----------------------

# Using blwh_krill_all
blwh_krill_all$presabs[blwh_krill_all$presabs == 1] <- "Presence"
blwh_krill_all$presabs[blwh_krill_all$presabs == 0] <- "Absence"

library(plyr)
mu_krill <- ddply(blwh_krill_all, "presabs", summarise, grp.mean=mean(krill))
mu_krill2 <- ddply(blwh_krill_all, "presabs", summarise, grp.mean=mean(krill2))
mu_krill3 <- ddply(blwh_krill_all, "presabs", summarise, grp.mean=mean(krill3))
mu_chl <- ddply(blwh_krill_all, "presabs", summarise, grp.mean=mean(chl,na.rm=TRUE))
mu_chl[1,2] <- log(mu_chl[1,2]); mu_chl[2,2] <- log(mu_chl[2,2])
mu_sst <- ddply(blwh_krill_all, "presabs", summarise, grp.mean=mean(sst_mean_0.1,na.rm=TRUE))


# Krill
ggplot(data=blwh_krill_all, aes(x=krill, fill=presabs)) + 
  geom_density(alpha=0.4) +
  geom_vline(data=mu_krill, aes(xintercept=grp.mean, color=presabs),
             linetype="dashed",linewidth=1) + 
  ggtitle("Krill Distribution for Presence and Absence Data")

# Krill2
ggplot(data=blwh_krill_all, aes(x=krill2, fill=presabs)) + 
  geom_density(alpha=0.4) +
  geom_vline(data=mu_krill2, aes(xintercept=grp.mean, color=presabs),
             linetype="dashed",linewidth=1) + 
  ggtitle("Krill Distribution for Presence and Absence Data")

# Krill3
ggplot(data=blwh_krill_all, aes(x=krill3, fill=presabs)) + 
  geom_density(alpha=0.4) +
  geom_vline(data=mu_krill3, aes(xintercept=grp.mean, color=presabs),
             linetype="dashed",linewidth=1) + 
  ggtitle("Krill Distribution for Presence and Absence Data")

# Chl
ggplot(data=blwh_krill_all, aes(x=log(chl), fill=presabs)) + 
  geom_density(alpha=0.4) +
  geom_vline(data=mu_chl, aes(xintercept=grp.mean, color=presabs),
             linetype="dashed",linewidth=1) + 
  ggtitle("Chlorophyll Distribution for Presence and Absence Data")

# SST
ggplot(data=blwh_krill_all, aes(x=sst_mean_0.1, fill=presabs)) + 
  geom_density(alpha=0.4) +
  geom_vline(data=mu_sst, aes(xintercept=grp.mean, color=presabs),
             linetype="dashed",linewidth=1) + 
  ggtitle("Chlorophyll Distribution for Presence and Absence Data")

  

# Plot
layer(sp.points(xy, pch=ifelse(xy$z1 < 0.5, 2, 3), cex=2, col=1), columns=1)
