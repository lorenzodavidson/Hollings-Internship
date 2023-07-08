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
}))


# Open files --------------------------------------------------------------

processed_path <- "~/Desktop/Hollingsinternship/Hollings-Internship/Processed/"
bwkr_all <- read.csv(paste0(processed_path,"/bwkr_all.csv"))
dates <- t(read.csv(paste0(processed_path,"/dates.csv")))

# Created stacked rasters for complete set

blwh_path <- paste0(processed_path,"blwh_all/")
krill_path <- paste0(processed_path,"krill_all/")
blwh_names <- list.files(blwh_path)
krill_names <- list.files(krill_path)
blwh_names <- list.files(blwh_path)[stringr::str_detect(blwh_names,".grd")]
krill_names <- list.files(krill_path)[stringr::str_detect(krill_names,".grd")]

blwh_all_raster <- stack() # Stack of 155 layers (31 years)
for (i in blwh_names) {
  raster <- raster(paste0(blwh_path,"/",i))
  blwh_all_raster <- stack(blwh_all_raster,raster)
}
blwh_all_raster <- setZ(blwh_all_raster, dates)

krill_all_raster <- stack()
for (i in krill_names) {
  raster <- raster(paste0(krill_path,"/",i))
  krill_all_raster <- stack(krill_all_raster,raster)
}
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

# Freq plot by month -----------------------------------------------------

months <- c('Apr', 'May', 'Jun', 'Jul', 'Aug')
bwkr_bymonth <- data.frame(group_by(bwkr_all,month))
bwkr_bymonth$month <- as.character(bwkr_bymonth$month)

ggplot(bwkr_bymonth, aes(blwh,color=month)) +
  geom_freqpoly(binwidth = 0.01,linewidth=1) +
  xlim(0, 1) +
  ylim(0,20000) +
  ggtitle("Monthly Blue Whale HS Frequency") +
  labs(x = "Blwh Habitat Suitability",
       y = "Count",
       color = "Months")

ggplot(bwkr_bymonth, aes(krill,color=month)) +
  geom_freqpoly(binwidth = 0.03,linewidth=1) +
  xlim(range(bwkr_bymonth$krill)[[1]], range(bwkr_bymonth$krill)[[2]]) +
  ggtitle("Monthly Krill CPUE Frequency") +
  labs(x = "Krill Catch per Unit Effort (CPUE)",
       y = "Count",
       color = "Months")

# Coastline shapefile for plots -------------------------------------------

fpath_ne <- "~/Desktop/Hollingsinternship/Hollings-Internship/Processed/Coastline/ne_10m_land/"
dfname_ne <- "ne_10m_land.shp"
neland_data <- read_sf(paste0(fpath_ne,dfname_ne))
world_land <- as(st_geometry(neland_data), Class="Spatial")

# Mean Blwh HS and Krill CPUE for all months -------------------------------
months <- c("April","May","June","July","August")

blwh_apr_mean <- mean(blwh_apr_raster)
blwh_may_mean <- mean(blwh_may_raster)
blwh_jun_mean <- mean(blwh_jun_raster)
blwh_jul_mean <- mean(blwh_jul_raster)
blwh_aug_mean <- mean(blwh_aug_raster)
blwh_all_mean <- stack(blwh_apr_mean,blwh_may_mean,blwh_jun_mean,blwh_jul_mean,blwh_aug_mean)
blwh_all_mean <- setZ(blwh_all_mean,months)
names(blwh_all_mean) <- months

krill_apr_mean <- mean(krill_apr_raster)
krill_may_mean <- mean(krill_may_raster)
krill_jun_mean <- mean(krill_jun_raster)
krill_jul_mean <- mean(krill_jul_raster)
krill_aug_mean <- mean(krill_aug_raster)
krill_all_mean <- stack(krill_apr_mean,krill_may_mean,krill_jun_mean,krill_jul_mean,krill_aug_mean)
names(krill_all_mean) <- months

# Breaks for legend based on percentiles
blwh_breaks <- unname(quantile(bwkr_all$blwh,probs = c(0.35,0.45,0.55,0.65,0.75,0.85,0.95,1)))
blwh_breaks <- c(0, blwh_breaks)
krill_breaks <- unname(quantile(bwkr_all$krill,probs = c(0.35,0.45,0.55,0.65,0.75,0.85,0.95,1)))
krill_breaks <- c(0,krill_breaks)

myTheme <- rasterTheme(region = brewer.pal(9,"YlOrRd"),panel.background = list(col = 'lightblue'))

# Blue Whale
plt <- levelplot(blwh_all_mean,main = "Mean Monthly Blue Whale HS",par.settings = myTheme,at = blwh_breaks)
plt + latticeExtra::layer(sp.polygons(world_land, fill="gray80"))

# Krill
plt <- levelplot(krill_all_mean,main = "Mean Monthly Krill CPUE",par.settings = myTheme,at = krill_breaks)
plt + latticeExtra::layer(sp.polygons(world_land, fill="gray80"))

# June, July and August blue whale and krill comparison --------------------

months <- c('Apr', 'May', 'Jun', 'Jul', 'Aug')
bwkr_bymonth <- data.frame(group_by(bwkr_all,month))
bwkr_bymonth$month <- as.character(bwkr_bymonth$month)

blwh_678_mean <- subset(blwh_all_mean,3:5)
krill_678_mean <- subset(krill_all_mean,3:5)

blwh_breaks <- unname(quantile(filter(bwkr_bymonth, month == "6" | month == "7" | month == "8")$blwh, 
                               probs = c(0.35,0.45,0.55,0.65,0.75,0.85,0.95,1)))
blwh_breaks <- c(0, blwh_breaks)

krill_breaks <- unname(quantile(filter(bwkr_bymonth, month == "6" | month == "7" | month == "8")$krill, 
                                probs = c(0.35,0.45,0.55,0.65,0.75,0.85,0.95,1)))
krill_breaks <- c(0,krill_breaks)

myTheme <- rasterTheme(region = brewer.pal(9,"YlOrRd"),panel.background = list(col = 'lightblue'))

library(gridExtra)
p1 <- levelplot(blwh_678_mean,main = "Mean Monthly Blue Whale HS",par.settings = myTheme,at = blwh_breaks) + 
      latticeExtra::layer(sp.polygons(world_land, fill="gray80"))
p2 <- levelplot(krill_678_mean,main = "Mean Monthly Krill CPUE",par.settings = myTheme,at = krill_breaks) + 
      latticeExtra::layer(sp.polygons(world_land, fill="gray80"))
grid.arrange(p1, p2, ncol=1)

# Reoccurrence of blue whale core habitat ---------------------------------

# Function to filter for months and threshold, group by xy, and calculate % of years
thresh_fn <- function(df, time, thresh) {
  df_count <- df %>%
    filter(month == time) %>% # time = month of interest ("6","7" or "8")
    filter(blwh > thresh) %>%
    group_by(lon,lat) %>%
    mutate(n=n()/0.31) %>%
    arrange(desc(n))
  output <- df_count[!duplicated(df_count[,c("lon","lat")]),]
}

# Creating rasters with proper extension
ext <- extent(-127.5, -115.5, 30, 48) # Extent of bwkr raster with edges
raster_fn <- function(df, ext) {
  extend(rasterFromXYZ(df[c(-3:-8)]),ext)
}

blwh_percentiles <- unname(quantile(filter(bwkr_bymonth, month == "6" | month == "7" | month == "8")$blwh, 
                probs = c(0.55,0.65,0.75,0.85,0.95)))

blwh_jun_count <- stack()
blwh_jul_count <- stack()
blwh_aug_count <- stack()
for (i in 1:length(blwh_percentiles)) {
  blwh_jun_prep <- raster_fn(thresh_fn(bwkr_bymonth,"6",blwh_percentiles[[i]]),ext)
  blwh_jun_count <- stack(blwh_jun_count,blwh_jun_prep)
  
  blwh_jul_prep <- raster_fn(thresh_fn(bwkr_bymonth,"7",blwh_percentiles[[i]]),ext)
  blwh_jul_count <- stack(blwh_jul_count,blwh_jul_prep)
  
  blwh_aug_prep <- raster_fn(thresh_fn(bwkr_bymonth,"8",blwh_percentiles[[i]]),ext)
  blwh_aug_count <- stack(blwh_aug_count,blwh_aug_prep)
}

thresholds <- c("55th Percentile","65th Percentile","75th Percentile","85th Percentile","95th Percentile")
blwh_jun_count <- setZ(blwh_jun_count, paste0(thresholds))
blwh_jul_count <- setZ(blwh_jul_count, paste0(thresholds))
blwh_aug_count <- setZ(blwh_aug_count, paste0(thresholds))

myTheme <- rasterTheme(region = brewer.pal(9,"YlOrRd"),panel.background = list(col = 'lightblue'))

# June
plt <- levelplot(blwh_jun_count,main = "Percent Junes above Blue Whale HS Threshold",
                 names.attr=as.character(thresholds), par.settings = myTheme)
plt + latticeExtra::layer(sp.polygons(world_land, fill="gray80"))

# July
plt <- levelplot(blwh_jul_count,main = "Percent July above Blue Whale HS Threshold",
                 names.attr=as.character(thresholds), par.settings = myTheme)
plt + latticeExtra::layer(sp.polygons(world_land, fill="gray80"))

# Aug
plt <- levelplot(blwh_aug_count,main = "Percent Augusts above Blue Whale HS Threshold",
                 names.attr=as.character(thresholds), par.settings = myTheme)
plt + latticeExtra::layer(sp.polygons(world_land, fill="gray80"))

# Plot persistent blue whale core habitat above thresholds -----------------

# Function that turns pixels with <50% to NA and >50% to threshold value
library(terra)

blwh_percentiles <- unname(quantile(filter(bwkr_bymonth, month == "6" | month == "7" | month == "8")$blwh, 
                                    probs = c(0.55,0.65,0.75,0.85,0.95)))

binary_fn <- function(raster,percentiles) {
  raster_binary <- as(raster, "SpatRaster")
  raster_binary <- app(raster_binary, fun=function(x){ x[x <= 50] <- NA; return(x)})
  raster_binary <- app(raster_binary, fun=function(x){ x[x > 50] <- 1; return(x)})
  
  for (i in seq_along(blwh_percentiles)) {
    raster_binary[[i]] <- app(raster_binary[[i]], fun=function(x){ x[x == 1] <- percentiles[[i]]; return(x)})
  }
  
  output <- cover(raster_binary[[5]],raster_binary[[4]])
  output <- cover(output,raster_binary[[3]])
  output <- cover(output,raster_binary[[2]])
  output <- cover(output,raster_binary[[1]])
  output <- raster(output)
}

# Creating 1 raster with persistent hot spots above each threshold

blwh_jun_binary <- binary_fn(blwh_jun_count,blwh_percentiles)
blwh_jul_binary <- binary_fn(blwh_jul_count,blwh_percentiles)
blwh_aug_binary <- binary_fn(blwh_aug_count,blwh_percentiles)
blwh_678_binary <- stack(blwh_jun_binary,blwh_jul_binary,blwh_aug_binary)

blwh_breaks <- c((blwh_percentiles-0.01),1)

plt <- levelplot(blwh_678_binary,main = "Persistent Core Blue Whale Habitat",
                  par.settings = myTheme, names.attr=as.character(c("June","July","August")), at = blwh_breaks)
plt + latticeExtra::layer(sp.polygons(world_land, fill="gray80"))

# Function to count pixels of persistent core habitat above each threshold

count_fn <- function(raster) {
  raster_binary <- as(raster, "SpatRaster")
  raster_binary <- app(raster_binary, fun=function(x){ x[x <= 50] <- NA; return(x)})
  raster_binary <- app(raster_binary, fun=function(x){ x[x > 50] <- 1; return(x)})
  
  count <- matrix(NA,1,5)
  for (i in 1:5) {
    raster_binary_layer <- raster(raster_binary[[i]])
    cnt <- data.frame(rasterToPoints(raster_binary_layer)) %>%
      summarize(count=n())
    count[i] <- unlist(cnt)
  }
  output <- data.frame(count)
}

blwh_jun_CHcount <- count_fn(blwh_jun_count)
blwh_jul_CHcount <- count_fn(blwh_jul_count)
blwh_aug_CHcount <- count_fn(blwh_aug_count)

# Finish plotting lines for each threshold

# Overlap Metrics Code ----------------------------------------------------

# Area Overlap
# for binary data
# measures proportion of an area where two species co-occur
area_overlapfn <- function(prey, pred, area){
  total_area <- sum(area, na.rm = T)
  sum(area[pred > 0 & prey > 0], na.rm = T)/total_area
}

# Range Overlap
# for binary data
# measures the proportion of one species range where the other co-occurs
range_overlapfn<-function(prey, pred, area){
  area_prey <- sum(area[prey > 0], na.rm = T)
  sum(area[pred > 0 & prey > 0], na.rm = T)/area_prey
}

# Schoener's D
# density or probability of occurrence data
# measures how equally predator and prey share available resources
schoeners_overlapfn <- function(prey, pred) {
  p_prey <- prey/sum(prey, na.rm = T)
  p_pred <- pred/sum(pred, na.rm = T)
  1 - 0.5 * (sum(abs(p_prey-p_pred), na.rm = T))
}

# Bhattacharyya's coefficient
# density or probability of occurrence data
# measures whether two species use space independently
bhatta_coeffn <- function(prey, pred) {
  p_prey <- prey/sum(prey, na.rm = T)
  p_pred <- pred/sum(pred, na.rm = T)
  sum(sqrt(p_prey*p_pred), na.rm = T)
}

# Center of gravity (to be developed)

# Area and Range Overlap Metrics for 9 different threshold combinations ----

blwh_thresh <- unname(quantile(bwkr_all$blwh, probs = c(.65, .75, .85)))
krill_thresh <- unname(quantile(bwkr_all$krill, probs = c(.65, .75, .85)))

bwkr_metrics_9combo <- bwkr_all
bwkr_metrics_9combo <- bwkr_all %>% 
  mutate(blwh_core_65 = ifelse(blwh >= blwh_thresh[[1]],1,0), 
         blwh_core_75 = ifelse(blwh >= blwh_thresh[[2]],1,0),
         blwh_core_85 = ifelse(blwh >= blwh_thresh[[3]],1,0),
         krill_core_65 = ifelse(krill >= krill_thresh[[1]],1,0),
         krill_core_75 = ifelse(krill >= krill_thresh[[2]],1,0),
         krill_core_85 = ifelse(krill >= krill_thresh[[3]],1,0), Area=1)

blwh_cores <- c("blwh_core_65","blwh_core_75","blwh_core_85")
krill_cores <- c("krill_core_65","krill_core_75","krill_core_85")
thresholds <- c("65","75","85")

# AO Plots
par(mfrow=c(3,3))
for (i in 1:3) { # blwh_tresh
  for (j in 1:3) { # krill_tresh
    x <- 8 + i # Column of blwh_core desired
    y <- 11 + j # Column of krill_core desired
    AO_plot <- bwkr_metrics_9combo[c(1:6,8,x,y,15)] %>%
      rename("blwh_core"=paste0(blwh_cores[[i]]),"krill_core"=paste0(krill_cores[[j]])) %>%
      group_by(year,month) %>% # organize by date
      summarise(AO=area_overlapfn(prey=krill_core,pred=blwh_core,area=Area)) %>% 
      pivot_longer(cols=c(AO),names_to="Metric_Name",values_to="Overlap_Metric")
    boxplot(Overlap_Metric~month,
            data=AO_plot,
            xlab="Month",
            ylab="Area Overlap",
            # main=paste0("blwh thresh: ", thresholds[[j]]," - ","krill thresh: ",thresholds[[i]]),
            ylim = c(0, 0.25))
  }
}

# RO Plots
par(mfrow=c(3,3))
for (i in 1:3) { # blwh_tresh
  for (j in 1:3) { # krill_tresh
    x <- 8 + i # Column of blwh_core desired
    y <- 11 + j # Column of krill_core desired
    RO_plot <- bwkr_metrics_9combo[c(1:6,8,x,y,15)] %>%
      rename("blwh_core"=paste0(blwh_cores[[i]]),"krill_core"=paste0(krill_cores[[j]])) %>%
      group_by(year,month) %>% # organize by date
      summarise(RO=range_overlapfn(prey=krill_core,pred=blwh_core,area=Area)) %>% 
      pivot_longer(cols=c(RO),names_to="Metric_Name",values_to="Overlap_Metric")
    boxplot(Overlap_Metric~month,
            data=RO_plot,
            xlab="Month",
            ylab="Range Overlap",
            # main=paste0(thresholds[[j]],"-",thresholds[[i]]),
            ylim = c(0, 1))
  }
}

dev.off()

# Metrics for thresh: blwh 0.44 + krill 75% -------------------------------

# Finding appropriate threshold for core habitat and filtering data
bwkr_bymonth <- data.frame(group_by(bwkr_all,month))
bwkr_bymonth$month <- as.character(bwkr_bymonth$month)

blwh_percentiles <- unname(quantile(bwkr_all$blwh,probs = c(0.35,0.45,0.55,0.65,0.75,0.85,0.95)))
krill_percentiles <- unname(quantile(bwkr_all$krill,probs = c(0.35,0.45,0.55,0.65,0.75,0.85,0.95)))
blwh_jun_percentiles <- unname(quantile(filter(bwkr_bymonth,month == 6)$blwh,probs = c(0.35,0.45,0.55,0.65,0.75,0.85,0.95)))
blwh_jul_percentiles <- unname(quantile(filter(bwkr_bymonth,month == 7)$blwh,probs = c(0.35,0.45,0.55,0.65,0.75,0.85,0.95)))
blwh_aug_percentiles <- unname(quantile(filter(bwkr_bymonth,month == 8)$blwh,probs = c(0.35,0.45,0.55,0.65,0.75,0.85,0.95)))
blwh_678_percentiles <- unname(quantile(filter(bwkr_bymonth,month == 6 | month == 7 | month == 8)$blwh,
                                        probs = c(0.35,0.45,0.55,0.65,0.75,0.85,0.95)))

plot(x=c(0.35,0.45,0.55,0.65,0.75,0.85,0.95),y=blwh_percentiles,main="Blwh HS vs Percentile",
     xlab="Percentile", ylab="Blwh HS", ylim=c(0,1),xlim=c(0.3,1))
lines(x=c(0.35,0.45,0.55,0.65,0.75,0.85,0.95),y=blwh_jun_percentiles,col="red")
lines(x=c(0.35,0.45,0.55,0.65,0.75,0.85,0.95),y=blwh_jul_percentiles,col="orange")
lines(x=c(0.35,0.45,0.55,0.65,0.75,0.85,0.95),y=blwh_aug_percentiles,col="green")
lines(x=c(0.35,0.45,0.55,0.65,0.75,0.85,0.95),y=blwh_678_percentiles,col="blue")
legend(0.3, 0.95, legend=c("All Months", "June", "July","August","Jun/Jul/Aug"),
       col=c("black","red","orange","green","blue"),lty=1:2, cex=0.8)


# Matching krill threshold to % blwh above 0.28 threshold and calculating metrics
blwh_044 <- filter(bwkr_all,bwkr_all$blwh>0.44) # top ~75th percentile of all data, top 60th Jun/Jul/Aug
krill_75p_thresh <- unname(quantile(bwkr_all$krill, 0.75))
krill_75p <- filter(bwkr_all,bwkr_all$krill>krill_75p_thresh) 

# For some reason not working anymore ---
bwkr_metrics <- bwkr_all %>% # make new dataframe called dailyzoom2015 which is a copy of sdm2015 and do the following below
  mutate(blwh_core = ifelse(blwh >= 0.44,1,0), krill_core = ifelse(krill >= krill_75p_thresh,1,0), Area=1) %>% # create a binary data column for range and area overlap of each species
  group_by(year,month) %>% # organize by date
  summarise(AO=area_overlapfn(prey=krill_core,pred=blwh_core,area=Area),
            RO=range_overlapfn(prey=krill_core,pred=blwh_core,area=Area),
            Schoener=schoeners_overlapfn(prey=krill,pred=blwh),
            Bhatty=bhatta_coeffn(prey=krill,pred=blwh))
pivot_longer(cols=c(AO,RO,Schoener,Bhatty),names_to="Metric_Name",values_to="Overlap_Metric") # change the dataframe into long format for easier plotting
# ---

processed_path <- "~/Desktop/Hollingsinternship/Hollings-Internship/Processed/"
dates <- t(read.csv(paste0(processed_path,"/dates.csv")))

bwkr_core <- bwkr_all %>% 
  mutate(blwh_core = ifelse(blwh >= 0.44,1,0), krill_core = ifelse(krill >= krill_75p_thresh,1,0), Area=1)
for (i in 1:155) {
  time <- dates[i]
  if (i == 1) {
    bwkr_metrics <- filter(bwkr_core,date==time) %>%
      summarize(AO=area_overlapfn(prey=krill_core,pred=blwh_core,area=Area),
                RO=range_overlapfn(prey=krill_core,pred=blwh_core,area=Area),
                Schoener=schoeners_overlapfn(prey=krill,pred=blwh),
                Bhatty=bhatta_coeffn(prey=krill,pred=blwh))
  } else {
    metrics <- filter(bwkr_core,date==time) %>%
      summarize(AO=area_overlapfn(prey=krill_core,pred=blwh_core,area=Area),
                RO=range_overlapfn(prey=krill_core,pred=blwh_core,area=Area),
                Schoener=schoeners_overlapfn(prey=krill,pred=blwh),
                Bhatty=bhatta_coeffn(prey=krill,pred=blwh))
    bwkr_metrics <- rbind(bwkr_metrics,metrics)
  }
}
bwkr_metrics <- cbind(t(dates),bwkr_metrics)
bwkr_metrics <- bwkr_metrics %>%
  separate(x,c('year','month'))
bwkr_metrics <- cbind(t(dates),bwkr_metrics)
colnames(bwkr_metrics)[1] <- "date"

# AO Plots
ggplot(bwkr_metrics,aes(x=date,y=AO,group = 1)) +
  geom_point() + 
  geom_line() +
  ggtitle("AO 1990-2020") +
  ylab("Area Overlap")
ggplot(bwkr_metrics,aes(x=month,y=AO)) +
  geom_boxplot() +
  ggtitle("AO grouped by month") +
  xlab("Month (March-July)") + 
  ylab("Area Overlap")

# RO Plots
ggplot(bwkr_metrics,aes(x=date,y=RO,group = 1)) +
  geom_point() + 
  geom_line() +
  ggtitle("RO 1990-2020") +
  ylab("Range Overlap")
ggplot(bwkr_metrics,aes(x=month,y=RO)) +
  geom_boxplot() +
  ggtitle("RO grouped by month") +
  xlab("Month (March-July)") + 
  ylab("Range Overlap")

# Schoener Plots
ggplot(bwkr_metrics,aes(x=date,y=Schoener,group = 1)) +
  geom_point() + 
  geom_line() +
  ggtitle("Schoener's D 1990-2020") +
  ylab("Schoener's D")
ggplot(bwkr_metrics,aes(x=month,y=Schoener,group=month)) +
  geom_boxplot() +
  ggtitle("Schoener's D grouped by month") +
  xlab("Month (March-July)") + 
  ylab("Schoener's D")

# Bhattacharyya Plots
ggplot(bwkr_metrics,aes(x=date,y=Bhatty,group = 1)) +
  geom_point() + 
  geom_line() +
  ggtitle("Bhattacharyya's Coefficient 1990-2020") +
  ylab("Bhattacharyya's Coefficient")
ggplot(bwkr_metrics,aes(x=month,y=Bhatty,group=month)) +
  geom_boxplot() +
  ggtitle("Bhattacharyya's Coefficient grouped by month") +
  xlab("Month (March-July)") + 
  ylab("Bhattacharyya's Coefficient")


# Visualizing blue whale and krill core habitat overlap -------------------

# Check for months with core habitat overlap
blwh_core <- bwkr_all %>%
  mutate(blwh_core = ifelse(blwh >= 0.44,1,0), krill_core = ifelse(krill >= krill_75p_thresh,1,0), Area=1) %>%
  filter(blwh_core == 1)

krill_core <- bwkr_all %>%
  mutate(blwh_core = ifelse(blwh >= 0.44,1,0), krill_core = ifelse(krill >= krill_75p_thresh,1,0), Area=1) %>%
  filter(krill_core == 1)

bwkr_core_overlap <- bwkr_all %>%
  mutate(blwh_core = ifelse(blwh >= 0.44,1,0), krill_core = ifelse(krill >= krill_75p_thresh,1,0), Area=1) %>%
  filter(blwh_core == 1 & krill_core == 1)

ggplot(bwkr_core_overlap, aes(x=month)) + 
  geom_histogram(binwidth=1) # Essentially only July and August

bwkr_core_overlap <- filter(bwkr_core_overlap, month == "7" | month == "8" )
dates_overlap <- unique(bwkr_core_overlap$date)

# Subset raster stack to jul/aug
dates_overlap_index <- data.frame(matrix(NA,1,62)) # Indexes of dates with overlap present
for (i in 1:62) {
  match <- grep(dates_overlap[i], dates, value = FALSE)
  if (length(match)==0) next
  dates_overlap_index[i] <- match
}

krill_julaug_raster <- subset(krill_all_raster,dates_overlap_index)
blwh_julaug_raster <- subset(blwh_all_raster,dates_overlap_index)

# Function to create categorical raster with krill, blwh and blwh+krill
library(terra)
overlap_binary_fn <- function(raster1,raster2,threshold1, threshold2) {
  
  # Using terra Spatraster to compute functions
  krill_binary <- as(raster1, "SpatRaster")
  krill_binary <- app(krill_binary, fun=function(x){ x[x < threshold1] <- NA; return(x)})
  krill_binary <- app(krill_binary, fun=function(x){ x[x >= threshold1] <- 1; return(x)})
  
  blwh_binary <- as(raster2, "SpatRaster")
  blwh_binary <- app(blwh_binary, fun=function(x){ x[x < threshold2] <- NA; return(x)})
  blwh_binary <- app(blwh_binary, fun=function(x){ x[x >= threshold2] <- 2; return(x)})
  
  bwkr_binary <- mask(krill_binary,blwh_binary)
  bwkr_binary <- app(bwkr_binary, fun=function(x){ x[x == 1] <- 3; return(x)})
  
  # Returning to raster
  krill_binary <- raster(krill_binary)
  blwh_binary <- raster(blwh_binary)
  bwkr_binary <- raster(bwkr_binary)
  
  output <- cover(bwkr_binary,krill_binary)
  output <- cover(output,blwh_binary)
  
  # Create categorical raster with 3 (or 2) classes
  output <- ratify(output)
  ratified <- levels(output)[[1]]
  if (length(ratified$ID) == 2) {
    ratified$CS <- c('blwh', 'krill + blwh')
  } else {
    ratified$CS <- c('krill', 'blwh', 'krill + blwh')
  }
  
  levels(output) <- ratified
  output <- output
  
}

# Create rasters and pngs for Jul/Aug monthly means and individual months

bwkr_julmean_raster <- overlap_binary_fn(krill_jul_mean,blwh_jul_mean,threshold1=krill_75p_thresh, threshold2=0.44)
bwkr_augmean_raster <- overlap_binary_fn(krill_aug_mean,blwh_aug_mean,threshold1=krill_75p_thresh, threshold2=0.44)

getwd()
png(file = paste0("core_overlap_julmean.png"), width=1800, height=1800, res=300)
myTheme <- rasterTheme(panel.background = list(col = 'lightblue'))
plt <- levelplot(bwkr_julmean_raster,col.regions=c('palegreen', 'midnightblue', 'indianred1'),
                 par.settings = myTheme, main = "July Mean") + 
  latticeExtra::layer(sp.polygons(world_land, fill="gray80"))
plot(plt) 
dev.off()

getwd()
png(file = paste0("core_overlap_augmean.png"), width=1800, height=1800, res=300)
myTheme <- rasterTheme(panel.background = list(col = 'lightblue'))
plt <- levelplot(bwkr_augmean_raster,col.regions=c('palegreen', 'midnightblue', 'indianred1'),
                 par.settings = myTheme, main = "August Mean") + 
  latticeExtra::layer(sp.polygons(world_land, fill="gray80"))
plot(plt) 
dev.off()

bwkr_overlap_raster <- stack()
for (i in 1:62) {
  a <- overlap_binary_fn(krill_julaug_raster[[i]],blwh_julaug_raster[[i]],threshold1=krill_75p_thresh, threshold2=0.44)
  bwkr_overlap_raster <- stack(bwkr_overlap_raster,a)
}

for (i in seq_along(dates_overlap)) {
  getwd()
  png(file = paste0("core_overlap_",dates_overlap[i],".png"), width=1800, height=1800, res=300)
  myTheme <- rasterTheme(panel.background = list(col = 'lightblue'))
  plt <- levelplot(bwkr_overlap_raster[[i]],col.regions=c('palegreen', 'midnightblue', 'indianred1'),
                   par.settings = myTheme, main = paste0(dates_overlap[i])) + 
    latticeExtra::layer(sp.polygons(world_land, fill="gray80"))
  plot(plt) 
  dev.off()
}

# Counting pixels of overlap

bwkr_julaug_count <- bwkr_core_overlap %>%
  group_by(year,month) %>%
  summarize(count = n())
bwkr_julaug_count$month[bwkr_julaug_count$month == 7] <- "July"
bwkr_julaug_count$month[bwkr_julaug_count$month == 8] <- "August"

ggplot(bwkr_julaug_count,aes(x=year,y=count,color=month)) +
  geom_line() +
  ggtitle("Pixels of Krill and Blue Whale Core Habitat Overlap")

