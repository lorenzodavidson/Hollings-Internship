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

blwh_all_raster <- raster(paste0(processed_path,"blwh_all_raster.grd"))
krill_all_raster <- raster(paste0(processed_path,"krill_all_raster.grd"))

blwh_apr_raster <- raster(paste0(processed_path,"blwh_apr_raster.grd"))
blwh_may_raster <- raster(paste0(processed_path,"blwh_may_raster.grd"))
blwh_jun_raster <- raster(paste0(processed_path,"blwh_jun_raster.grd"))
blwh_jul_raster <- raster(paste0(processed_path,"blwh_jul_raster.grd"))
blwh_aug_raster <- raster(paste0(processed_path,"blwh_aug_raster.grd"))

krill_apr_raster <- raster(paste0(processed_path,"krill_apr_raster.grd"))
krill_may_raster <- raster(paste0(processed_path,"krill_may_raster.grd"))
krill_jun_raster <- raster(paste0(processed_path,"krill_jun_raster.grd"))
krill_jul_raster <- raster(paste0(processed_path,"krill_jul_raster.grd"))
krill_aug_raster <- raster(paste0(processed_path,"krill_aug_raster.grd"))

# Freq plot by month -----------------------------------------------------

months <- c('Apr', 'May', 'Jun', 'Jul', 'Aug')
bwkr_bymonth <- data.frame(group_by(bwkr_all,monthnum))
bwkr_bymonth$monthnum <- as.character(bwkr_bymonth$monthnum)

ggplot(bwkr_bymonth, aes(blwh,color=monthnum)) +
  geom_freqpoly(binwidth = 0.01,linewidth=1) +
  xlim(0, 1) +
  ylim(0,20000) +
  ggtitle("Monthly Blue Whale HS Frequency") +
  labs(x = "Blwh Habitat Suitability",
       y = "Count",
       color = "Months")

ggplot(bwkr_bymonth, aes(krill,color=monthnum)) +
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

myTheme <- rasterTheme(region = brewer.pal(9,"YlOrRd"),panel.background = list(col = 'grey95'))

# Blue Whale
plt <- levelplot(blwh_all_mean,main = "Mean Monthly Blue Whale HS",par.settings = myTheme,at = blwh_breaks)
plt + latticeExtra::layer(sp.polygons(world_land, fill="gray80"))

# Krill
plt <- levelplot(krill_all_mean,main = "Mean Monthly Krill CPUE",par.settings = myTheme,at = krill_breaks)
plt + latticeExtra::layer(sp.polygons(world_land, fill="gray80"))


# June, July and August blue whale and krill comparison --------------------

blwh_678_mean <- subset(blwh_all_mean,3:5)
krill_678_mean <- subset(krill_all_mean,3:5)

blwh_breaks <- unname(quantile(filter(bwkr_bymonth, monthnum == "6" | monthnum == "7" | monthnum == "8")$blwh, 
                               probs = c(0.35,0.45,0.55,0.65,0.75,0.85,0.95,1)))
blwh_breaks <- c(0, blwh_breaks)

krill_breaks <- unname(quantile(filter(bwkr_bymonth, monthnum == "6" | monthnum == "7" | monthnum == "8")$krill, 
                                probs = c(0.35,0.45,0.55,0.65,0.75,0.85,0.95,1)))
krill_breaks <- c(0,krill_breaks)

library(gridExtra)
p1 <- levelplot(blwh_678_mean,main = "Mean Monthly Blue Whale HS",par.settings = myTheme,at = blwh_breaks) + 
      latticeExtra::layer(sp.polygons(world_land, fill="gray80"))
p2 <- levelplot(krill_678_mean,main = "Mean Monthly Krill CPUE",par.settings = myTheme,at = krill_breaks) + 
      latticeExtra::layer(sp.polygons(world_land, fill="gray80"))
grid.arrange(p1, p2, ncol=1)

# Investigating persistence of hot spots ----------------------------------

# Function to filter for months and threshold, group by xy, and calculate % of years
thresh_fn <- function(df, time, thresh) {
  df_count <- df %>%
    filter(monthnum == time) %>% # time = month of interest ("6","7" or "8")
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

blwh_percentiles <- unname(quantile(filter(bwkr_bymonth, monthnum == "6" | monthnum == "7" | monthnum == "8")$blwh, 
                probs = c(0.35,0.45,0.55,0.65,0.75,0.85,0.95,1)))

blwh_jun_count <- stack()
blwh_jul_count <- stack()
blwh_aug_count <- stack()
for (i in 3:7) {
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
names(blwh_jun_count) <- paste0(thresholds)
names(blwh_jul_count) <- paste0(thresholds)
names(blwh_aug_count) <- paste0(thresholds)

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
      group_by(year,monthnum) %>% # organize by date
      summarise(AO=area_overlapfn(prey=krill_core,pred=blwh_core,area=Area)) %>% 
      pivot_longer(cols=c(AO),names_to="Metric_Name",values_to="Overlap_Metric")
    boxplot(Overlap_Metric~monthnum,
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
      group_by(year,monthnum) %>% # organize by date
      summarise(RO=range_overlapfn(prey=krill_core,pred=blwh_core,area=Area)) %>% 
      pivot_longer(cols=c(RO),names_to="Metric_Name",values_to="Overlap_Metric")
    boxplot(Overlap_Metric~monthnum,
            data=RO_plot,
            xlab="Month",
            ylab="Range Overlap",
            # main=paste0(thresholds[[j]],"-",thresholds[[i]]),
            ylim = c(0, 1))
  }
}

dev.off()

# Metrics for thresh: blwh 0.28 + krill 75% -------------------------------

# Matching krill threshold to % blwh above 0.28 threshold and calculating metrics
blwh_028 <- filter(bwkr_all,bwkr_all$blwh>0.28) # 23% of data
krill_75p_thresh <- unname(quantile(bwkr_all$krill, 0.75))
krill_75p <- filter(bwkr_all,bwkr_all$krill>krill_75p_thresh) 

bwkr_metrics <- bwkr_all %>% # make new dataframe called dailyzoom2015 which is a copy of sdm2015 and do the following below
  mutate(blwh_core = ifelse(blwh >= 0.28,1,0), krill_core = ifelse(krill >= krill_75p_thresh,1,0), Area=1) %>% # create a binary data column for range and area overlap of each species
  group_by(date,year,month) %>% # organize by date
  summarise(AO=area_overlapfn(prey=krill_core,pred=blwh_core,area=Area),
            RO=range_overlapfn(prey=krill_core,pred=blwh_core,area=Area),
            Schoener=schoeners_overlapfn(prey=krill,pred=blwh),
            Bhatty=bhatta_coeffn(prey=krill,pred=blwh)) %>% 
  pivot_longer(cols=c(AO,RO,Schoener,Bhatty),names_to="Metric_Name",values_to="Overlap_Metric") # change the dataframe into long format for easier plotting

AO_all <- filter(bwkr_metrics,Metric_Name=="AO")
RO_all <- filter(bwkr_metrics,Metric_Name=="RO")
Schoener_all <- filter(bwkr_metrics,Metric_Name=="Schoener")
Bhatty_all <- filter(bwkr_metrics,Metric_Name=="Bhatty")

# AO Plots
ggplot(AO_all,aes(x=date,y=Overlap_Metric,group = 1)) +
  geom_point() + 
  geom_line() +
  ggtitle("AO 1990-2020") +
  ylab("Area Overlap")
ggplot(AO_all,aes(x=month,y=Overlap_Metric,group=month)) +
  geom_boxplot() +
  ggtitle("AO grouped by month") +
  xlab("Month (March-July)") + 
  ylab("Area Overlap")

# RO Plots
ggplot(RO_all,aes(x=date,y=Overlap_Metric,group = 1)) +
  geom_point() + 
  geom_line() +
  ggtitle("RO 1990-2020") +
  ylab("Range Overlap")
ggplot(RO_all,aes(x=month,y=Overlap_Metric,group=month)) +
  geom_boxplot() +
  ggtitle("RO grouped by month") +
  xlab("Month (March-July)") + 
  ylab("Range Overlap")

# Schoener Plots
ggplot(Schoener_all,aes(x=date,y=Overlap_Metric,group = 1)) +
  geom_point() + 
  geom_line() +
  ggtitle("Schoener's D 1990-2020") +
  ylab("Schoener's D")
ggplot(Schoener_all,aes(x=month,y=Overlap_Metric,group=month)) +
  geom_boxplot() +
  ggtitle("Schoener's D grouped by month") +
  xlab("Month (March-July)") + 
  ylab("Schoener's D")

# Bhattacharyya Plots
ggplot(Bhatty_all,aes(x=date,y=Overlap_Metric,group = 1)) +
  geom_point() + 
  geom_line() +
  ggtitle("Bhattacharyya's Coefficient 1990-2020") +
  ylab("Bhattacharyya's Coefficient")
ggplot(Bhatty_all,aes(x=month,y=Overlap_Metric,group=month)) +
  geom_boxplot() +
  ggtitle("Bhattacharyya's Coefficient grouped by month") +
  xlab("Month (March-July)") + 
  ylab("Bhattacharyya's Coefficient")