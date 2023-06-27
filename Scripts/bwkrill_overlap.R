#Libraries ---------------------------------------------------------------
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
    rename("x" = "X1","y" = "X2","krill" = "X3")
  #krill_raster <- rasterFromXYZ(krill_rasterprep,crs=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  #krill_monthly_raster[[t]] <- krill_raster
}

# Creating bwkr_monthly_xyz dataset
bwkr_monthly_xyz <- list()
for (t in seq_along(time)) {
  bwkr_monthly_xyz[[t]] <- inner_join(blwh_monthly_xyz[[t]],krill_monthly_xyz[[t]]) %>%
    rename("blwh" = "layer")
}

# Creating bwkr_monthly_raster dataset (layer 1 = blwh, layer 2 = krill)
ext <- extent(-127.45, -115.55, 30.05, 47.95) # Extent of bwkr raster with edges

bwkr_monthly_raster <- list()
for (t in seq_along(time)) {
  bwkr_rasterprep <- bwkr_monthly_xyz[[t]][,-5:-7]
  blwh_raster <- extend(rasterFromXYZ(bwkr_rasterprep[,-4]),ext)
  krill_raster <- extend(rasterFromXYZ(bwkr_rasterprep[,-3]),ext)
  bwkr_monthly_raster[[t]] <- stack(blwh_raster,krill_raster)                  
}

# Creating bwkr_all 
months <- c("03","04","05","06","07")
years <- c(1990:2020)
dates <- vector("character", 155)
year_month <- data.frame(matrix(NA,155,2)) %>%
  rename("year"= "X1") %>%
  rename("month"= "X2")

for (i in seq_along(years)) {
  for (j in seq_along(months)) {
    x <- i*5 - 5 + j
    if (i == 1 & j == 1) {
      dates[[x]] <- as.character(paste0(years[[i]],"-",months[[j]]))
      year_month[[x,1]] <- as.character(paste0(years[[i]]))
      year_month[[x,2]] <- as.character(paste0(months[[j]]))
    } else {
      dates[[x]] <- rbind(as.character(paste0(years[[i]],"-",months[[j]])))
      year_month[[x,1]] <- as.character(paste0(years[[i]]))
      year_month[[x,2]] <- as.character(paste0(months[[j]]))
    }
  }
}

bwkr_all <- data.frame(bwkr_monthly[[1]])
bwkr_all$date <- dates[[1]]
bwkr_all$year <- year_month$year[[1]]
bwkr_all$month <- year_month$month[[1]]
for (i in 2:155) {
  bwkr_add <- bwkr_monthly[[i]]
  bwkr_add$date <- dates[[i]]
  bwkr_add$year <- year_month$year[[i]]
  bwkr_add$month <- year_month$month[[i]]
  bwkr_all <- rbind(bwkr_all, bwkr_add)
}

# Save all 7 datasets as csv or raster files
fpath_bwkr <- "~/Dropbox/blwh_krill/bwkr_monthly"
fpath_blwh <- "~/Dropbox/blwh_krill/blwh_monthly"
fpath_krill <- "~/Dropbox/blwh_krill/krill_monthly"
for (i in seq_along(years)) {
  for (j in seq_along(months)) {
    x <- i*5 - 5 + j
    bwkr_name <- paste0("bwkr_", years[[i]], "_", months[[j]])
    write.csv(bwkr_monthly[[x]],paste0(fpath_bwkr,"/",bwkr_name,".csv"),row.names = FALSE)
    writeRaster(bwkr_monthly_raster[[x]],paste0(fpath_bwkr,"/",bwkr_name,".grd"),format = "raster")
    
    blwh_name <- paste0("blwh_", years[[i]], "_", months[[j]])
    write.csv(blwh_monthly_xyz[[x]],paste0(fpath_blwh,"/",blwh_name,".csv"),row.names = FALSE)
    writeRaster(blwh_monthly_raster[[x]],paste0(fpath_blwh,"/",blwh_name,".grd"),format = "raster")

    krill_name <- paste0("krill_", years[[i]], "_", months[[j]])
    write.csv(krill_monthly_xyz[[x]],paste0(fpath_krill,"/",krill_name,".csv"),row.names = FALSE)
    writeRaster(krill_monthly_raster[[x]],paste0(fpath_krill,"/",krill_name,".grd"),format = "raster")
  }
}

write.csv(bwkr_all,paste0(fpath_bwkr,"/bwkr_all.csv"),row.names = FALSE) # bwkr_all

# Opening whale and krill data --------------------------------------------
fpath_bwkr <- "~/Dropbox/blwh_krill/bwkr_monthly"
fpath_blwh <- "~/Dropbox/blwh_krill/blwh_monthly"
fpath_krill <- "~/Dropbox/blwh_krill/krill_monthly"

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

bwkr_monthly_raster <- list()
for(i in dfnames_bwkr[stringr::str_detect(dfnames_bwkr, ".grd")]){
  tt <- raster(paste0(fpath_bwkr,"/",i))
  bwkr_monthly_raster[[i]] <- tt
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

bwkr_monthly_xyz <- list()
for(i in dfnames_bwkr[stringr::str_detect(dfnames_bwkr, ".csv")]){
  tt <- read.csv(paste0(fpath_bwkr,"/",i))
  bwkr_monthly_xyz[[i]] <- tt
}

bwkr_all <- read.csv(paste0(fpath_bwkr,"/bwkr_all.csv"))

# Freq plot by month -----------------------------------------------------
bwkr_march <- filter(bwkr_all, month == 3)
bwkr_april <- filter(bwkr_all, month == 4)
bwkr_may <- filter(bwkr_all, month == 5)
bwkr_june <- filter(bwkr_all, month == 6)
bwkr_july <- filter(bwkr_all, month == 7)

bwkr_bymonth <- data.frame(group_by(bwkr_all,month))
bwkr_bymonth$month <- as.character(bwkr_bymonth$month) # Convert to character for grouping in ggplot

ggplot(bwkr_bymonth, aes(blwh,color=month)) +
  geom_freqpoly(binwidth = 0.01,linewidth=1) +
  xlim(0, 1) +
  ylim(0,20000) +
  ggtitle("Monthly Blue Whale HS Frequency")

ggplot(bwkr_bymonth, aes(krill,color=month)) +
  geom_freqpoly(binwidth = 0.02,linewidth=1) +
  xlim(range(bwkr_bymonth$krill)[[1]], range(bwkr_bymonth$krill)[[2]]) +
  # ylim(0,20000) +
  ggtitle("Monthly Krill HS Frequency")

# Overlap Metrics Code ----------------------------------------------------

# Center of Gravity


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
# Initial investigation of threshold values --------------------------------

# Area and Range Overlap vs threshold value
krill_threshold_t3 <- unname(quantile(bwkr_monthly[[3]]$krill, probs = seq(.01, 1, by = .01)))
krill_threshold_t5 <- unname(quantile(bwkr_monthly[[5]]$krill, probs = seq(.01, 1, by = .01)))
blwh_threshold <- c(0.15,0.25,0.35,0.45)

AO_bwkr_t3 <- data.frame(matrix(NA,length(krill_threshold_t3),length(blwh_threshold)))
colnames(AO_bwkr_t3) <- c("blwh_15","blwh_25","blwh_35","blwh_45")
RO_bwkr_t3 <- AO_bwkr_t3
for (j in 1:4) {
  for (i in seq_along(krill_threshold_t3)) {
    AO_bwkr_t3[[i,j]] <- bwkr_monthly[[3]] %>%
      mutate(blwh_core = ifelse(blwh >= blwh_threshold[[j]],1,0), krill_core = ifelse(krill >= krill_threshold_t3[[i]],1,0), Area=1) %>%
      summarise(AO=area_overlapfn(prey=krill_core,pred=blwh_core,area=Area)) %>%
      unlist()
    RO_bwkr_t3[[i,j]] <- bwkr_monthly[[3]] %>%
      mutate(blwh_core = ifelse(blwh >= blwh_threshold[[j]],1,0), krill_core = ifelse(krill >= krill_threshold_t3[[i]],1,0), Area=1) %>%
      summarise(RO=range_overlapfn(prey=krill_core,pred=blwh_core,area=Area)) %>%
      unlist()
  }
}

AO_bwkr_t5 <- data.frame(matrix(NA,length(krill_threshold_t5),length(blwh_threshold)))
colnames(AO_bwkr_t5) <- c("blwh_15","blwh_25","blwh_35","blwh_45")
RO_bwkr_t5 <- AO_bwkr_t5
for (j in 1:4) {
  for (i in seq_along(krill_threshold_t5)) {
    AO_bwkr_t5[[i,j]] <- bwkr_monthly[[5]] %>%
      mutate(blwh_core = ifelse(blwh >= blwh_threshold[[j]],1,0), krill_core = ifelse(krill >= krill_threshold_t5[[i]],1,0), Area=1) %>%
      summarise(AO=area_overlapfn(prey=krill_core,pred=blwh_core,area=Area)) %>%
      unlist()
    RO_bwkr_t5[[i,j]] <- bwkr_monthly[[5]] %>%
      mutate(blwh_core = ifelse(blwh >= blwh_threshold[[j]],1,0), krill_core = ifelse(krill >= krill_threshold_t5[[i]],1,0), Area=1) %>%
      summarise(RO=range_overlapfn(prey=krill_core,pred=blwh_core,area=Area)) %>%
      unlist()
  }
}

# May 1990
plot(rep(krill_threshold_t3[75],4),AO_bwkr_t3[75,],main="May 1990 AO vs Krill Threshold",
     xlab="Krill Threshold", ylab="Area Overlap", ylim=c(0,1),xlim=c(2,11))
lines(krill_threshold_t3,AO_bwkr_t3$blwh_15,col="red")
lines(krill_threshold_t3,AO_bwkr_t3$blwh_25,col="orange")
lines(krill_threshold_t3,AO_bwkr_t3$blwh_35,col="green")
lines(krill_threshold_t3,AO_bwkr_t3$blwh_45,col="blue")
legend(7.5, 0.95, legend=c("75th Percentile", "bw_tresh: 0.15", "bw_tresh: 0.25","bw_tresh: 0.35","bw_tresh: 0.45"),
       col=c("black","red","orange","green","blue"),lty=1:2, cex=0.8)

plot(rep(krill_threshold_t3[75],4),RO_bwkr_t3[75,],main="May 1990 RO vs Krill Threshold",
     xlab="Krill Threshold", ylab="Range Overlap", ylim=c(0,1),xlim=c(2,11))
lines(krill_threshold_t3,RO_bwkr_t3$blwh_15,col="red")
lines(krill_threshold_t3,RO_bwkr_t3$blwh_25,col="orange")
lines(krill_threshold_t3,RO_bwkr_t3$blwh_35,col="green")
lines(krill_threshold_t3,RO_bwkr_t3$blwh_45,col="blue")
legend(7.5, 0.95, legend=c("75th Percentile", "bw_tresh: 0.15", "bw_tresh: 0.25","bw_tresh: 0.35","bw_tresh: 0.45"),
       col=c("black","red","orange","green","blue"),lty=1:2, cex=0.8)

# July 1990
plot(rep(krill_threshold_t5[75],4),AO_bwkr_t5[75,],main="July 1990 AO vs Krill Threshold",
     xlab="Krill Threshold", ylab="Area Overlap", ylim=c(0,1),xlim=c(2,11))
lines(krill_threshold_t5,AO_bwkr_t5$blwh_15,col="red")
lines(krill_threshold_t5,AO_bwkr_t5$blwh_25,col="orange")
lines(krill_threshold_t5,AO_bwkr_t5$blwh_35,col="green")
lines(krill_threshold_t5,AO_bwkr_t5$blwh_45,col="blue")
legend(7.5, 0.95, legend=c("75th Percentile", "bw_tresh: 0.15", "bw_tresh: 0.25","bw_tresh: 0.35","bw_tresh: 0.45"),
       col=c("black","red","orange","green","blue"),lty=1:2, cex=0.8)

plot(rep(krill_threshold_t5[75],4),RO_bwkr_t5[75,],main="July 1990 RO vs Krill Threshold",
     xlab="Krill Threshold", ylab="Range Overlap", ylim=c(0,1),xlim=c(2,11))
lines(krill_threshold_t5,RO_bwkr_t5$blwh_15,col="red")
lines(krill_threshold_t5,RO_bwkr_t5$blwh_25,col="orange")
lines(krill_threshold_t5,RO_bwkr_t5$blwh_35,col="green")
lines(krill_threshold_t5,RO_bwkr_t5$blwh_45,col="blue")
legend(2, 0.4, legend=c("75th Percentile", "bw_tresh: 0.15", "bw_tresh: 0.25","bw_tresh: 0.35","bw_tresh: 0.45"),
       col=c("black","red","orange","green","blue"),lty=1:2, cex=0.8)


# Overlap Metrics for 9 different threshold combinations -------------------

krill_85p_thresh <- unname(quantile(bwkr_all$krill, 0.85))
krill_75p_thresh <- unname(quantile(bwkr_all$krill, 0.75))
krill_65p_thresh <- unname(quantile(bwkr_all$krill, 0.65))
blwh_85p_thresh <- unname(quantile(bwkr_all$blwh, 0.85))
blwh_75p_thresh <- unname(quantile(bwkr_all$blwh, 0.75))
blwh_65p_thresh <- unname(quantile(bwkr_all$blwh, 0.65))
krill_thresh <- c(krill_65p_thresh,krill_75p_thresh,krill_85p_thresh)
blwh_thresh <- c(blwh_65p_thresh,blwh_75p_thresh,blwh_85p_thresh)

bwkr_metrics_9combo <- bwkr_all
bwkr_metrics_9combo <- bwkr_metrics_9combo %>% 
  mutate(blwh_core_65 = ifelse(blwh >= blwh_tresh[[1]],1,0), 
         blwh_core_75 = ifelse(blwh >= blwh_tresh[[2]],1,0),
         blwh_core_85 = ifelse(blwh >= blwh_tresh[[3]],1,0),
         krill_core_65 = ifelse(krill >= krill_tresh[[1]],1,0),
         krill_core_75 = ifelse(krill >= krill_tresh[[2]],1,0),
         krill_core_85 = ifelse(krill >= krill_tresh[[3]],1,0), Area=1)

blwh_cores <- c("blwh_core_65","blwh_core_75","blwh_core_85")
krill_cores <- c("krill_core_65","krill_core_75","krill_core_85")
thresholds <- c("65","75","85")

# AO Plots
par(mfrow=c(3,3))
for (i in 1:3) { # blwh_tresh
  for (j in 1:3) { # krill_tresh
    x <- 7 + i # Column of blwh_core desired
    y <- 10 + j # Column of krill_core desired
    AO_plot <- bwkr_metrics_9combo[c(1:7,x,y,14)] %>%
      rename("blwh_core"=paste0(blwh_cores[[i]]),"krill_core"=paste0(krill_cores[[j]])) %>%
      group_by(date,year,month) %>% # organize by date
      summarise(AO=area_overlapfn(prey=krill_core,pred=blwh_core,area=Area)) %>% 
      pivot_longer(cols=c(AO),names_to="Metric_Name",values_to="Overlap_Metric")
    boxplot(Overlap_Metric~month,
            data=AO_plot,
            xlab="Month Number",
            ylab="Area Overlap",
            # main=paste0(thresholds[[j]],"-",thresholds[[i]]),
            ylim = c(0, 0.30))
  }
}

# RO Plots
par(mfrow=c(3,3))
for (i in 1:3) { # blwh_tresh
  for (j in 1:3) { # krill_tresh
    x <- 7 + i # Column of blwh_core desired
    y <- 10 + j # Column of krill_core desired
    
    RO_plot <- bwkr_metrics_9combo[c(1:7,x,y,14)] %>%
      rename("blwh_core"=paste0(blwh_cores[[i]]),"krill_core"=paste0(krill_cores[[j]])) %>%
      group_by(year,month) %>% # organize by date
      summarise(RO=range_overlapfn(prey=krill_core,pred=blwh_core,area=Area)) %>% 
      pivot_longer(cols=c(RO),names_to="Metric_Name",values_to="Overlap_Metric")
    boxplot(Overlap_Metric~month,
            data=RO_plot,
            xlab="Month Number",
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

# Creating pngs/rasters of area overlap and high blwh HS -------------------

# High blwh HS and AO plots + png
bwkr_028filter <- data.frame(bwkr_all) %>%
  mutate(Area=1) %>% 
  filter(blwh >= 0.28)

AO_028filter <- data.frame(bwkr_all) %>% 
  mutate(Area=1) %>% 
  filter(blwh >= 0.28 & krill >= krill_75p_thresh)

ggplot(bwkr_028filter, aes(x=month)) + 
  geom_histogram(binwidth=1) # Majority are in June/July but ~10% in March/April/May
ggplot(AO_028filter, aes(x=month)) + 
  geom_histogram(binwidth=1) # Almost all are in June/July

bwkr_junjul_filter <- filter(bwkr_028filter, month == "6" | month == "7" )
AO_junjul_filter <- filter(AO_028filter, month == "6" | month == "7" )
dates_filtered <- unique(AO_junjul_filter$date)

# Fix rasters to be same size as original
ext <- extent(-133.95, -115.55, 30.05, 47.95)

bwkr_junjul_raster <- list()
AO_junjul_raster <- list()
for (i in 1:length(dates_filtered)) {
  
  bwkr_rasterprep <- filter(bwkr_junjul_filter, date == dates_filtered[[i]]) 
  bwkr_rasterprep <- rasterFromXYZ(bwkr_rasterprep[,-3:-7]) 
  bwkr_junjul_raster[[i]] <- extend(bwkr_rasterprep, ext)
  
  AO_rasterprep <- filter(AO_junjul_filter, date == dates_filtered[[i]])
  AO_rasterprep <- rasterFromXYZ(AO_rasterprep[,-3:-7]) 
  AO_junjul_raster[[i]] <- extend(AO_rasterprep, ext) # For some reason [[5]] becomes 178x185
}

plot(bwkr_junjul_raster[[48]])
plot(AO_junjul_raster[[48]])

# Save bwkr_all, bwkr_junjul_filter, and AO_junjul_filter as pngs and rasters

fpath_pngs <- "~/Dropbox/blwh_krill/AO_HS_junjul/pngs/"
fpath_rasters <- "~/Dropbox/blwh_krill/AO_HS_junjul/rasters/"
for (i in seq_along(dates_filtered)) {
  bwkr_file <- bwkr_junjul_raster[[i]]
  AO_file <- AO_junjul_raster[[i]]
  name <- paste0("blwhHS_AO_", dates_filtered[[i]])

  png(file = paste0(fpath_pngs,name,".png"))
  
  plot(bwkr_file, xlab=as.character(dates_filtered[[i]]), main="Blue Whale HS > 0.28 and Blue Whale-Krill Overlap (Red)", col="#fdbb84")
  plot(AO_file, add=T, breaks=c(0.3,1), legend=F, col="#e34a33")
  dev.off()
}

# Blwh HS for June and July -------------------------------------
fpath_bwkr <- "~/Dropbox/blwh_krill/bwkr_monthly"
dfnames_bwkr <- list.files(fpath_bwkr)
bwkr_all <- read.csv(paste0(fpath_bwkr,"/bwkr_all.csv"))
bwkr_junjul <- data.frame(bwkr_all) %>%
  filter(month == 6 | month == 7) %>%
  group_by(month)
bwkr_junjul$month <- as.character(bwkr_junjul$month) # For grouping purposes

ggplot(bwkr_junjul, aes(blwh,color=month)) +
  geom_freqpoly(binwidth = 0.01,linewidth=1) +
  xlim(0, 1) +
  ylim(0,20000) +
  ggtitle("Monthly Blue Whale HS Frequency")

# Visualizing Mean
bwkr_monthly_raster <- list()
for(i in dfnames_bwkr[stringr::str_detect(dfnames_bwkr, ".grd")]){
  tt <- raster(paste0(fpath_bwkr,"/",i))
  bwkr_monthly_raster[[i]] <- tt
}

bwkr_jun_stacked <- stack()
bwkr_jul_stacked <- stack()
bwkr_junjul_stacked <- stack()
for (i in 1:31) {
  jun <- 5*i-1
  jul <- 5*i
  jun_raster <- bwkr_monthly_raster[[jun]]
  jul_raster <- bwkr_monthly_raster[[jul]]
  bwkr_jun_stacked <- stack(bwkr_jun_stacked,jun_raster)
  bwkr_jul_stacked <- stack(bwkr_jul_stacked,jul_raster)
  bwkr_junjul_stacked <- stack(bwkr_junjul_stacked,jun_raster,jul_raster)
}

bwkr_jun_mean <- mean(bwkr_jun_stacked)
bwkr_jul_mean <- mean(bwkr_jul_stacked)
bwkr_junjul_mean <- mean(bwkr_junjul_stacked)

cols <- colorRampPalette(brewer.pal(9,"YlOrRd"))
levelplot(bwkr_jun_mean,main = "Mean June Blue Whale HS",col.regions=cols,at =  c(0, 0.25, 0.34, 0.44, 0.54, 0.64, 0.74, 1))
levelplot(bwkr_jul_mean,main = "Mean July Blue Whale HS",col.regions=cols,at =  c(0, 0.25, 0.34, 0.44, 0.54, 0.64, 0.74, 1))
levelplot(bwkr_junjul_mean,main = "Mean June/July Blue Whale HS",col.regions=cols,at =  c(0, 0.25, 0.45, 0.55, 0.65, 0.75, 1))

# Investigating persistence of hotspots -----------------------------------
fpath_bwkr <- "~/Dropbox/blwh_krill/bwkr_monthly"
bwkr_all <- read.csv(paste0(fpath_bwkr,"/bwkr_all.csv"))
bwkr_junjul <- data.frame(bwkr_all) %>%
  filter(month == 6 | month == 7) %>%
  group_by(month)

# Function to filter for months and threshold, group by xy, and count pixels
thresh_fn <- function(df, time, thresh) {
  df_count <- df %>%
    filter(month == time) %>% # time = month of interest (6 and/or 7)
    filter(blwh >= thresh) %>%
    group_by(x,y) %>%
    mutate(n=n()) %>%
    arrange(desc(n))
  output <- df_count[!duplicated(df_count[,c("x","y")]),]
}

# Creating rasters with proper extension
ext <- extent(-127.45, -115.55, 30.05, 47.95) # Extent of bwkr raster with edges
raster_fn <- function(df, ext) {
  extend(rasterFromXYZ(df[c(-3:-7)]),ext)
}

junjul_percentiles <- unname(quantile(bwkr_junjul$blwh, probs = c(0.55,0.65,0.75,0.85,0.95)))
jun_percentiles <- quantile(filter(bwkr_junjul,month ==6)$blwh, probs = c(0.65,0.75,0.85,0.95))
jul_percentiles <- quantile(filter(bwkr_junjul,month ==7)$blwh, probs = c(0.65,0.75,0.85,0.95))

# Creating rasters for jun jul at different thresholds
thresh <- 0.35 # 55%
bwkr_035count_jun <- thresh_fn(bwkr_junjul,6,thresh)
bwkr_035count_jun_raster <- raster_fn(bwkr_035count_jun,ext)
bwkr_035count_jul <- thresh_fn(bwkr_junjul,7,thresh)
bwkr_035count_jul_raster <- raster_fn(bwkr_035count_jul,ext)

thresh <- 0.45 # 65%
bwkr_045count_jun <- thresh_fn(bwkr_junjul,6,thresh)
bwkr_045count_jun_raster <- raster_fn(bwkr_045count_jun,ext)
bwkr_045count_jul <- thresh_fn(bwkr_junjul,7,thresh)
bwkr_045count_jul_raster <- raster_fn(bwkr_045count_jul,ext)

thresh <- 0.55 # 75%
bwkr_055count_jun <- thresh_fn(bwkr_junjul,6,thresh)
bwkr_055count_jun_raster <- raster_fn(bwkr_055count_jun,ext)
bwkr_055count_jul <- thresh_fn(bwkr_junjul,7,thresh)
bwkr_055count_jul_raster <- raster_fn(bwkr_055count_jul,ext)

thresh <- 0.65 # 85%
bwkr_065count_jun <- thresh_fn(bwkr_junjul,6,thresh)
bwkr_065count_jun_raster <- raster_fn(bwkr_065count_jun,ext)
bwkr_065count_jul <- thresh_fn(bwkr_junjul,7,thresh)
bwkr_065count_jul_raster <- raster_fn(bwkr_065count_jul,ext)

thresh <- 0.75 # 95%
bwkr_075count_jun <- thresh_fn(bwkr_junjul,6,thresh)
bwkr_075count_jun_raster <- raster_fn(bwkr_075count_jun,ext)
bwkr_075count_jul <- thresh_fn(bwkr_junjul,7,thresh)
bwkr_075count_jul_raster <- raster_fn(bwkr_075count_jul,ext)

jun_stacked_thresh <- stack(bwkr_035count_jun_raster,bwkr_045count_jun_raster,bwkr_055count_jun_raster,
                            bwkr_065count_jun_raster,bwkr_075count_jun_raster)
jul_stacked_thresh <- stack(bwkr_035count_jul_raster,bwkr_045count_jul_raster,bwkr_055count_jul_raster,
                            bwkr_065count_jul_raster,bwkr_075count_jul_raster)

# Add coastline for plotting
fpath_ne <- "~/Dropbox/blwh_krill/ne_10m_land/"
dfname_ne <- "ne_10m_land.shp"
necoast_data <- readOGR(dsn = fpath_ne, 
                        layer = file_path_sans_ext(dfname_ne))

necoast_data2 <- read_sf(paste0(fpath_ne,dfname_ne))
world_outline <- as(st_geometry(necoast_data2), Class="Spatial")

# Pixel count at each threshold

library(colorspace)

thresholds <- c("Thresh 0.35","Thresh 0.45","Thresh 0.55","Thresh 0.65","Thresh 0.75")
jun_stacked_thresh <- setZ(jun_stacked_thresh, thresholds)
jul_stacked_thresh <- setZ(jul_stacked_thresh, thresholds)
names(jun_stacked_thresh) <- thresholds
names(jul_stacked_thresh) <- thresholds

myTheme <- rasterTheme(region = brewer.pal(9,"YlOrRd"),
                       panel.background = list(col = 'grey90'))

plt <- levelplot(jun_stacked_thresh,main = "June Pixel Count above Threshold",par.settings = myTheme)
plt + latticeExtra::layer(sp.lines(world_outline, col="black", lwd=0.5))

plt <- levelplot(jul_stacked_thresh,main = "July Pixel Count above Threshold",par.settings = myTheme)
plt + latticeExtra::layer(sp.lines(world_outline, col="black", lwd=0.5))

# June rasters that plot pixels with n>15 at each threshold
library(terra)

ext <- extent(-127.45, -115.50, 30.05, 47.95)
bwkr_035_jun_binary <- app(rast(bwkr_035count_jun[,-3:-7]), fun=function(x){ x[x <= 15] <- NA; return(x)})
bwkr_035_jun_binary <- extend(bwkr_035_jun_binary,ext)
bwkr_035_jun_binary <- app(bwkr_035_jun_binary, fun=function(x){ x[x > 15] <- 0.35; return(x)})

bwkr_045_jun_binary <- app(rast(bwkr_045count_jun[,-3:-7]), fun=function(x){ x[x <= 15] <- NA; return(x)})
bwkr_045_jun_binary <- extend(bwkr_045_jun_binary,ext)
bwkr_045_jun_binary <- app(bwkr_045_jun_binary, fun=function(x){ x[x > 15] <- 0.45; return(x)})

bwkr_055_jun_binary <- app(rast(bwkr_055count_jun[,-3:-7]), fun=function(x){ x[x <= 15] <- NA; return(x)})
bwkr_055_jun_binary <- extend(bwkr_055_jun_binary,ext)
bwkr_055_jun_binary <- app(bwkr_055_jun_binary, fun=function(x){ x[x > 15] <- 0.55; return(x)})

bwkr_065_jun_binary <- app(rast(bwkr_065count_jun[,-3:-7]), fun=function(x){ x[x <= 15] <- NA; return(x)})
bwkr_065_jun_binary <- extend(bwkr_065_jun_binary,ext)
bwkr_065_jun_binary <- app(bwkr_065_jun_binary, fun=function(x){ x[x > 15] <- 0.65; return(x)})

ext <- extent(-127.45, -115.55, 30.05, 47.95)
bwkr_075_jun_binary <- app(rast(bwkr_075count_jun[,-3:-7]), fun=function(x){ x[x <= 15] <- NA; return(x)})
bwkr_075_jun_binary <- extend(bwkr_075_jun_binary,ext)
bwkr_075_jun_binary <- app(bwkr_075_jun_binary, fun=function(x){ x[x > 15] <- 0.75; return(x)})

bwkr_jun_binary <- cover(bwkr_075_jun_binary,bwkr_065_jun_binary)
bwkr_jun_binary <- cover(bwkr_jun_binary,bwkr_055_jun_binary)
bwkr_jun_binary <- cover(bwkr_jun_binary,bwkr_045_jun_binary)
bwkr_jun_binary <- cover(bwkr_jun_binary,bwkr_035_jun_binary)
bwkr_jun_binary <- raster(bwkr_jun_binary)

plt <- levelplot(bwkr_jun_binary,main = "Persistently High Blue Whale HS Pixels (June)",
                 par.settings = myTheme, at =  c(0, 0.25, 0.44, 0.54, 0.64, 0.74, 1))
plt + latticeExtra::layer(sp.lines(world_outline, col="black", lwd=0.5))

# July
ext <- extent(-127.45, -115.50, 30.05, 47.95)
bwkr_035_jul_binary <- app(rast(bwkr_035count_jul[,-3:-7]), fun=function(x){ x[x <= 15] <- NA; return(x)})
bwkr_035_jul_binary <- extend(bwkr_035_jul_binary,ext)
bwkr_035_jul_binary <- app(bwkr_035_jul_binary, fun=function(x){ x[x > 15] <- 0.35; return(x)})

bwkr_045_jul_binary <- app(rast(bwkr_045count_jul[,-3:-7]), fun=function(x){ x[x <= 15] <- NA; return(x)})
bwkr_045_jul_binary <- extend(bwkr_045_jul_binary,ext)
bwkr_045_jul_binary <- app(bwkr_045_jul_binary, fun=function(x){ x[x > 15] <- 0.45; return(x)})

bwkr_055_jul_binary <- app(rast(bwkr_055count_jul[,-3:-7]), fun=function(x){ x[x <= 15] <- NA; return(x)})
bwkr_055_jul_binary <- extend(bwkr_055_jul_binary,ext)
bwkr_055_jul_binary <- app(bwkr_055_jul_binary, fun=function(x){ x[x > 15] <- 0.55; return(x)})

bwkr_065_jul_binary <- app(rast(bwkr_065count_jul[,-3:-7]), fun=function(x){ x[x <= 15] <- NA; return(x)})
bwkr_065_jul_binary <- extend(bwkr_065_jul_binary,ext)
bwkr_065_jul_binary <- app(bwkr_065_jul_binary, fun=function(x){ x[x > 15] <- 0.65; return(x)})

bwkr_075_jul_binary <- app(rast(bwkr_075count_jul[,-3:-7]), fun=function(x){ x[x <= 15] <- NA; return(x)})
bwkr_075_jul_binary <- extend(bwkr_075_jul_binary,ext)
bwkr_075_jul_binary <- app(bwkr_075_jul_binary, fun=function(x){ x[x > 15] <- 0.75; return(x)})

bwkr_jul_binary <- cover(bwkr_075_jul_binary,bwkr_065_jul_binary)
bwkr_jul_binary <- cover(bwkr_jul_binary,bwkr_055_jul_binary)
bwkr_jul_binary <- cover(bwkr_jul_binary,bwkr_045_jul_binary)
bwkr_jul_binary <- cover(bwkr_jul_binary,bwkr_035_jul_binary)
bwkr_jul_binary <- raster(bwkr_jul_binary)

plt <- levelplot(bwkr_jul_binary,main = "Persistently High Blue Whale HS Pixels (July)",
                 par.settings = myTheme, at =  c(0, 0.25, 0.44, 0.54, 0.64, 0.74, 1))
plt + latticeExtra::layer(sp.lines(world_outline, col="black", lwd=0.5))

# Visualizing relationship between blwh and krill -------------------------
fpath_bwkr <- "~/Dropbox/blwh_krill/bwkr_monthly"
bwkr_all <- read.csv(paste0(fpath_bwkr,"/bwkr_all.csv"))
bwkr_junjul <- data.frame(bwkr_all) %>%
  filter(month == 6 | month == 7) %>%
  group_by(month)
bwkr_jun <- data.frame(bwkr_all) %>%
  filter(month == 6)
bwkr_jul <- data.frame(bwkr_all) %>%
  filter(month == 7)

# Scatter Plots
library("ggpubr")
ggscatter(bwkr_jun, x = "krill", y = "blwh", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Krill (CPUE)", ylab = "Blue Whale (HS)") +
  ggtitle("Blue Whale vs Krill June")

ggscatter(bwkr_jul, x = "krill", y = "blwh", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Krill (CPUE)", ylab = "Blue Whale (HS)") +
  ggtitle("Blue Whale vs Krill July")

ggscatter(filter(bwkr_jul,blwh >=0.35), x = "krill", y = "blwh", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Krill (CPUE)", ylab = "Blue Whale (HS)") +
  ylim(0, 1) +
  ggtitle("Blue Whale vs Krill (blwh > 0.35) July")

ggscatter(filter(bwkr_jul,blwh >=0.55), x = "krill", y = "blwh", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Krill (CPUE)", ylab = "Blue Whale (HS)") +
  ylim(0, 1) +
  ggtitle("Blue Whale vs Krill (blwh > 0.55) July")

ggscatter(filter(bwkr_jul,blwh >=0.75), x = "krill", y = "blwh", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Krill (CPUE)", ylab = "Blue Whale (HS)") +
  ylim(0, 1) +
  ggtitle("Blue Whale vs Krill (blwh > 0.75) July")

# Plotting Pearson's Coefficient for all 155 months
pcoef <- data.frame(matrix(NA,1,155))
for (i in 1:155) {
  pcoef[[i]] <- cor(bwkr_monthly[[i]]$blwh, bwkr_monthly[[i]]$krill, method = 'pearson')
}
plot(1:155,pcoef)
lines(1:155,pcoef)
