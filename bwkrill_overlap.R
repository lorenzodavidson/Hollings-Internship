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

# Creating bwkr_monthly dataset
bwkr_monthly <- list()
for (t in seq_along(time)) {
  bwkr_monthly[[t]] <- inner_join(blwh_monthly_xyz[[t]],krill_monthly_xyz[[t]]) %>%
    rename("blwh" = "layer")
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

# Save all 6 datasets as csv or raster files
fpath_bwkr <- "~/Dropbox/blwh_sst_monthly/blwh_krill/bwkr_monthly"
fpath_blwh <- "~/Dropbox/blwh_sst_monthly/blwh_krill/blwh_monthly"
fpath_krill <- "~/Dropbox/blwh_sst_monthly/blwh_krill/krill_monthly"
for (i in seq_along(years)) {
  for (j in seq_along(months)) {
    x <- i*5 - 5 + j
    bwkr_name <- paste0("bwkr_", years[[i]], "_", months[[j]],".csv")
    write.csv(bwkr_monthly[[x]],paste0(fpath_bwkr,"/",bwkr_name),row.names = FALSE)
    
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

# Scatter Plots
library("ggpubr")
ggscatter(bwkr_t3, x = "krill", y = "blwh", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Krill (CPUE)", ylab = "Blue Whale (HS)") +
  ggtitle("Blue Whale vs Krill 1990 May")

ggscatter(bwkr_t5, x = "krill", y = "blwh", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Krill (CPUE)", ylab = "Blue Whale (HS)") +
  ggtitle("Blue Whale vs Krill 1990 July")

# Plotting Pearson's Coefficient for all 155 months
pcoef <- data.frame(matrix(NA,1,155))
for (i in 1:155) {
  pcoef[[i]] <- cor(bwkr_monthly[[i]]$blwh, bwkr_monthly[[i]]$krill, method = 'pearson')
}
plot(1:155,pcoef)
lines(1:155,pcoef)

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

# Visualize area overlap --------------------------------------------------

AO_filtered <- data.frame(bwkr_all) %>% # make new dataframe called dailyzoom2015 which is a copy of sdm2015 and do the following below
  mutate(blwh_core = ifelse(blwh >= 0.28,1,0), krill_core = ifelse(krill >= krill_75p_thresh,1,0), Area=1) %>% # create a binary data column for range and area overlap of each species
  filter(blwh_core == 1 & krill_core == 1)

ggplot(AO_filtered, aes(x=month)) + 
  geom_histogram(binwidth=1) # Almost all are in June/July

AO_junejuly <- AO_filtered[,-8:-9] %>%
  filter(month == "6" | month == "7" )
dates_filtered <- unique(AO_junejuly$date)

ext <- extent(-133.95, -115.55, 30.05, 47.95)
AO_junejuly_raster <- list()
for (i in 1:length(dates_filtered)) {
  AO_rasterprep <- filter(AO_junejuly, date == dates_filtered[[i]])
  AO_rasterprep <- rasterFromXYZ(AO_rasterprep[,-3:-7]) 
  AO_junejuly_raster[[i]] <- extend(AO_rasterprep, ext)
}


AO_filtered_t1_ext <- projectRaster(AO_junejuly[[1]], temp, method="bilinear")

AO_filtered_t1_ext <- 
plot(AO_junejuly_raster[[1]])
plot(AO_filtered_t1_ext)





AO_monthly_raster <- list()
for (i in 1:length(dates_filtered)) {
  AO_rasterprep <- filter(AO_filtered, date == dates_filtered[[i]])
  AO_rasterprep <- AO_rasterprep[,-3:-7] 
  AO_monthly_raster[[i]] <- rasterFromXYZ(AO_rasterprep)
}

plot(blwh_monthly_raster[[5]])
plot(AO_monthly_raster[[1]])

for (i in 1:nrow(sst.rasters)) {
  sst.raster <- raster(sst.rasters$file_paths[i])
  blwh.raster <- raster(whale.rasters$file_paths[i])
  png(file = paste0("C:/Users/kaila/Dropbox/Monthly Hindcast Blue Whale Data 1981-2021/combined/", as.character(sst.rasters$year[i]), ".", as.character(sst.rasters$month[i]), ".png"))
  plot(sst.raster, xlab=as.character(sst.rasters$year[i]), main="Sea Surface Temperature and Blue Whale Habitat Suitability > 0.3 (Green)", col=grey(1:50/50))
  plot(blwh.raster, alpha=0.5, add=T, breaks=c(0.3,1), legend=F, col=terrain.colors(3))
  dev.off()
}

r <- raster(xmn=-150, xmx=-120, ymx=60, ymn=30, ncol=36, nrow=18)
values(r) <- 1:ncell(r)

re <- extend(r, e)

# extend with a number of rows and columns (at each side)
re2 <- extend(r, c(2,10))

# Extent object
e <- extent(r)
e
extend(e, 10)
extend(e, 10, -10, 0, 20)
e + 10
e * 2
