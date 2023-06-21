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

# Create vector with all dates
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

# All values for 155 months
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

# Freq plots for 1990 -----------------------------------------------------
bwkr_t1 <- bwkr_monthly[[1]]
bwkr_t2 <- bwkr_monthly[[2]]
bwkr_t3 <- bwkr_monthly[[3]]
bwkr_t4 <- bwkr_monthly[[4]]
bwkr_t5 <- bwkr_monthly[[5]]

ggplot(bwkr_t1, aes(blwh)) +
  geom_freqpoly(binwidth = 0.01) +
  xlim(0, 1) +
  ylim(0, 250) +
  ggtitle("1990 March")

ggplot(bwkr_t2, aes(blwh)) +
  geom_freqpoly(binwidth = 0.01) +
  xlim(0, 1) +
  ylim(0, 250) +
  ggtitle("1990 April")

ggplot(bwkr_t3, aes(blwh)) +
  geom_freqpoly(binwidth = 0.01) +
  xlim(0, 1) +
  ylim(0, 250) +
  ggtitle("1990 May")

ggplot(bwkr_t4, aes(blwh)) +
  geom_freqpoly(binwidth = 0.01) +
  xlim(0, 1) +
  ylim(0, 250) +
  ggtitle("1990 June")

ggplot(bwkr_t5, aes(blwh)) +
  geom_freqpoly(binwidth = 0.01) +
  xlim(0, 1) +
  ylim(0, 250) +
  ggtitle("1990 July")

ggplot(bwkr_t1, aes(krill)) +
  geom_freqpoly(binwidth = 0.05) +
  xlim(0, 10.5) +
  ylim(0, 100) +
  ggtitle("Krill CPUE 1990 March")

ggplot(bwkr_t2, aes(krill)) +
  geom_freqpoly(binwidth = 0.05) +
  xlim(0, 10.5) +
  ylim(0, 100) +
  ggtitle("Krill CPUE 1990 April")

ggplot(bwkr_t3, aes(krill)) +
  geom_freqpoly(binwidth = 0.05) +
  xlim(0, 10.5) +
  ylim(0, 100) +
  ggtitle("Krill CPUE 1990 May")

ggplot(bwkr_t4, aes(krill)) +
  geom_freqpoly(binwidth = 0.05) +
  xlim(0, 10.5) +
  ylim(0, 100) +
  ggtitle("Krill CPUE 1990 June")

ggplot(bwkr_t5, aes(krill)) +
  geom_freqpoly(binwidth = 0.05) +
  xlim(0, 10.5) +
  ylim(0, 100) +
  ggtitle("Krill CPUE 1990 July")

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
# Investigating threshold values ------------------------------------------

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
AO_bymonth <- group_by(AO_all,month)
RO_all <- filter(bwkr_metrics,Metric_Name=="RO")
RO_bymonth <- group_by(RO_all,month)
Schoener_all <- filter(bwkr_metrics,Metric_Name=="Schoener")
Schoener_bymonth <- group_by(Schoener_all,month)
Bhatty_all <- filter(bwkr_metrics,Metric_Name=="Bhatty")
Bhatty_bymonth <- group_by(Bhatty_all,month)

# AO Plots
ggplot(AO_all,aes(x=date,y=Overlap_Metric,group = 1)) +
  geom_point() + 
  geom_line() +
  ggtitle("AO 1990-2020") +
  ylab("Area Overlap")
ggplot(AO_bymonth,aes(x=month,y=Overlap_Metric)) +
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
ggplot(RO_bymonth,aes(x=month,y=Overlap_Metric)) +
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
ggplot(Schoener_bymonth,aes(x=month,y=Overlap_Metric)) +
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
ggplot(Bhatty_bymonth,aes(x=month,y=Overlap_Metric)) +
  geom_boxplot() +
  ggtitle("Bhattacharyya's Coefficient grouped by month") +
  xlab("Month (March-July)") + 
  ylab("Bhattacharyya's Coefficient")
