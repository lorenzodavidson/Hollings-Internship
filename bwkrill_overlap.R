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

# Area Overlap vs threshold value
krill_threshold_t3 <- unname(quantile(bwkr_monthly[[3]]$krill, probs = seq(.01, 1, by = .01)))
krill_threshold_t5 <- unname(quantile(bwkr_monthly[[5]]$krill, probs = seq(.01, 1, by = .01)))
blwh_threshold <- c(0.15,0.25,0.35,0.45)

AO_bwkr_t3 <- data.frame(matrix(NA,length(krill_threshold_t3),length(blwh_threshold)+1)) %>%
  rename("krill_t" = "X1") %>%
  rename("blwh_t15" = "X2") %>%
  rename("blwh_t25" = "X3") %>%
  rename("blwh_t35" = "X4") %>%
  rename("blwh_t45" = "X5")
for (j in 1:4) {
  for (i in seq_along(krill_threshold_t3)) {
    AO_bwkr_t3[[i,1]] <- krill_threshold_t3[[i]]
    AO_bwkr_t3[[i,j+1]] <- bwkr_monthly[[3]] %>%
      mutate(blwh_core = ifelse(blwh >= blwh_threshold[[j]],1,0), krill_core = ifelse(krill >= krill_threshold_t3[[i]],1,0), Area=1) %>%
      summarise(AO=area_overlapfn(prey=krill_core,pred=blwh_core,area=Area)) %>%
      unlist()
  }
}

AO_bwkr_t5 <- data.frame(matrix(NA,length(krill_threshold_t5),length(blwh_threshold)+1)) %>%
  rename("krill_t" = "X1") %>%
  rename("blwh_t15" = "X2") %>%
  rename("blwh_t25" = "X3") %>%
  rename("blwh_t35" = "X4") %>%
  rename("blwh_t45" = "X5")
for (j in 1:4) {
  for (i in seq_along(krill_threshold_t5)) {
    AO_bwkr_t5[[i,1]] <- krill_threshold_t5[[i]]
    AO_bwkr_t5[[i,j+1]] <- bwkr_monthly[[5]] %>%
      mutate(blwh_core = ifelse(blwh >= blwh_threshold[[j]],1,0), krill_core = ifelse(krill >= krill_threshold_t5[[i]],1,0), Area=1) %>%
      summarise(AO=area_overlapfn(prey=krill_core,pred=blwh_core,area=Area)) %>%
      unlist()
  }
}

plot(rep(AO_bwkr_t3$krill_t[75],5),AO_bwkr_t3[75,],main="May 1990 AO vs Krill Threshold",
     xlab="Krill Threshold", ylab="Area Overlap", ylim=c(0,1),xlim=c(2,11))
lines(AO_bwkr_t3$krill_t,AO_bwkr_t3$blwh_t15,col="red")
lines(AO_bwkr_t3$krill_t,AO_bwkr_t3$blwh_t25,col="orange")
lines(AO_bwkr_t3$krill_t,AO_bwkr_t3$blwh_t35,col="green")
lines(AO_bwkr_t3$krill_t,AO_bwkr_t3$blwh_t45,col="blue")
legend(7.5, 0.95, legend=c("75th Percentile", "bw_tresh: 0.15", "bw_tresh: 0.25","bw_tresh: 0.35","bw_tresh: 0.45"),
       col=c("black","red","orange","green","blue"),lty=1:2, cex=0.8)

plot(rep(AO_bwkr_t5$krill_t[75],5),AO_bwkr_t5[75,],main="July 1990 AO vs Krill Threshold",
     xlab="Krill Threshold", ylab="Area Overlap",ylim=c(0,1),xlim=c(2,11))
lines(AO_bwkr_t5$krill_t,AO_bwkr_t5$blwh_t15,col="red")
lines(AO_bwkr_t5$krill_t,AO_bwkr_t5$blwh_t25,col="orange")
lines(AO_bwkr_t5$krill_t,AO_bwkr_t5$blwh_t35,col="green")
lines(AO_bwkr_t5$krill_t,AO_bwkr_t5$blwh_t45,col="blue")
legend(7.5, 0.95, legend=c("75th Percentile","bw_tresh: 0.15", "bw_tresh: 0.25","bw_tresh: 0.35","bw_tresh: 0.45"),
       col=c("black","red","orange","green","blue"),lty=1:2, cex=0.8)

# All 4 metrics

Overlap_bwkr_t5 <- data.frame(matrix(NA,length(krill_threshold_t5),17))
colnames(Overlap_bwkr_t5) <- c("krill_t", "AO_bw15","AO_bw25","AO_bw35","AO_bw45","RO_bw15","RO_bw25","RO_bw35","RO_bw45",
                    "Schoener_bw15","Schoener_bw25","Schoener_bw35","Schoener_bw45",
                    "Bhatty_bw15","Bhatty_bw25","Bhatty_bw35","Bhatty_bw45")
for (j in 1:4) {
  for (i in seq_along(krill_threshold_t5)) {
    Overlap_bwkr_t5[[i,1]] <- krill_threshold_t5[[i]]
    AO_bwkr_t5[[i,j+1]] <- bwkr_monthly[[5]] %>%
      mutate(blwh_core = ifelse(blwh >= blwh_threshold[[j]],1,0), krill_core = ifelse(krill >= krill_threshold_t5[[i]],1,0), Area=1) %>%
      summarise(AO=area_overlapfn(prey=krill_core,pred=blwh_core,area=Area),
                RO=range_overlapfn(prey=krill_core,pred=blwh_core,area=Area),
                Schoener=schoeners_overlapfn(prey=krill,pred=blwh),
                Bhatty=bhatta_coeffn(prey=krill,pred=blwh)) %>%
      unlist()
  }

# Applying overlap metrics to actual data ---------------------------------




percentile75 <- data.frame(matrix(NA,155,2))
for (i in 1:155) {
  percentile75[[i,1]] <- quantile(bwkr_monthly[[i]]$blwh, c(0.75))
  percentile75[[i,2]] <- quantile(bwkr_monthly[[i]]$krill, c(0.75))             
}
percentile75_avg <- data.frame(matrix(NA,5,2))
for (i in 1:31) {
  percentile75_avg[[i,1]] <- mean(percentile75[[]])
  percentile75_avg[[i,2]] <- quantile(bwkr_monthly[[i]]$krill, c(0.75))             
}
lines(1:155,percentile75$X2)




#####

overlap_bwkr_t1 <- bwkr_t1 %>%
  mutate(blwh_core = ifelse(blwh >= 0.28,1,0), krill_core = ifelse(krill >= 0.5,1,0), Area=1) %>%
  # filter(Lat>=34 & Lat < 35,Lon <= -119 & Lon > -120) %>%
  summarise(AO=area_overlapfn(prey=krill_core,pred=blwh_core,area=Area),
            RO=range_overlapfn(prey=krill_core,pred=blwh_core,area=Area),
            Schoener=schoeners_overlapfn(prey=krill,pred=blwh),
            Bhatty=bhatta_coeffn(prey=krill,pred=blwh)) %>%
  #pivot_longer(cols=c(AO,RO,Schoener,Bhatty),names_to="Metric_Name",values_to="Overlap_Metric")


plot(blwh_t2_cut)
plot(krill_monthly_raster[[1]])

daily2015zoom <- sdm2015 %>% # make new dataframe called dailyzoom2015 which is a copy of sdm2015 and do the following below
  mutate(Humpback_Core = ifelse(Humpback_HS >= 0.28,1,0),Anchovy_Core = ifelse(Anchovy_HS >= 0.5,1,0),Area=1) %>% # create a binary data column for range and area overlap of each species
  filter(Lat>=34 & Lat < 35,Lon <= -119 & Lon > -120) %>% # retain values only within this lat-lon range
  group_by(Date) %>% # organize by date
  summarise(AO=area_overlapfn(prey=Anchovy_Core,pred=Humpback_Core,area=Area),RO=range_overlapfn(prey=Anchovy_Core,pred=Humpback_Core,area=Area), # apply the overlap metrics and output as new dataframe
            Schoener=schoeners_overlapfn(prey=Anchovy_HS,pred=Humpback_HS), Bhatty=bhatta_coeffn(prey=Anchovy_HS,pred=Humpback_HS)) %>% 
  pivot_longer(cols=c(AO,RO,Schoener,Bhatty),names_to="Metric_Name",values_to="Overlap_Metric") # change the dataframe into long format for easier plotting



blwh_t1_cut <- rasterFromXYZ(select(bwkr_t1,x,y,blwh),crs=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
krill_t1_cut <- rasterFromXYZ(select(bwkr_t1,x,y,krill),crs=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
blwh_t2_cut <- rasterFromXYZ(select(bwkr_t2,x,y,blwh),crs=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
krill_t2_cut <- rasterFromXYZ(select(bwkr_t2,x,y,krill),crs=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

