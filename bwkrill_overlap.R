library(anytime)


# Creating Yearly Whale+Krill Files ---------------------------------------


file_names <- c()

fpath <- "~/Dropbox/blwh_sst_monthly/blwh/"
dfnames <- list.files(fpath)[stringr::str_detect(dfnames,".grd")]

keywords1 <- c("_03_","_04_","_05_","_06_","_07_")
years <- c(1990:2020)
keywords2 <- as.character(years)
# Filter for years and months
dfnames <- dfnames[Reduce("|", lapply(keywords1, function(x) grepl(x, dfnames)))]
dfnames <- dfnames[Reduce("|", lapply(keywords2, function(x) grepl(x, dfnames)))]

# Create a list for each year containing dataframes for each of the 5 months
for(i in dfnames[stringr::str_detect(dfnames, ".grd")]){
  tt <- raster(paste0(fpath,"/",i))
  
  blwh_monthly[[i]] <- tt
}

cl
blwh_t1 <- data.frame(rasterToPoints(blwh_monthly[[1]]))
krill_t1test <- data.frame(rasterToPoints(krill_monthly[[1]]))

blwh_t1test <- inner_join(blwh_t1,krill_t1test, by = c("x", y))

krill_wtf <- krill_t1test %>%
  filter(x == -127.05 | y == 47.95)

blwh_t1[[1,1]]

# Opening whale and krill data into year files --------------------------------------------
fpath <- "~/Dropbox/blwh_sst_monthly/blwh/" 
dfnames <- list.files(fpath)
blwh_monthly<- list()
for(i in dfnames[stringr::str_detect(dfnames, ".grd")]){
  tt <- raster(paste0(fpath,"/",i))
  blwh_monthly[[i]] <- tt
}

totalkrill <- nc_open('TotalKril_CPUE.nc')
lon <- ncvar_get(totalkrill, "Longitude")
lat <- ncvar_get(totalkrill, "Latitude")
time <- anytime(ncvar_get(totalkrill, "Time"))
TotalKrill_CPUE <- ncvar_get(totalkrill, "TotalKrill_CPUE")

krill_t1 <- TotalKrill_CPUE[,,1]
krill_t1_transpose <- t(krill_t1)
krill_t1_rasterprep <- matrix(NA,33300,3)


krill_monthly_long<- list()
krill_monthly_raster<- list()
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
  krill_monthly_long[[t]]<-krill_rasterprep
  #krill_raster <- rasterFromXYZ(krill_rasterprep,crs=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  #krill_monthly_raster[[t]] <- krill_raster
}

writeRaster(krill_monthly, "~/Desktop/Hollingsinternship/Code/Hollings-Internship/krill_monthly.nc", format="CDF")

write.csv(krill_monthly, "~/Desktop/Hollingsinternship/Code/Hollings-Internship/krill_t1.csv", row.names = FALSE)


krill_t1_rasterprep1 <- c(lon[[1]],lat[[1]],krill_t1[[1,1]])
new_row <- c(lon[[1]],lat[[2]],krill_t1[[1,2]])
krill_t1_rasterprep1 <- rbind(krill_t1_rasterprep1, new_row)
krill_t1_raster <- rasterFromXYZ(krill_t1_rasterprep1,crs=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

plot(krill_t1_raster)

krill_monthly<- list()
for (i in seq_along(time)) {
  x1 <- TotalKrill_CPUE[,,i]
  x2 <- t(x1)
  krill_monthly[[i]] <- raster(x2,
                               xmn=range(lon[1]), xmx=range(lon[2])
                               ymn=range(lat[1]), ymx=range(lat[2]),
                               crs=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
                               )
}

dat1=list()
dat1$x=lon
dat1$y=lat
dat1$z=t(tmp.array.day)

r <-raster(
  dat1$z,
  xmn=range(dat1$x)[1], xmx=range(dat1$x)[2],
  ymn=range(dat1$y)[1], ymx=range(dat1$y)[2], 
  crs=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
)


krill_t1 <- TotalKrill_CPUE[,,1]
krill_t1_transpose <- t(krill_t1)
plot(krill_monthly[[2]])

# Flipping Krill, Lat and Lon Arrays for Visualization
krill_flipped <- array(dim = c(185,180,155))
for (t in 1:155) {
  for (i in 1:185) {
    for (j in 1:180) {
      krill_flipped[[i,181-j,t]] <- krill_monthly[[i,j,t]]
    }
  }
}

lon_flipped <- array(NA,185)
for (i in 1:185) {
  lon_flipped[[186-i]] <- -lon[[i]]
}

lat_flipped <- array(NA,180)
for (j in 1:180) {
  lat_flipped[[181-j]] <- lat[[j]]
}

# Rasterizing Krill Data
krill_t1 <- krill_monthly[,,1]
krill_t1_flipped <- krill_flipped[,,1]

rotate <- function(x) t(apply(x, 2, rev))

krill_t1_transpose <- t(krill_t1)

plot(raster(krill_t1))
plot(raster(krill_t1_flipped))
plot(raster(krill_t1_rotated))
plot(raster(krill_t1_transpose))



krill_t1_raster <- data_frame()
raster()
for (i in seq_along(lon)) {
  for (j in seq_along(lat)) {
    
  }
}



for (i in seq_along(t)) {
  
}




plot(krill_flipped[,,6])

image(lon_flipped,lat_flipped,krill_flipped[,,6])






lon_cal <- matrix(NA,185,180)
lat_cal <- matrix(NA,185,180)
for (i in 1:185) {
  for (j in 1:180) {
    lon_cal[[i,j]] <- lon[[i]]
    lat_cal[[i,j]] <- lat[[j]]
  }
}

write.csv(krill_t1, "~/Desktop/Hollingsinternship/Code/Hollings-Internship/krill_t1.csv", row.names = FALSE)
write.csv(lon_cal, "~/Desktop/Hollingsinternship/Code/Hollings-Internship/lon_cal.csv", row.names = FALSE)
write.csv(lat_cal, "~/Desktop/Hollingsinternship/Code/Hollings-Internship/lat_cal.csv", row.names = FALSE)


himage(lon,lat,krill_t1)

plot(krill_t1)

data(wrld_simpl)
plot(wrld_simpl, xlim=c(-120,-80), ylim=c(0,60), axes=TRUE, col="light yellow")
points(non_anomalous, pch=16, col=rgb(0, 0, 1, alpha=0.5), cex=0.3)
points(anomalous, pch=16, col=rgb(1, 0, 0, alpha=0.5), cex=0.3)