library(anytime)

# Opening whale and krill data --------------------------------------------
fpath <- "~/Dropbox/blwh_sst_monthly/blwh/" 
dfnames <- list.files(fpath)
dflist <- list()
for(i in dfnames[stringr::str_detect(dfnames, ".grd")]){
  tt <- raster(paste0(fpath,"/",i))
  dflist[[i]] <- tt
}

totalkrill <- nc_open('TotalKril_CPUE.nc')

lon <- ncvar_get(totalkrill, "Longitude")
lat <- ncvar_get(totalkrill, "Latitude", verbose = F)
t <- anytime(ncvar_get(totalkrill, "Time"))
totalkrillCPUE <- ncvar_get(totalkrill, "TotalKrill_CPUE")

krill_flipped <- array(dim = c(185,180,155))
for (t in 1:155) {
    for (i in 1:185) {
      for (j in 1:180) {
        krill_flipped[[i,181-j,t]] <- totalkrillCPUE[[i,j,t]]
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