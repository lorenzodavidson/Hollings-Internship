## create ROMs RASTER
#rm(list=ls())
library(raster) # ------------------- > load these directly in parent script (lines 2-3)
library(ncdf4)
library(mondate)

create_ROMS_RASTER=function(nc,dname,month_nc,template=temp){  # fast verion, use this
  nc.data=nc_open(nc)
  #print("grabbing dimensions and variables")
  lat <- ncvar_get(nc.data,'Latitude')
  lon <- ncvar_get(nc.data,'Longitude')
  nrows <- length(lat); ncols <- length(lon)
  tim <- format(as.Date(mondate(anytime(ncvar_get(nc.data,'Time')))), "%Y-%m")
  tmp.array <- ncvar_get(nc.data, dname)
  fillvalue <- ncatt_get(nc.data, dname, "_FillValue")
  #print("creating matrix for month of interest")
  names=unlist(lapply(tim,function(x)gsub(" UTC","",x)))
  index=grep(month_nc,names)
  if(length(index)==0){
    #print("Variable not available for date of interest")
    r2=NA
    return(r2)
  } else {
    
    tmp.array.day=tmp.array[,,index]
    tmp.array.day[tmp.array.day==fillvalue$value]=NA #setting fill value
    
    dat1=list()
    dat1$x=lon
    dat1$y=lat
    dat1$z=t(tmp.array.day)
    
    r <-raster(
      dat1$z,
      xmn=range(dat1$x)[1], xmx=range(dat1$x)[2],
      ymn=range(dat1$y)[1], ymx=range(dat1$y)[2], 
    )
    r2=flip(r,2)
    #???plot(r2)
    #Force to template res
    r3 <- raster::resample(r2, template, method="bilinear")
  }
}


months <- c("03","04","05","06","07")
years <- c(1990:2020)
dates <- vector("character", 155)
for (i in seq_along(years)) {
  for (j in seq_along(months)) {
    x <- i*5 - 5 + j
    if (i == 1 & j == 1) {
      dates[[x]] <- as.character(paste0(years[[i]],"-",months[[j]]))
    } else {
      dates[[x]] <- rbind(as.character(paste0(years[[i]],"-",months[[j]])))
    }
  }
}

krill_monthly_xyz <- list()
krill_monthly_raster <- list()
krill_monthly_xyz_2 <- list()
krill_monthly_raster_2 <- list()

nc = 'TotalKril_CPUE.nc'
dname = 'TotalKrill_CPUE'
temp <- raster(ncol=185, nrow=180, xmn=-133.95, xmx=-115.55, ymn=30.05, ymx=47.95)
template <- blwh_monthly_raster[[1]]
for (t in seq_along(dates)) {
  month_nc <- dates[[t]]
  
  krill_monthly_raster_2[[t]] <- create_ROMS_RASTER(nc,dname,month_nc)
  
  krill_monthly_xyz_2[[t]] <- data.frame(na.omit(rasterToPoints(krill_monthly_raster_2[[t]]))) %>%
    rename("lon" = "x") %>%
    rename("lat" = "y") %>%
    rename("krill" = "layer")
}

#### demo of run
#library(ncdf4)
#nc="D:\\Daily\\bf_daily_roms_gfdl_1980_2100.nc"
#dname="bf"  ## find d name by running nc=nc_open("filename.nc"); print(nc)
#day_nc="1980-02-01"  ## format yyyy-mm-dd  (as.character)
#template= raster("C:/Users/nereo/Dropbox (Personal)/Eco-ROMS/ROMS & Bathym Data/Bathymetry ETOPO1/template.grd")
#windows(20,20)
#raster=create_ROMS_RASTER(nc=nc,dname=dname,day_nc = day_nc,template=template)
#plot(raster)
#raster
