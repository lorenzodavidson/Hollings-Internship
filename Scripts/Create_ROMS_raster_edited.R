## create ROMs RASTER
#rm(list=ls())
library(raster) # ------------------- > load these directly in parent script (lines 2-3)
library(ncdf4)
library(mondate)

create_ROMS_RASTER=function(nc,dname,month_nc,template){
  nc.data=nc_open(nc)
  lat <- ncvar_get(nc.data,'Latitude')
  lon <- ncvar_get(nc.data,'Longitude')
  nrows <- length(lat); ncols <- length(lon)
  tim <- format(as.Date(mondate(anytime(ncvar_get(nc.data,'Time')))), "%Y-%m")
  tmp.array <- ncvar_get(nc.data, dname)
  fillvalue <- ncatt_get(nc.data, dname, "_FillValue")
  names=unlist(lapply(tim,function(x)gsub(" UTC","",x)))
  index=grep(month_nc,names)
  
  # Create rasters
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
  r <- raster::resample(r, template, method="bilinear")
}

# Example for krill data
nc = "~/Dropbox/krill_raw/TotalKril_CPUE.nc"
dname = 'TotalKrill_CPUE'
temp <- raster(ncol=185, nrow=180, xmn=-134, xmx=-115.5, ymn=30, ymx=48)

krill_monthly_raster <- list()
for (t in seq_along(dates)) {
  month_nc <- dates[[t]]
  
  krill_monthly_raster[[t]] <- create_ROMS_RASTER(nc,dname,month_nc,temp)
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
