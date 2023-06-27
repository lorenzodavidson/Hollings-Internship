## create ROMs RASTER
#rm(list=ls())
library(raster) # ------------------- > load these directly in parent script (lines 2-3)
library(ncdf4)

create_ROMS_RASTER=function(nc,dname,day_nc,template=template){  # fast verion, use this
  nc.data=nc_open(nc)
  print("grabbing dimensions and variables")
  lat <- ncvar_get(nc.data,'lat')
  lon <- ncvar_get(nc.data,'lon')
  nrows <- length(lat); ncols <- length(lon)
  yr <- ncvar_get(nc.data,'year'); mth <- ncvar_get(nc.data,'month'); day <- ncvar_get(nc.data,'day')
  tim <- as.POSIXct(paste(yr,mth,day,sep='-'),tz='UTC')
  tmp.array <- ncvar_get(nc.data, dname)
  fillvalue <- ncatt_get(nc.data, dname, "_FillValue")
  print("creating matrix for day of interest")
  names=unlist(lapply(tim,function(x)gsub(" UTC","",x)))
  index=grep(day_nc,names)
  if(length(index)==0){
    print("Variable not available for date of interest")
    r3=NA
    return(r3)
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
      crs=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
    )
    r2=flip(r,2)
    #???plot(r2)
    #Force to template res
    r3 <- raster::resample(r2, template, method="bilinear")  
    
    return(r3)
  }
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
