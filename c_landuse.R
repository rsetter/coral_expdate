library(raster)
library(ncdf4)
library(doParallel)
library(velox)
library(sf)
library(foreach)
library(rgis)
library(dplyr)
library(tibble)


#empirical
#use HYDE 3.2
#metadata: https://dataportaal.pbl.nl/downloads/HYDE/HYDE3.2/readme_release_HYDE3.2.1.txt


#2005 data

setwd("land_hurtt/empirical")

#import needed files
areahyde <- read.asciigrid("garea_cr.asc") #total gridcell area in km2
  areah <- raster(areahyde)
  #cropland
  crop <- read.asciigrid("cropland2005AD.asc") #total cropland area, in km2 per grid cell
  crops <- raster(crop)
  c <- crops / areah #calculate the percentage of each cell that is cropland
  #pasture
  pasture <- read.asciigrid("pasture2005AD.asc") #total pasture area, in km2 per grid cell
  pastures <- raster(pasture)
  p <- pastures / areah #calculate the percentage of each cell that is pasture
  #urban
  urban <- read.asciigrid("uopp_2005AD.asc") #total built-up area, such as towns, cities, etc, in km2 per grid cell
  urbans <- raster(urban)
  u <- urbans / areah #calculate the percentage of each cell that is urban

#calculate total unsuitable land
tu <- overlay(c,p,u,fun=function(c,p,u,na.rm=T){return(c+p+u)},filename="emp_unsuitland.tif")
  
####
#extract reef suit vals
#find maximum unsuitable land % within 50km radius of reef site - save max % unsuitability as raster
#classify as suitable or not - save as .csv
###

polyfiles <- list.files("reef_polygon",pattern="\\.shp$",full.names=T)
nfiles <- length(polyfiles)
reefs <- st_read("currentcoral/reefs.shp")
IDmatch <- read.csv("land_hurtt/polygonIDmatch.csv",header=T)


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)
  land <- tu
  df.extract1 <- foreach::foreach(i=1:5, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract2 <- foreach::foreach(i=6:10, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract3 <- foreach::foreach(i=11:15, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract4 <- foreach::foreach(i=16:20, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract5 <- foreach::foreach(i=20:25, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract6 <- foreach::foreach(i=26:30, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  dff <- rbind(df.extract1,df.extract2)
  dff1 <- rbind(dff,df.extract3)
  dff2 <- rbind(dff1,df.extract4)
  dff3 <- rbind(dff2,df.extract5)
  dff4 <- rbind(dff3,df.extract6)
  dfext <- as.data.frame(dff4)
  dfext <- subset(dfext,select=-c(grid_code, BUFF_DIST, ORIG_FID,Shape_Leng,Shape_Area,geometry))
  names(dfext) <- c("pointid","zone","maxunsuit")
  result <- transform(dfext,pointid=as.numeric(pointid),maxunsuit=as.numeric(maxunsuit))
  result$suit <- with(result,ifelse(maxunsuit <0.5, 1, 0)) #if the percentage of unsuitable land is less than 0.5, it is suitable (1)
  result$suit <- with(result, ifelse(is.na(maxunsuit),1,suit)) #if land is NA, then land is too far away to return value so is suitable (1)
  b<-merge(IDmatch,result,by="pointid",all.x=T,all.y=F)
  c <- b[!duplicated(b$ID),]
  d <- c[base::order(c$ID),]
  write.csv(d,"emp_2005_suit.csv", row.names=FALSE)
  gc()
  
stopCluster(cl)
finish <- Sys.time()
finish-start

#get count of suitability
  b <- d[!(d$iscoral==0),]
  count(b, vars=suit)






#LUH model

#metadata: https://daac.ornl.gov/VEGETATION/guides/Land_Use_Harmonization_V1.html
#
#historic = LUHa.t1u2.v1
#rcp26 = IMAGE
#rcp45 = GCAM
#rcp85 = MESSAGE
#proportion landcover of urban+cropland+pasture is unsuitable (primary and secondary and pasture land are suitable)
  #gcrop, gurbn, gpast



#historic

setwd("land_hurtt/historic")

#input unsuitable-land data. keep just 1851-2005 (this dataset is 1700-2005)
c <- brick("land_hurtt/LAND_USE_HARMONIZATION_V1_1248/data/LUHa_u2t1.v1_gcrop.nc4")
c <- dropLayer(c,c(0:151))
p <- brick("land_hurtt/LAND_USE_HARMONIZATION_V1_1248/data/LUHa_u2t1.v1_gpast.nc4")
p <- dropLayer(c,c(0:151))
u <- brick("land_hurtt/LAND_USE_HARMONIZATION_V1_1248/data/LUHa_u2t1.v1_gurbn.nc4")
u <- dropLayer(u,c(0:151))

#calculate total unsuitable land
tu <- overlay(c,p,u,fun=function(c,p,u,na.rm=T){return(c+p+u)},filename="hist_unsuitland.tif")

#calculate 10 year averages
for (i in 1:16) {
  index <- c(1:nlayers(tu))
  ten <- c((i*10-9):(i*10))
  window <- intersect(index,ten)
  win <- dropLayer(tu,c(0:ten[1]-1,ten[length(ten)]+1:nlayers(tu)))
  new10 <- calc(win,fun=mean,na.rm=T)
  tu10 <- stack(tu10,new10)
}
writeRaster(tu10,filename="hist_unsuitland10.tif")



####
#extract reef suit vals
#find maximum unsuitable land % within 50km radius of reef site - save max % unsuitability as raster
#classify as suitable or not - save as .csv
###

polyfiles <- list.files("reef_polygon",pattern="\\.shp$",full.names=T)
nfiles <- length(polyfiles)
reefs <- st_read("currentcoral/reefs.shp")
IDmatch <- read.csv("land_hurtt/polygonIDmatch.csv",header=T)
histsuitcsvout <- paste("hist_",1850+(1:nlayers(tu10)-1)*10+5,"_suit.csv",sep="")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)
for (i in 1:nlayers(tu10)){
  land <- tu10[[i]]
  df.extract1 <- foreach::foreach(i=1:5, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract2 <- foreach::foreach(i=6:10, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract3 <- foreach::foreach(i=11:15, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract4 <- foreach::foreach(i=16:20, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract5 <- foreach::foreach(i=20:25, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract6 <- foreach::foreach(i=26:30, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  dff <- rbind(df.extract1,df.extract2)
  dff1 <- rbind(dff,df.extract3)
  dff2 <- rbind(dff1,df.extract4)
  dff3 <- rbind(dff2,df.extract5)
  dff4 <- rbind(dff3,df.extract6)
  dfext <- as.data.frame(dff4)
  dfext <- subset(dfext,select=-c(grid_code, BUFF_DIST, ORIG_FID,Shape_Leng,Shape_Area,geometry))
  names(dfext) <- c("pointid","zone","maxunsuit")
  result <- transform(dfext,pointid=as.numeric(pointid),maxunsuit=as.numeric(maxunsuit))
  result$suit <- with(result,ifelse(maxunsuit <0.5, 1, 0)) #if the percentage of unsuitable land is less than 0.5, it is suitable (1)
  result$suit <- with(result, ifelse(is.na(maxunsuit),1,suit)) #if land is NA, then land is too far away to return value so is suitable (1)
  landpt <- inner_join(reef,result,by="pointid")
  landdf <- as.data.frame(landpt)
  write.csv(landdf,histsuitcsvout[i], row.names=FALSE)
  gc()
}
stopCluster(cl)
finish <- Sys.time()
finish-start

#get count of suitability
names(IDmatch) <- c("ObjectID","iscoral","ID","pointid")
csvfiles <- list.files("land_hurtt/historic/with_pasture",pattern="_suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
histsuitIDcsvout <- paste("hist_",1845+(1:ncsvfiles)*10,"_suit.csv",sep="")
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  names(a) <- c("pointid","grid_code","zone","maxunsuit","suit","geometryx","geometryy")
  b<-merge(IDmatch,a,by="pointid",all.x=T,all.y=F)
  c <- b[!duplicated(b$ID),]
  d <- c[base::order(c$ID),]
  e <- d[!(d$iscoral==0),]
  f <- count(e, vars=suit)
  y <- add_column(f,year=paste0(i*10+1845),.before=T)
  count <- rbind(y,count)
  write.csv(d,histsuitIDcsvout[i], row.names=FALSE)
  
} 







#rcp26 
# = IMAGE

setwd("land_hurtt/rcp26")

#input unsuitable-land data
c <- brick("land_hurtt/LAND_USE_HARMONIZATION_V1_1248/data/LUHa_u2t1.v1_image.v1.1_gcrop.nc4")
p <- brick("land_hurtt/LAND_USE_HARMONIZATION_V1_1248/data/LUHa_u2t1.v1_image.v1.1_gpast.nc4")
u <- brick("land_hurtt/LAND_USE_HARMONIZATION_V1_1248/data/LUHa_u2t1.v1_image.v1.1_gurbn.nc4")

#calculate total unsuitable land
tu <- overlay(c,p,u,fun=function(c,p,u,na.rm=T){return(c+p+u)},filename="rcp26_unsuitland.tif")

#calculate 10 year averages
tu10<-raster()
for (i in 1:20) {
  index <- c(1:nlayers(tu))
  five <- c(((i*5-4)-2):((i*5-4)+2))
  window <- intersect(index,five)
  win <- dropLayer(tu,c(0:five[1]-1,five[length(five)]+1:nlayers(tu)))
  new10 <- calc(win,fun=mean,na.rm=T)
  tu10 <- stack(tu10,new10)
}
writeRaster(tu10,filename="rcp26_unsuitland10.tif")



####
#extract reef suit vals
#find maximum unsuitable land % within 50km radius of reef site - save max % unsuitability as raster
#classify as suitable or not - save as .csv
###

polyfiles <- list.files("reef_polygon",pattern="\\.shp$",full.names=T)
nfiles <- length(polyfiles)
rcp26suitcsvout <- paste("rcp26_",2000+(1:nlayers(tu10))*5,"_suit.csv",sep="")
tu10 <- brick("rcp26_unsuitland10.tif")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)
for (i in 1:nlayers(tu10)){
  land <- tu10[[i]]
  df.extract1 <- foreach::foreach(i=1:5, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract2 <- foreach::foreach(i=6:10, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract3 <- foreach::foreach(i=11:15, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract4 <- foreach::foreach(i=16:20, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract5 <- foreach::foreach(i=20:25, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract6 <- foreach::foreach(i=26:30, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  dff <- rbind(df.extract1,df.extract2)
  dff1 <- rbind(dff,df.extract3)
  dff2 <- rbind(dff1,df.extract4)
  dff3 <- rbind(dff2,df.extract5)
  dff4 <- rbind(dff3,df.extract6)
  dfext <- as.data.frame(dff4)
  dfext <- subset(dfext,select=-c(grid_code, BUFF_DIST, ORIG_FID,Shape_Leng,Shape_Area,geometry))
  names(dfext) <- c("pointid","zone","maxunsuit")
  result <- transform(dfext,pointid=as.numeric(pointid),maxunsuit=as.numeric(maxunsuit))
  result$suit <- with(result,ifelse(maxunsuit <0.5, 1, 0)) #if the percentage of unsuitable land is less than 0.5, it is suitable (1)
  result$suit <- with(result, ifelse(is.na(maxunsuit),1,suit)) #if land is NA, then land is too far away to return value so is suitable (1)
  b<-merge(IDmatch,result,by="pointid",all.x=T,all.y=F)
  c <- b[!duplicated(b$ID),]
  d <- c[base::order(c$ID),]
  write.csv(d,rcp26suitcsvout[i], row.names=FALSE)
  gc()
}
stopCluster(cl)
finish <- Sys.time()
finish-start

#get count of suitability
csvfiles <- list.files("land_hurtt/rcp26",pattern="_suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  b <- a[!(a$iscoral==0),]
  c <- count(b, vars=suit)
  y <- add_column(c,year=paste0(i*5+2000),.before=T)
  count <- rbind(y,count)
} 








#rcp45 
# = GCAM/minicam

setwd("land_hurtt/rcp45")

#input unsuitable-land data. keep just 1851-2005 (this dataset is 1700-2005)
c <- brick("land_hurtt/LAND_USE_HARMONIZATION_V1_1248/data/LUHa_u2t1.v1_minicam.v1_gcrop.nc4")
p <- brick("land_hurtt/LAND_USE_HARMONIZATION_V1_1248/data/LUHa_u2t1.v1_minicam.v1_gpast.nc4")
u <- brick("land_hurtt/LAND_USE_HARMONIZATION_V1_1248/data/LUHa_u2t1.v1_minicam.v1_gurbn.nc4")

#calculate total unsuitable land
tu <- overlay(c,p,u,fun=function(c,p,u,na.rm=T){return(c+p+u)},filename="rcp45_unsuitland.tif")

#calculate 10 year averages
tu10<-raster()
for (i in 1:20) {
  index <- c(1:nlayers(tu))
  five <- c(((i*5-4)-2):((i*5-4)+2))
  window <- intersect(index,five)
  win <- dropLayer(tu,c(0:five[1]-1,five[length(five)]+1:nlayers(tu)))
  new10 <- calc(win,fun=mean,na.rm=T)
  tu10 <- stack(tu10,new10)
}
writeRaster(tu10,filename="rcp45_unsuitland10.tif")



####
#extract reef suit vals
#find maximum unsuitable land % within 50km radius of reef site - save max % unsuitability as raster
#classify as suitable or not - save as .csv
###

polyfiles <- list.files("reef_polygon",pattern="\\.shp$",full.names=T)
nfiles <- length(polyfiles)
rcp45suitcsvout <- paste("rcp45_",2000+(1:nlayers(tu10))*5,"_suit.csv",sep="")
tu10 <- brick("rcp45_unsuitland10.tif")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)
for (i in 1:nlayers(tu10)){
  land <- tu10[[i]]
  df.extract1 <- foreach::foreach(i=1:5, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract2 <- foreach::foreach(i=6:10, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract3 <- foreach::foreach(i=11:15, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract4 <- foreach::foreach(i=16:20, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract5 <- foreach::foreach(i=20:25, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract6 <- foreach::foreach(i=26:30, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  dff <- rbind(df.extract1,df.extract2)
  dff1 <- rbind(dff,df.extract3)
  dff2 <- rbind(dff1,df.extract4)
  dff3 <- rbind(dff2,df.extract5)
  dff4 <- rbind(dff3,df.extract6)
  dfext <- as.data.frame(dff4)
  dfext <- subset(dfext,select=-c(grid_code, BUFF_DIST, ORIG_FID,Shape_Leng,Shape_Area,geometry))
  names(dfext) <- c("pointid","zone","maxunsuit")
  result <- transform(dfext,pointid=as.numeric(pointid),maxunsuit=as.numeric(maxunsuit))
  result$suit <- with(result,ifelse(maxunsuit <0.5, 1, 0)) #if the percentage of unsuitable land is less than 0.5, it is suitable (1)
  result$suit <- with(result, ifelse(is.na(maxunsuit),1,suit)) #if land is NA, then land is too far away to return value so is suitable (1)
  b<-merge(IDmatch,result,by="pointid",all.x=T,all.y=F)
  c <- b[!duplicated(b$ID),]
  d <- c[base::order(c$ID),]
  write.csv(d,rcp45suitcsvout[i], row.names=FALSE)
  gc()
}
stopCluster(cl)
finish <- Sys.time()
finish-start

#get count of suitability
csvfiles <- list.files("land_hurtt/rcp45",pattern="_suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  b <- a[!(a$iscoral==0),]
  c <- count(b, vars=suit)
  y <- add_column(c,year=paste0(i*5+2000),.before=T)
  count <- rbind(y,count)
} 






#rcp85 
# = MESSAGE

setwd("land_hurtt/rcp85")

#input unsuitable-land data. keep just 1851-2005 (this dataset is 1700-2005)
c <- brick("land_hurtt/LAND_USE_HARMONIZATION_V1_1248/data/LUHa_u2t1.v1_message.v1_gcrop.nc4")
p <- brick("land_hurtt/LAND_USE_HARMONIZATION_V1_1248/data/LUHa_u2t1.v1_message.v1_gpast.nc4")
u <- brick("land_hurtt/LAND_USE_HARMONIZATION_V1_1248/data/LUHa_u2t1.v1_message.v1_gurbn.nc4")

#calculate total unsuitable land
tu <- overlay(c,p,u,fun=function(c,p,u,na.rm=T){return(c+p+u)},filename="rcp85_unsuitland.tif")

#calculate 10 year averages
tu10<-raster()
for (i in 1:20) {
  index <- c(1:nlayers(tu))
  five <- c(((i*5-4)-2):((i*5-4)+2))
  window <- intersect(index,five)
  win <- dropLayer(tu,c(0:five[1]-1,five[length(five)]+1:nlayers(tu)))
  new10 <- calc(win,fun=mean,na.rm=T)
  tu10 <- stack(tu10,new10)
}
writeRaster(tu10,filename="rcp85_unsuitland10.tif")



####
#extract reef suit vals
#find maximum unsuitable land % within 50km radius of reef site - save max % unsuitability as raster
#classify as suitable or not - save as .csv
###

polyfiles <- list.files("reef_polygon",pattern="\\.shp$",full.names=T)
nfiles <- length(polyfiles)
rcp85suitcsvout <- paste("rcp85_",2000+(1:nlayers(tu10))*5,"_suit.csv",sep="")
tu10 <- brick("rcp85_unsuitland10.tif")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)
for (i in 1:nlayers(tu10)){
  land <- tu10[[i]]
  df.extract1 <- foreach::foreach(i=1:5, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract2 <- foreach::foreach(i=6:10, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract3 <- foreach::foreach(i=11:15, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract4 <- foreach::foreach(i=16:20, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract5 <- foreach::foreach(i=20:25, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract6 <- foreach::foreach(i=26:30, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  dff <- rbind(df.extract1,df.extract2)
  dff1 <- rbind(dff,df.extract3)
  dff2 <- rbind(dff1,df.extract4)
  dff3 <- rbind(dff2,df.extract5)
  dff4 <- rbind(dff3,df.extract6)
  dfext <- as.data.frame(dff4)
  dfext <- subset(dfext,select=-c(grid_code, BUFF_DIST, ORIG_FID,Shape_Leng,Shape_Area,geometry))
  names(dfext) <- c("pointid","zone","maxunsuit")
  result <- transform(dfext,pointid=as.numeric(pointid),maxunsuit=as.numeric(maxunsuit))
  result$suit <- with(result,ifelse(maxunsuit <0.5, 1, 0)) #if the percentage of unsuitable land is less than 0.5, it is suitable (1)
  result$suit <- with(result, ifelse(is.na(maxunsuit),1,suit)) #if land is NA, then land is too far away to return value so is suitable (1)
  b<-merge(IDmatch,result,by="pointid",all.x=T,all.y=F)
  c <- b[!duplicated(b$ID),]
  d <- c[base::order(c$ID),]
  write.csv(d,rcp85suitcsvout[i], row.names=FALSE)
  gc()
}
stopCluster(cl)
finish <- Sys.time()
finish-start

#get count of suitability
csvfiles <- list.files("land_hurtt/rcp85",pattern="_suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  b <- a[!(a$iscoral==0),]
  c <- count(b, vars=suit)
  y <- add_column(c,year=paste0(i*5+2000),.before=T)
  count <- rbind(y,count)
} 
