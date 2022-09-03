library(raster)
library(ncdf4)
library(doParallel)
library(velox)
library(sf)
library(foreach)
library(rgis)
library(dplyr)
library(tibble)



###
#historic model - hyde 3.1
###

setwd("pop_hr/hist_hyde")

popd1850 <- raster("popd_1850AD.asc")


#calculate 5 year averages
histhydefiles <- list.files("pop_hr/hist_hyde", pattern="popd_", full.names=T)
histhyde5yrnames <- paste("histhyde_",1845+10*(1:15),"_dens.tif",sep="")

for (i in 1:16){
  x <- raster(histhydefiles[i])
  n <- i+1
  y <- raster(histhydefiles[n])
  xy5 <- overlay(x,y,fun=function(x){ mean(x,na.rm=T)},filename=histhyde5yrnames[i])
}




####
#extract reef suit vals
#find maximum population density within 50km radius of reef site - save max pop as raster
#classify as suitable or not - save as .csv
###

polyfiles <- list.files("reef_polygon",pattern="\\.shp$",full.names=T)
nfiles <- length(polyfiles)
reefs <- st_read("currentcoral/reefs.shp")
IDmatch <- read.csv("land_hurtt/polygonIDmatch.csv",header=T)
names(IDmatch) <- c("ObjectID","iscoral","ID","pointid")
hydefiles <- list.files("pop_hr/hist_hyde",pattern="_densbc", full.names=T)
nhydefiles <- length(hydefiles)
hydesuitcsvout <- paste("histhyde_",1845+10*(1:16),"_bc_suit.csv",sep="")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)
for (i in 1:nhydefiles){
  pop <- raster(hydefiles[i])
  df.extract1 <- foreach::foreach(i=1:5, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract2 <- foreach::foreach(i=6:10, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract3 <- foreach::foreach(i=11:15, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract4 <- foreach::foreach(i=16:20, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract5 <- foreach::foreach(i=20:25, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract6 <- foreach::foreach(i=26:30, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  dff <- rbind(df.extract1,df.extract2)
  dff1 <- rbind(dff,df.extract3)
  dff2 <- rbind(dff1,df.extract4)
  dff3 <- rbind(dff2,df.extract5)
  dff4 <- rbind(dff3,df.extract6)
  dfext <- as.data.frame(dff4)
  dfext <- subset(dfext,select=-c(grid_code, BUFF_DIST, ORIG_FID,Shape_Leng,Shape_Area,geometry))
  names(dfext) <- c("pointid","zone","maxpop")
  result <- transform(dfext,pointid=as.numeric(pointid),maxpop=as.numeric(maxpop))
  result$logpop <- with(result, ifelse(maxpop==0 | is.na(maxpop),0,log10(maxpop)))
  result$suit <- with(result, ifelse(logpop <2, 1, 0))
  result$suit <- with(result, ifelse(is.na(logpop),1,suit))
  b<-merge(IDmatch,result,by="pointid",all.x=T,all.y=F)
  c <- b[!duplicated(b$ID),]
  d <- c[base::order(c$ID),]
  write.csv(d,hydesuitcsvout[i], row.names=FALSE)
  gc()
}
stopCluster(cl)
finish <- Sys.time()
finish-start

#get count of suitability
csvfiles <- list.files("pop_hr/hist_hyde",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  b <- a[!(a$iscoral==0),]
  c <- count(b, vars=suit)
  y <- add_column(c,year=paste0(1845+10*i),.before=T)
  count <- rbind(y,count)
} 




###
#historic empirical
###

setwd("pop_hr/hist_gpw")

histgpw <- raster("pop_hr/hist_gpw/gpw_v4_population_density_rev11_2005_30_sec.tif")

####
#extract reef suit vals
#find maximum population density within 50km radius of reef site - save max pop as raster
#classify as suitable or not - save as .csv
###

polyfiles <- list.files("reef_polygon",pattern="\\.shp$",full.names=T)
nfiles <- length(polyfiles)
reefs <- st_read("currentcoral/reefs.shp")
IDmatch <- read.csv("land_hurtt/polygonIDmatch.csv",header=T)
names(IDmatch) <- c("ObjectID","iscoral","ID","pointid")



start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

  pop <- histgpw
  crs(pop) <- crs(reefs)
  df.extract1 <- foreach::foreach(i=1:5, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract2 <- foreach::foreach(i=6:10, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract3 <- foreach::foreach(i=11:15, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract4 <- foreach::foreach(i=16:20, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract5 <- foreach::foreach(i=20:25, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract6 <- foreach::foreach(i=26:30, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  dff <- rbind(df.extract1,df.extract2)
  dff1 <- rbind(dff,df.extract3)
  dff2 <- rbind(dff1,df.extract4)
  dff3 <- rbind(dff2,df.extract5)
  dff4 <- rbind(dff3,df.extract6)
  dfext <- as.data.frame(dff4)
  dfext <- subset(dfext,select=-c(grid_code, BUFF_DIST, ORIG_FID,Shape_Leng,Shape_Area,geometry))
  names(dfext) <- c("pointid","zone","maxpop")
  result <- transform(dfext,pointid=as.numeric(pointid),maxpop=as.numeric(maxpop))
  result$logpop <- with(result, ifelse(maxpop==0 | is.na(maxpop),0,log10(maxpop)))
  result$suit <- with(result, ifelse(logpop <2, 1, 0))
  result$suit <- with(result, ifelse(is.na(logpop),1,suit))
  b<-merge(IDmatch,result,by="pointid",all.x=T,all.y=F)
  c <- b[!duplicated(b$ID),]
  d <- c[base::order(c$ID),]
  write.csv(d,"histgpw_v1_2005_suit.csv", row.names=FALSE)
  gc()

stopCluster(cl)
finish <- Sys.time()
finish-start

#get count of suitability
  suitcor <- d[!(d$iscoral==0),]
  count(suitcor, vars="suit")


#make coarse empirical dataset for comparison against historic model
hyde2005 <- raster("pop_hr/hist_hyde/popd_2005_dens.asc")
gpw4v2005 <- raster("pop_hr/hist_gpw/gpw_v4_population_density_rev11_2005_30_sec.tif")
gpw4v2005a <- aggregate(gpw4v2005,fact=10,fun=mean,na.rm=T,filename="pop_hr/hist_gpw/gpw_v4_population_density_083.tif")




### bias calculation
library(terra)

#open datasets
#open empirical GPW 2000-2015
emp.files <- c("pop_hr/hist_gpw/gpw_v4_population_density_rev11_2000_30_sec.tif",
               "pop_hr/hist_gpw/gpw_v4_population_density_rev11_2005_30_sec.tif",
               "pop_hr/hist_gpw/gpw_v4_population_density_rev11_2010_30_sec.tif",
               "pop_hr/hist_gpw/gpw_v4_population_density_rev11_2015_30_sec.tif")
emppop <- rast(emp.files)
#open historical scenario 2000-2015
hist.files <- c("pop_hr/hist_hyde/popd_2000AD.asc","pop_hr/hist_hyde/popd_2005AD.asc",
                "pop_hr/hist_hyde/popd_2010AD.asc","pop_hr/hist_hyde/popd_2015AD.asc")
histpop <- rast(hist.files)
#open projected scenarios 2000-2015
rcp26.files <- c("pop_hr/base_2000_dens.tif","pop_hr/ssp1/ssp1_2005_dens.tif",
             "pop_hr/ssp1/ssp1_2010_dens.tif","pop_hr/ssp1/ssp1_2015_dens.tif")
rcp45.files <- c("pop_hr/base_2000_dens.tif","pop_hr/ssp1/ssp2_2005_dens.tif",
                 "pop_hr/ssp1/ssp2_2010_dens.tif","pop_hr/ssp1/ssp2_2015_dens.tif")
rcp85.files <- c("pop_hr/base_2000_dens.tif","pop_hr/ssp1/ssp5_2005_dens.tif",
                 "pop_hr/ssp1/ssp5_2010_dens.tif","pop_hr/ssp1/ssp5_2015_dens.tif")
rcp26pop <- rast(rcp26.files)
rcp45pop <- rast(rcp45.files)
rcp85pop <- rast(rcp85.files)


#calculate mean of each dataset
#empirical
empm <- mean(emppop, na.rm=T)
#historical
histm <- mean(histpop, na.rm=T)
#projected
rcp26m <- mean(rcp26pop, na.rm=T)
rcp45m <- mean(rcp45pop, na.rm=T)
rcp85m <- mean(rcp85pop, na.rm=T)

#calculate bias (difference between model and gpw)
#historical
emphist <- sds(empm,histm)
histbias <- lapp(emphist,fun=function(x,y){x-y},filename="pop_hr/hist_hyde/histpopbias.tif")
#projected
emprcp26 <- sds(empm,rcp26m)
emprcp45 <- sds(empm,rcp45m)
emprcp85 <- sds(empm,rcp85m)
rcp26bias <- lapp(emprcp26,fun=function(x,y){x-y},filename="pop_hr/ssp1/ssp1popbias.tif")
rcp45bias <- lapp(emprcp45,fun=function(x,y){x-y},filename="pop_hr/ssp2/ssp2popbias.tif")
rcp85bias <- lapp(emprcp85,fun=function(x,y){x-y},filename="pop_hr/ssp5/ssp5popbias.tif")

#apply bias correction to datasets
year <- c(1845+10*(1:16),2005+5*(1:19))
#historical
histfiles <- list.files("pop_hr/hist_hyde",pattern="_dens.tif$",full.names=T)
for(i in 1:length(histfiles)){
  f <- c("pop_hr/hist_hyde/histpopbias.tif",i)
  r <- rast(f)
  lapp(r,fun=sum,filename=paste0("pop_hr/hist_hyde/hist_popd_",year[i],"_densbc.tif",sep=""))
}
#projected
rcp26files <- list.files("pop_hr/ssp1",pattern="_dens.tif$",full.names=T)
rcp45files <- list.files("pop_hr/ssp2",pattern="_dens.tif$",full.names=T)
rcp85files <- list.files("pop_hr/ssp5",pattern="_dens.tif$",full.names=T)
for(i in 1:length(rcp26files)){
  f <- c("pop_hr/ssp1/ssp1popbias.tif",i)
  r <- rast(f)
  lapp(r,fun=sum,filename=paste0("pop_hr/ssp1/ssp1_popd_",year[i+15],"_densbc.tif",sep=""))
}
for(i in 1:length(rcp45files)){
  f <- c("pop_hr/ssp2/ssp2popbias.tif",i)
  r <- rast(f)
  lapp(r,fun=sum,filename=paste0("pop_hr/ssp2/ssp2_popd_",year[i+15],"_densbc.tif",sep=""))
}
for(i in 1:length(rcp85files)){
  f <- c("pop_hr/ssp5/ssp5popbias.tif",i)
  r <- rast(f)
  lapp(r,fun=sum,filename=paste0("pop_hr/ssp5/ssp5_popd_",year[i+15],"_densbc.tif",sep=""))
}

for(i in 17:length(rcp85files)){
  a<- rast(rcp85files[i])
  b <- rast("pop_hr/ssp5/ssp5popbias.tif")
  ext(b) <- ext(a)
  r <- sds(a,b)
  t <- lapp(r,fun=function(x,y){x+y})
  t[t<0] <-0
  writeRaster(t,filename=paste0("pop_hr/ssp5/ssp5_popd_",year[i+15],"_densbc.tif",sep=""))
}




###
#SSP1
###

setwd("pop_hr/ssp1")

#make area raster
r1.2010 <- raster("human_population/population_projection/ssp1/SSP1_1km/ssp1_total_2010.nc4")
area <- area(reefbuff, filename="area_1km.tif", na.rm=FALSE) #raster layer where Cell values represent the size of the cell in km2

#crop to only the coastal areas near reef sites. then calculate pop density for all years
path.ssp1 <- "human_population/population_projection/ssp1/SSP1_1km/"
ssp1files <- list.files(path.ssp1, pattern="_total_", full.names=T)
nfiles <- length(ssp1files)
ssp1_list  <- paste("ssp1.20",1:nfiles,"0",sep="")
ssp1_dens_list <- paste("ssp1.20",1:nfiles,"0d",sep="")
reefbuff <- raster("pop_hr/reeflandbuffer.tif")
ssp1densfiles <- paste("ssp1_20",1:nfiles,"0_dens.tif",sep="")

for (i in 1:nfiles){
  
  x <- raster(ssp1files[i]) #open each raster
  assign(ssp1_list[i],x)
  
  y <- extend(x,extent(reefbuff),value=NA) #change extent so matches reef buffer 
  extent(y) <- extent(reefbuff)
  z <- overlay(y,reefbuff,fun=function(x,y){(x*y)}) #get rid of inland areas we aren't interested in
  
  d <- overlay(z,area,fun=function(x,y,na.rm=T){return(x/y)},filename=ssp1densfiles[i]) #calculate density
  assign(ssp1_dens_list[i],d)
}

#calculate 5 year averages
ssp1dfiles <- list.files("pop_hr/ssp1", pattern="_dens.tif", full.names=T)
ssp15yrnames <- paste("ssp1_20",1:nfiles,"5_dens.tif",sep="")

for (i in 1:nfiles){
  x <- raster(ssp1dfiles[i])
  n <- i+1
  y <- raster(ssp1dfiles[n])
  xy5 <- overlay(x,y,fun=function(x){ mean(x,na.rm=T)},filename=ssp15yrnames[i])
}


#calculate 2005 pop density using the 2000 base year
#calculate 2000 pop density, then calculate average with 2010
base2000 <- raster("pop_hr/baseYr_total_2000.tif")
area <- raster("pop_hr/area_1km.tif")
reefbuff <- raster("pop_hr/reeflandbuffer.tif")
y <- extend(base2000,extent(reefbuff),value=NA) #change extent so matches reef buffer 
extent(y) <- extent(reefbuff)
z <- overlay(y,reefbuff,fun=function(x,y){(x*y)}) #get rid of inland areas we aren't interested in
dens2000 <- overlay(z,area,fun=function(x,y,na.rm=T){return(x/y)},filename="pop_hr/base_2000_dens.tif") #calculate density
#calculate 2005 density
ssp1dens2010 <- raster("ssp1_2010_dens.tif")
ssp1dens2005 <- overlay(ssp1dens2010,dens2000,fun=function(x){ mean(x,na.rm=T)},filename="ssp1_2005_dens.tif") #calculate 2005 pop dens



####
#extract reef suit vals
#find maximum population density within 50km radius of reef site - save max pop as raster
#classify as suitable or not - save as .csv
###

polyfiles <- list.files("reef_polygon",pattern="\\.shp$",full.names=T)
nfiles <- length(polyfiles)
reefs <- st_read("currentcoral/reefs.shp")
IDmatch <- read.csv("land_hurtt/polygonIDmatch.csv",header=T)
names(IDmatch) <- c("ObjectID","iscoral","ID","pointid")
ssp1files <- list.files("pop_hr/ssp1",pattern="_densbc", full.names=T)
nssp1files <- length(ssp1files)
ssp1suitcsvout <- paste("ssp1_",2000+(1:nssp1files)*5,"_bc_suit.csv",sep="")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)
for (i in 1:nssp1files){
  pop <- raster(ssp1files[i])
  df.extract1 <- foreach::foreach(i=1:5, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract2 <- foreach::foreach(i=6:10, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract3 <- foreach::foreach(i=11:15, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract4 <- foreach::foreach(i=16:20, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract5 <- foreach::foreach(i=20:25, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract6 <- foreach::foreach(i=26:30, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  dff <- rbind(df.extract1,df.extract2)
  dff1 <- rbind(dff,df.extract3)
  dff2 <- rbind(dff1,df.extract4)
  dff3 <- rbind(dff2,df.extract5)
  dff4 <- rbind(dff3,df.extract6)
  dfext <- as.data.frame(dff4)
  dfext <- subset(dfext,select=-c(grid_code, BUFF_DIST, ORIG_FID,Shape_Leng,Shape_Area,geometry))
  names(dfext) <- c("pointid","zone","maxpop")
  result <- transform(dfext,pointid=as.numeric(pointid),maxpop=as.numeric(maxpop))
  result$logpop <- with(result, ifelse(maxpop==0 | is.na(maxpop),0,log10(maxpop)))
  result$suit <- with(result, ifelse(logpop <2, 1, 0))
  result$suit <- with(result, ifelse(is.na(logpop),1,suit))
  b<-merge(IDmatch,result,by="pointid",all.x=T,all.y=F)
  c <- b[!duplicated(b$ID),]
  d <- c[base::order(c$ID),]
  write.csv(d,ssp1suitcsvout[i], row.names=FALSE)
  gc()
}
stopCluster(cl)
finish <- Sys.time()
finish-start

#get count of suitability
csvfiles <- list.files("pop_hr/ssp1",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  b <- a[!(a$iscoral==0),]
  c <- count(b, vars=suit)
  y <- add_column(c,year=paste0(i*5+2000),.before=T)
  count <- rbind(y,count)
} 
 


#####
#####
####
#SSP2
####

setwd("pop_hr/ssp2")

area <- raster("pop_hr/area_1km.tif")

#crop to only the coastal areas near reef sites. then calculate pop density for all years
path.ssp2 <- "human_population/population_projection/ssp2/SSP2_1km/"
ssp2files <- list.files(path.ssp2, pattern="_total_", full.names=T)
nfiles <- length(ssp2files)
ssp2_list  <- paste("ssp2.20",1:nfiles,"0",sep="")
ssp2_dens_list <- paste("ssp2.20",1:nfiles,"0d",sep="")
reefbuff <- raster("pop_hr/reeflandbuffer.tif")
ssp2densfiles <- paste("ssp2_20",1:nfiles,"0_dens.tif",sep="")

for (i in 1:nfiles){
  
  x <- raster(ssp2files[i]) #open each raster
  assign(ssp2_list[i],x)
  
  y <- extend(x,extent(reefbuff),value=NA) #change extent so matches reef buffer 
  extent(y) <- extent(reefbuff)
  z <- overlay(y,reefbuff,fun=function(x,y){(x*y)}) #get rid of inland areas we aren't interested in
  
  d <- overlay(z,area,fun=function(x,y,na.rm=T){return(x/y)},filename=ssp2densfiles[i]) #calculate density
  assign(ssp2_dens_list[i],d)
}

#calculate 5 year averages
ssp2dfiles <- list.files("pop_hr/ssp2", pattern="_dens.tif", full.names=T)
ssp25yrnames <- paste("ssp2_20",1:nfiles,"5_dens.tif",sep="")

for (i in 1:nfiles){
  x <- raster(ssp2dfiles[i])
  n <- i+1
  y <- raster(ssp2dfiles[n])
  xy5 <- overlay(x,y,fun=function(x){ mean(x,na.rm=T)},filename=ssp25yrnames[i])
}

#calculate 2005 density
dens2000 <- raster("pop_hr/base_2000_dens.tif")
ssp2dens2010 <- raster("ssp2_2010_dens.tif")
ssp2dens2005 <- overlay(ssp2dens2010,dens2000,fun=function(x){ mean(x,na.rm=T)},filename="ssp2_2005_dens.tif") #calculate 2005 pop dens


####
#extract reef suit vals
#find maximum population density within 50km radius of reef site - save max pop as raster
#classify as suitable or not - save as .csv
###
setwd("pop_hr/ssp2")

polyfiles <- list.files("reef_polygon",pattern="\\.shp$",full.names=T)
nfiles <- length(polyfiles)
reefs <- st_read("reefs.shp")
IDmatch <- read.csv("land_hurtt/polygonIDmatch.csv",header=T)
names(IDmatch) <- c("ObjectID","iscoral","ID","pointid")
ssp2files <- list.files("pop_hr/ssp2",pattern="_densbc", full.names=T)
nssp2files <- length(ssp2files)
ssp2suitcsvout <- paste("ssp2_",2000+(1:nssp2files)*5,"_bc_suit.csv",sep="")



start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)
for (i in 1:nssp2files){
  pop <- raster(ssp2files[i])
  df.extract1 <- foreach::foreach(i=1:5, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract2 <- foreach::foreach(i=6:10, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract3 <- foreach::foreach(i=11:15, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract4 <- foreach::foreach(i=16:20, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract5 <- foreach::foreach(i=20:25, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract6 <- foreach::foreach(i=26:30, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  dff <- rbind(df.extract1,df.extract2)
  dff1 <- rbind(dff,df.extract3)
  dff2 <- rbind(dff1,df.extract4)
  dff3 <- rbind(dff2,df.extract5)
  dff4 <- rbind(dff3,df.extract6)
  dfext <- as.data.frame(dff4)
  dfext <- subset(dfext,select=-c(grid_code, BUFF_DIST, ORIG_FID,Shape_Leng,Shape_Area,geometry))
  names(dfext) <- c("pointid","zone","maxpop")
  result <- transform(dfext,pointid=as.numeric(pointid),maxpop=as.numeric(maxpop))
  result$logpop <- with(result, ifelse(maxpop==0 | is.na(maxpop),0,log10(maxpop)))
  result$suit <- with(result, ifelse(logpop <2, 1, 0))
  result$suit <- with(result, ifelse(is.na(logpop),1,suit))
  b<-merge(IDmatch,result,by="pointid",all.x=T,all.y=F)
  c <- b[!duplicated(b$ID),]
  d <- c[base::order(c$ID),]
  write.csv(d,ssp2suitcsvout[i], row.names=FALSE)
  gc()
}
stopCluster(cl)
finish <- Sys.time()
finish-start

#get count of suitability
csvfiles <- list.files("pop_hr/ssp2",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  b <- a[!(a$iscoral==0),]
  c <- count(b, vars=suit)
  y <- add_column(c,year=paste0(i*5+2000),.before=T)
  count <- rbind(y,count)
} 



#####
#####
####
#SSP5
####

setwd("pop_hr/ssp5")

area <- raster("pop_hr/area_1km.tif")

#crop to only the coastal areas near reef sites. then calculate pop density for all years
path.ssp5 <- "human_population/population_projection/ssp5/SSP5_1km/"
ssp5files <- list.files(path.ssp5, pattern="_total_", full.names=T)
nfiles <- length(ssp5files)
ssp5_list  <- paste("ssp5.20",1:nfiles,"0",sep="")
ssp5_dens_list <- paste("ssp5.20",1:nfiles,"0d",sep="")
reefbuff <- raster("pop_hr/reeflandbuffer.tif")
ssp5densfiles <- paste("ssp5_20",1:nfiles,"0_dens.tif",sep="")

for (i in 1:nfiles){
  
  x <- raster(ssp5files[i]) #open each raster
  assign(ssp5_list[i],x)
  
  y <- extend(x,extent(reefbuff),value=NA) #change extent so matches reef buffer 
  extent(y) <- extent(reefbuff)
  z <- overlay(y,reefbuff,fun=function(x,y){(x*y)}) #get rid of inland areas we aren't interested in
  
  d <- overlay(z,area,fun=function(x,y,na.rm=T){return(x/y)},filename=ssp5densfiles[i]) #calculate density
  assign(ssp5_dens_list[i],d)
}

#calculate 5 year averages
ssp5dfiles <- list.files("pop_hr/ssp5", pattern="_dens.tif", full.names=T)
ssp55yrnames <- paste("ssp5_20",1:nfiles,"5_dens.tif",sep="")

for (i in 1:nfiles){
  x <- raster(ssp5dfiles[i])
  n <- i+1
  y <- raster(ssp5dfiles[n])
  xy5 <- overlay(x,y,fun=function(x){ mean(x,na.rm=T)},filename=ssp55yrnames[i])
}


#calculate 2005 density
ssp5dens2010 <- raster("ssp5_2010_dens.tif")
ssp5dens2005 <- overlay(ssp5dens2010,dens2000,fun=function(x){ mean(x,na.rm=T)},filename="ssp5_2005_dens.tif") #calculate 2005 pop dens



####
#extract reef suit vals
#find maximum population density within 50km radius of reef site - save max pop as raster
#classify as suitable or not - save as .csv
###
setwd("pop_hr/ssp5")

polyfiles <- list.files("reef_polygon",pattern="\\.shp$",full.names=T)
nfiles <- length(polyfiles)
reefs <- st_read("currentcoral/reefs.shp")
IDmatch <- read.csv("land_hurtt/polygonIDmatch.csv",header=T)
names(IDmatch) <- c("ObjectID","iscoral","ID","pointid")
ssp5files <- list.files("pop_hr/ssp5",pattern="_densbc", full.names=T)
nssp5files <- length(ssp5files)
ssp5suitcsvout <- paste("ssp5_",2000+(1:nssp5files)*5,"_bc_suit.csv",sep="")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)
for (i in 1:nssp5files){
  pop <- raster(ssp5files[i])
  df.extract1 <- foreach::foreach(i=1:5, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract2 <- foreach::foreach(i=6:10, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract3 <- foreach::foreach(i=11:15, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract4 <- foreach::foreach(i=16:20, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract5 <- foreach::foreach(i=20:25, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract6 <- foreach::foreach(i=26:30, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  dff <- rbind(df.extract1,df.extract2)
  dff1 <- rbind(dff,df.extract3)
  dff2 <- rbind(dff1,df.extract4)
  dff3 <- rbind(dff2,df.extract5)
  dff4 <- rbind(dff3,df.extract6)
  dfext <- as.data.frame(dff4)
  dfext <- subset(dfext,select=-c(grid_code, BUFF_DIST, ORIG_FID,Shape_Leng,Shape_Area,geometry))
  names(dfext) <- c("pointid","zone","maxpop")
  result <- transform(dfext,pointid=as.numeric(pointid),maxpop=as.numeric(maxpop))
  result$logpop <- with(result, ifelse(maxpop==0 | is.na(maxpop),0,log10(maxpop)))
  result$suit <- with(result, ifelse(logpop <2, 1, 0))
  result$suit <- with(result, ifelse(is.na(logpop),1,suit))
  b<-merge(IDmatch,result,by="pointid",all.x=T,all.y=F)
  c <- b[!duplicated(b$ID),]
  d <- c[base::order(c$ID),]
  write.csv(d,ssp5suitcsvout[i], row.names=FALSE)
  gc()
}
stopCluster(cl)
finish <- Sys.time()
finish-start


#get count of suitability
csvfiles <- list.files("pop_hr/ssp5",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  b <- a[!(a$iscoral==0),]
  c <- count(b, vars=suit)
  y <- add_column(c,year=paste0(i*5+2000),.before=T)
  count <- rbind(y,count)
}  




