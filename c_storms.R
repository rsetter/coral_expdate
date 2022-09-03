library(raster)
library(ncdf4)
library(doParallel)
library(velox)
library(sf)
library(foreach)
library(rgis)
library(dplyr)
library(tibble)
library(plyr)
library(plotrix)
library(tmap)
library(rgdal)
library(fasterize)
library(rasterVis)
library(rworldmap)




setwd("")

#######
#use empirical IBTrACS dataset to calibrate the power dissipation datasets
########

#import polygon storm tracks (storm tracks with buffer of maximum wind speed radius)
#rasterize wind speed value in knots. sum overlapping areas.
stormpolys <- list.files(path="ibtracs_years/buffer",pattern='*.shp$',full.names=T)
stormlines <- list.files(path="ibtracs_years",pattern='*.shp$',full.names=T)
ibtracout <- paste("ibtracs_",1984+(1:21),"max.tif",sep="")
y <- raster("histmodel_annual.tif")
for (i in 1:length(stormpolys)) {
  poly <- st_read(stormpolys[i])
  line <- st_read(stormlines[i+137])
  a<-rasterize(poly,y,field="USA_WIND",fun=max)
  b<-rasterize(line,y,field="USA_WIND",fun=max)
  merge(b,a,filename=ibtracout[i],overwrite=T)
}
#compile all rasters into single brick file
rastfiles <- list.files(path="empirical",pattern='*max.tif$',full.names=T)
ibtracsall <- lapply(rastfiles, raster::raster)
empibtracs <- brick(ibtracsall)
writeRaster(empibtracs,filename="ibtracs_1985_2005max.tif",overwrite=T)

#create convex hull to only analyze sites within coral points. add a 100km buffer
reefs <- st_read("coralpoint.shp")
cohull <- st_convex_hull(st_union(reefs))
qtm(cohull)
cohull100 <- st_buffer(cohull,dist=1)
st_write(cohull100,dsn="convexhullreef.shp")
chull <- readOGR("convexhullreef.shp")
empibtracsc <- mask(ibtracsmax,chull)

#plot ibtracs & hist model mean 1985-2005
ibtracsmax <- brick("ibtracs_1985_2005max.tif")
ibtracsmax[is.na(ibtracsmax)] <- 0
ibtracsmaxm<-calc(ibtracsmax,fun=mean,na.rm=T)
ibtracsmaxm[is.na(ibtracsmaxm)] <- 0
p <- levelplot(ibtracsmaxm,main='IBTrACS Mean Windspeed 1985-2005', par.settings = viridisTheme)
p + layer(sp.lines(coastsCoarse, lwd=0.8, col='darkgray'))

modpdi <- brick("storms/historic/histmodel_annual.tif")
modelpdi <- dropLayer(modpdi,c(1:35)) #we only want to compare the same temporal range - 1985-2005
modelmean <- calc(modelpdi,fun=mean,na.rm=T)
q <- levelplot(modelmean, main='Historical Model Mean PDI 1985-2005', par.settings = viridisTheme)
data("coastsCoarse")
q + layer(sp.lines(coastsCoarse,lwd=0.8,col='darkgray'))

#define a severe storm using quantile calibration method
empibtracsr<-reclassify(ibtracsmax,c(-1,113,0,
                                     113,5000,1)) #severe storms are those with >113 knots
empibtracsrm<-calc(empibtracsr,fun=mean,na.rm=T)
p <- levelplot(empibtracsrm, main='ibtracs severe storms', par.settings = viridisTheme)
p + layer(sp.lines(coastsCoarse, lwd=0.8, col='darkgray'))
ibtracsmaxc <- mask(ibtracsmax,chull)
ibtracsmaxd<-dropLayer(ibtracsmaxc,c(1:10))
vals <- getValues(ibtracsmaxd)
vals[is.na(vals)] <- 0 #make all NA values 0
hist(vals)
quant<-quantile(vals, na.rm=T)
qu <- ecdf(vals)
qu(113) #percentile for category 4 storms is 0.9958365 


modpdi <- brick("storms/historic/histmodel_annual.tif")
modelpdi <- dropLayer(modpdi,c(1:35)) #we only want to compare the same temporal range - 1985-2005
modelpdic <- mask(modelpdi,chull)
valsm <- getValues(modelpdic)
valsm[is.na(valsm)] <- 0
hist(valsm)
quantile(valsm,na.rm=T,c(0.9958365)) #category 4 is same quantile as 911389592 pdi
modpdir <- reclassify(modpdi,c(-1,911389592,0, 
                               911389592,382752358400,1)) #severe storms are those with >911389592 pdi
modpdirm<-calc(modpdir,fun=mean,na.rm=T)
q <- levelplot(modpdirm, main='historical model severe storms', par.settings = viridisTheme)
q + layer(sp.lines(coastsCoarse, lwd=0.8, col='darkgray'))






###
#historic empirical
###

setwd("storms/empirical")

ibtracsmax <- brick("ibtracs_1985_2005max.tif")
ibtracsmax[is.na(ibtracsmax)] <- 0

#calculate the number of years between repeat storms occurring in each cell in a ten year window (5 years before and 5 years after)
#threshold: severe storms are those with windspeed > 113. sites with return years less than 5 years are unsuitable
#function find the minimum number of years between severe storms in a given pixel
fn <- function(x,na.rm=T) {
  v<-which(x>113)
  ifelse(1>=length(v),return(NA),return(min(abs(diff(v)))))
}
returnrast <- paste("ibtracs_",1:nlayers(empibtracs)+1984,"_maxreturnyr.tif",sep="")
for (i in 1:nlayers(ibtracsmax)) {
  index <- c(1:nlayers(ibtracsmax))
  ten <- c((i-5):(i+5))
  window <- intersect(index,ten)
  empw <- dropLayer(ibtracsmax,c(0:ten[1],ten[length(ten)]+1:nlayers(ibtracsmax)))
  empwr <- calc(empw,fun=fn,na.rm=T,filename=returnrast[i],overwrite=T)
}
empyfiles <- list.files(path="storms/empirical",pattern='*maxreturnyr.tif$',full.names=T)
empyr <- lapply(empyfiles, raster::raster)
empyb <- brick(empyr)
writeRaster(empyb,filename="ibtracs_all_maxreturnyr.tif",overwrite=T)


#make 5 year raster
empy5files <- empyfiles[seq(1,length(empyfiles),5)]
emp5yr <- lapply(empy5files, raster::raster)
emp5yb <- brick(emp5yr)
writeRaster(emp5yb,filename="ibtracs_all_return5yr.tif")



####
#extract reef suit vals
#find number of years between severe storms for each point
#classify as suitable or not - save as .csv
###


reefs <- st_read("reefs.shp")
storm <- brick("ibtracs_all_return5yr.tif")
empsuitcsvout <- paste("emp_",1980+(1:5)*5,"_suit.csv",sep="")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

reefext <- rgis::fast_extract(sf=reefs,ras=storm,funct=function(x){min(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
for (i in 1:5){
  x <- paste0("storm_ibtracs_all_return5yr.",i)
  dfi <- subset(dfext,select=c("pointid","grid_code", "iscoral","ID",x))
  names(dfi) <- c("pointid","grid_code", "iscoral","ID", "storm")
  dfi$suit <- with(dfi, ifelse(storm <= 5,0, 1)) #if there are less than 5 years between severe storms, site is unsuitable
  dfi$suit <- with(dfi, ifelse(is.na(storm),1,suit)) #if no windspeed is recorded for the pixel, then it's suitable
  write.csv(dfi,empsuitcsvout[i], row.names=FALSE)
}
stopCluster(cl)
finish <- Sys.time()
finish-start


#get count of suitability
csvfiles <- list.files(path="storms/empirical",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  a$storm_coral <- a$suit * a$iscoral
  c <- count(a, vars="storm_coral")
  y <- add_column(c,year=paste0(i*5+1980),.before=T)
  count <- rbind(y,count)
} 








#historical model


setwd("storms/historic")
histmdf <- read.csv(file="storms/projected/storm_multimodelmean_rcp26.csv", header=TRUE, sep=",")
histmr <- rasterFromXYZ(histmdf,crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ")
histmr <- dropLayer(histmr,c(57:151))
histmr <- rotate(histmr,filename="storms/historic/histmodel_annual.tif",overwrite=T)



#calculate the number of years between repeat storms occurring in each cell in a ten year window (5 years before and 5 years after)
#threshold: severe storms are those with pdi > 911389592 (category 4 and up) sites with return years less than 5 years are unsuitable
#function find the minimum number of years between severe storms in a given pixel
fn <- function(x,na.rm=T) {
  v<-which(x>911389592)
  ifelse(1>=length(v),return(NA),return(min(abs(diff(v)))))
}
returnrast <- paste("histmod_",1:nlayers(histmr)+1949,"_returnyr.tif",sep="")
for (i in seq(from=1, to=(nlayers(histmr)),by=5)) {
  index <- c(1:nlayers(histmr))
  ten <- c((i-5):(i+5))
  window <- intersect(index,ten)
  modw <- dropLayer(histmr,c(0:ten[1],ten[length(ten)]+1:nlayers(histmr)))
  modwr <- calc(modw,fun=fn,na.rm=T,filename=returnrast[i],overwrite=T) 
}
modyfiles <- list.files(path="storms/historic",pattern='*_returnyr.tif',full.names=T)
modyr <- lapply(modyfiles, raster::raster)
modyb <- brick(modyr)
writeRaster(modyb,filename="histmod_all_returnyr.tif",overwrite=T)
q <- levelplot(modyb[[46]], main='historic model return year - 1995', par.settings = RdBuTheme)
q + layer(sp.lines(coastsCoarse, lwd=0.8, col='darkgray'))


####
#extract reef suit vals
#find number of years between severe storms for each point
#classify as suitable or not - save as .csv
###


reefs <- st_read("reefs.shp")
storm <- brick("histmod_all_returnyr.tif")
suitcsvout <- paste("hist_",1945+(1:12)*5,"_suit.csv",sep="")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

reefext <- rgis::fast_extract(sf=reefs,ras=storm,funct=function(x){min(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
for (i in 1:nlayers(storm)){
  x <- paste0("storm_histmod_all_returnyr.",i)
  dfi <- subset(dfext,select=c("pointid","grid_code", "iscoral","ID",x))
  names(dfi) <- c("pointid","grid_code", "iscoral","ID", "storm")
  dfi$suit <- with(dfi, ifelse(storm <= 5,0, 1)) #if there are less than 5 years between severe storms, site is unsuitable
  dfi$suit <- with(dfi, ifelse(is.na(storm),1,suit)) #if no windspeed is recorded for the pixel, then it's suitable
  #stormpt <- inner_join(reef,dfi,by="pointid")
  #stormdf <- as.data.frame(stormpt)
  write.csv(dfi,suitcsvout[i], row.names=FALSE)
  #rasterize(dhwpt,r008333,field="dhw",fun=max,na.rm=T,filename=rcp85rout[i],overwrite=T)
  #rasterize(dhwpt,r008333,field="suit",fun=max,na.rm=T,filename=rcp85suitout[i],overwrite=T)
}
stopCluster(cl)
finish <- Sys.time()
finish-start


#get count of suitability
csvfiles <- list.files(path="storms/historic",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  a$storm_coral <- a$suit * a$iscoral
  c <- count(a, vars="storm_coral")
  y <- add_column(c,year=paste0(i*5+1945),.before=T)
  count <- rbind(y,count)
} 




###RCP26

####
#extract reef suit vals
#find storm wind for each point
#classify as suitable or not - save as .csv
###

setwd("storms/rcp26")


storm26df <- read.csv(file="storms/projected/storm_multimodelmean_rcp26.csv", header=TRUE, sep=",")
storm26xyz <- rasterFromXYZ(storm26df,crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ")
storm26 <- dropLayer(storm26xyz,c(1:50))
storm26r <- rotate(storm26,filename="rcp26_annual.tif",overwrite=T)

#calculate the number of years between repeat storms occurring in each cell in a ten year window (5 years before and 5 years after)
#threshold: severe storms are those with pdi > 911389592 (category 4 and up) sites with return years less than 5 years are unsuitable
#function find the minimum number of years between severe storms in a given pixel
fn <- function(x,na.rm=T) {
  v<-which(x>911389592)
  ifelse(1>=length(v),return(NA),return(min(abs(diff(v)))))
}
returnrast <- paste("storms/rcp26/rcp26_",1:(nlayers(storm26r)-5)+2004,"_returnyr.tif",sep="")
for (i in seq(from=1, to=(nlayers(storm26r)-5),by=5)) {
  index <- c(1:nlayers(storm26r))
  ten <- c((i-5):(i+5))
  window <- intersect(index,ten)
  modw <- dropLayer(storm26r,c(0:ten[1],ten[length(ten)]+1:nlayers(storm26r)))
  modwr <- calc(modw,fun=fn,na.rm=T,filename=returnrast[i],overwrite=T) 
}
modyfiles <- list.files(path="storms/rcp26",pattern='*_returnyr.tif',full.names=T)
modyr <- lapply(modyfiles, raster::raster)
modyb <- brick(modyr)
writeRaster(modyb,filename="rcp26_returnyr.tif",overwrite=T)



####
#extract reef suit vals
#find number of years between severe storms for each point
#classify as suitable or not - save as .csv
###


reefs <- st_read("currentcoral/reefs.shp")
storm <- brick("rcp26_returnyr.tif")
suitcsvout <- paste("rcp26_",2000+(1:20)*5,"_suit.csv",sep="")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

reefext <- rgis::fast_extract(sf=reefs,ras=storm,funct=function(x){min(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
for (i in 1:nlayers(storm)){
  x <- paste0("storm_rcp26_returnyr.",i)
  dfi <- subset(dfext,select=c("pointid","grid_code", "iscoral","ID",x))
  names(dfi) <- c("pointid","grid_code", "iscoral","ID", "storm")
  dfi$suit <- with(dfi, ifelse(storm <= 5,0, 1)) #if there are less than 5 years between severe storms, site is unsuitable
  dfi$suit <- with(dfi, ifelse(is.na(storm),1,suit)) #if no windspeed is recorded for the pixel, then it's suitable
  write.csv(dfi,suitcsvout[i], row.names=FALSE)
}
stopCluster(cl)
finish <- Sys.time()
finish-start


#get count of suitability
csvfiles <- list.files(path="storms/rcp26",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  a$storm_coral <- a$suit * a$iscoral
  c <- count(a, vars="storm_coral")
  y <- add_column(c,year=paste0(i*5+2000),.before=T)
  count <- rbind(y,count)
} 










###RCP45

####
#extract reef suit vals
#find storm wind for each point
#classify as suitable or not - save as .csv
###

setwd("storms/rcp45")

storm45df <- read.csv(file="storms/projected/storm_multimodelmean_rcp60.csv", header=TRUE, sep=",")
storm45xyz <- rasterFromXYZ(storm45df,crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ")
storm45 <- dropLayer(storm45xyz,c(1:50))
storm45r <- rotate(storm45,filename="rcp45_annual.tif",overwrite=T)

#calculate the number of years between repeat storms occurring in each cell in a ten year window (5 years before and 5 years after)
#threshold: severe storms are those with pdi > 911389592 (category 4 and up) sites with return years less than 5 years are unsuitable
#function find the minimum number of years between severe storms in a given pixel
fn <- function(x,na.rm=T) {
  v<-which(x>911389592)
  ifelse(1>=length(v),return(NA),return(min(abs(diff(v)))))
}
returnrast <- paste("storms/rcp45/rcp45_",1:(nlayers(storm45r)-5)+2004,"_returnyr.tif",sep="")
for (i in seq(from=1, to=(nlayers(storm45r)-5),by=5)) {
  index <- c(1:nlayers(storm45r))
  ten <- c((i-5):(i+5))
  window <- intersect(index,ten)
  modw <- dropLayer(storm45r,c(0:ten[1],ten[length(ten)]+1:nlayers(storm45r)))
  modwr <- calc(modw,fun=fn,na.rm=T,filename=returnrast[i],overwrite=T) 
}
modyfiles <- list.files(path="storms/rcp45",pattern='*_returnyr.tif',full.names=T)
modyr <- lapply(modyfiles, raster::raster)
modyb <- brick(modyr)
writeRaster(modyb,filename="rcp45_returnyr.tif",overwrite=T)



####
#extract reef suit vals
#find number of years between severe storms for each point
#classify as suitable or not - save as .csv
###


reefs <- st_read("currentcoral/reefs.shp")
storm <- brick("rcp45_returnyr.tif")
suitcsvout <- paste("rcp45_",2000+(1:20)*5,"_suit.csv",sep="")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

reefext <- rgis::fast_extract(sf=reefs,ras=storm,funct=function(x){min(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
for (i in 1:nlayers(storm)){
  x <- paste0("storm_rcp45_returnyr.",i)
  dfi <- subset(dfext,select=c("pointid","grid_code", "iscoral","ID",x))
  names(dfi) <- c("pointid","grid_code", "iscoral","ID", "storm")
  dfi$suit <- with(dfi, ifelse(storm <= 5,0, 1)) #if there are less than 5 years between severe storms, site is unsuitable
  dfi$suit <- with(dfi, ifelse(is.na(storm),1,suit)) #if no windspeed is recorded for the pixel, then it's suitable
  write.csv(dfi,suitcsvout[i], row.names=FALSE)
}
stopCluster(cl)
finish <- Sys.time()
finish-start


#get count of suitability
csvfiles <- list.files(path="storms/rcp45",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  a$storm_coral <- a$suit * a$iscoral
  c <- count(a, vars="storm_coral")
  y <- add_column(c,year=paste0(i*5+2000),.before=T)
  count <- rbind(y,count)
} 








###RCP85

####
#extract reef suit vals
#find storm wind for each point
#classify as suitable or not - save as .csv
###

setwd("storms/rcp85")


reefs <- st_read("reefs.shp")
storm85df <- read.csv(file="storms/projected/storm_multimodelmean_rcp85.csv", header=TRUE, sep=",")
storm85xyz <- rasterFromXYZ(storm85df,crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ")
storm85 <- dropLayer(storm85xyz,c(1:50))
storm85r <- rotate(storm85,filename="rcp85_annual.tif",overwrite=T)

#calculate the number of years between repeat storms occurring in each cell in a ten year window (5 years before and 5 years after)
#threshold: severe storms are those with pdi > 911389592 (category 4 and up) sites with return years less than 5 years are unsuitable
#function find the minimum number of years between severe storms in a given pixel
fn <- function(x,na.rm=T) {
  v<-which(x>911389592)
  ifelse(1>=length(v),return(NA),return(min(abs(diff(v)))))
}
returnrast <- paste("storms/rcp85/rcp85_",1:(nlayers(storm85r)-5)+2004,"_returnyr.tif",sep="")
for (i in seq(from=1, to=(nlayers(storm85r)-5),by=5)) {
  index <- c(1:nlayers(storm85r))
  ten <- c((i-5):(i+5))
  window <- intersect(index,ten)
  modw <- dropLayer(storm85r,c(0:ten[1],ten[length(ten)]+1:nlayers(storm85r)))
  modwr <- calc(modw,fun=fn,na.rm=T,filename=returnrast[i],overwrite=T) 
}
modyfiles <- list.files(path="storms/rcp85",pattern='*_returnyr.tif',full.names=T)
modyr <- lapply(modyfiles, raster::raster)
modyb <- brick(modyr)
writeRaster(modyb,filename="rcp85_returnyr.tif",overwrite=T)



####
#extract reef suit vals
#find number of years between severe storms for each point
#classify as suitable or not - save as .csv
###


reefs <- st_read("reefs.shp")
storm <- brick("rcp85_returnyr.tif")
suitcsvout <- paste("rcp85_",2000+(1:20)*5,"_suit.csv",sep="")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

reefext <- rgis::fast_extract(sf=reefs,ras=storm,funct=function(x){min(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
for (i in 1:nlayers(storm)){
  x <- paste0("storm_rcp85_returnyr.",i)
  dfi <- subset(dfext,select=c("pointid","grid_code", "iscoral","ID",x))
  names(dfi) <- c("pointid","grid_code", "iscoral","ID", "storm")
  dfi$suit <- with(dfi, ifelse(storm <= 5,0, 1)) #if there are less than 5 years between severe storms, site is unsuitable
  dfi$suit <- with(dfi, ifelse(is.na(storm),1,suit)) #if no windspeed is recorded for the pixel, then it's suitable
  write.csv(dfi,suitcsvout[i], row.names=FALSE)
}
stopCluster(cl)
finish <- Sys.time()
finish-start


#get count of suitability
csvfiles <- list.files(path="storms/rcp85",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  a$storm_coral <- a$suit * a$iscoral
  c <- count(a, vars="storm_coral")
  y <- add_column(c,year=paste0(i*5+2000),.before=T)
  count <- rbind(y,count)
} 


