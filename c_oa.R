library(raster)
library(ncdf4)
library(doParallel)
library(velox)
library(sf)
library(foreach)
library(rgis)
library(dplyr)
library(tibble)
library(beepr)

###empirical

####
#extract reef suit vals
#find omega arag for each point
#classify as suitable or not - save as .csv
###

setwd("oa/empirical")

reefs <- st_read("reefs.shp")
oa <- raster("oa/empirical/OAempiricalfocalcoast.tif")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="")
registerDoParallel(cl)

reefext <- rgis::fast_extract(sf=reefs,ras=oa,funct=function(x){max(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)

  dfi <- subset(dfext,select=c("pointid","grid_code","iscoral","ID","oa_OAempiricalfocalcoast"))
  names(dfi) <- c("pointid","grid_code","iscoral","ID","oa")
  dfi$suit <- with(dfi, ifelse(oa >= 3.3, 1, 0)) #if omega arag is > 3.3, it is suitable(1). otherwise not (0)
  write.csv(dfi,"emp_2005new_suit.csv", row.names=FALSE)

stopCluster(cl)
finish <- Sys.time()
finish-start
beep()


#get count of suitability
b <- dfi[!(dfi$iscoral==0),]
countemp <- count(b, vars="suit")

empcsv <- read.csv(file="oa/empirical/emp_2005new_suit.csv",row.names=NULL,header=T)
  b <- empcsv[!(empcsv$iscoral==0),]
  countemp <- count(b, vars=suit)

beep()






###historic model

####
#extract reef suit vals
#find omega arag for each point
#classify as suitable or not - save as .csv
###

setwd("oa/hist")

reefs <- st_read("currentcoral/reefs.shp")
oa <- brick("oa/co2sys_input/historical/OAhist_modelmedian.tif") #this is model median with all available models
histsuitcsvout <- paste("hist_",1845+(1:16)*10,"_suit.csv",sep="")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

reefext <- rgis::fast_extract(sf=reefs,ras=oa,funct=function(x){max(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
for (i in 1:16){
  x <- paste0("oa_OA_hist_modelmedian.",i)
  dfi <- subset(dfext,select=c("pointid","grid_code", "iscoral","ID",x))
  names(dfi) <- c("pointid","grid_code", "iscoral","ID","oa")
  dfi$suit <- with(dfi, ifelse(oa >= 3.3, 1, 0)) #if omega arag is > 3.3, it is suitable(1). otherwise not (0)
  write.csv(dfi,histsuitcsvout[i], row.names=FALSE)
}

stopCluster(cl)
finish <- Sys.time()
finish-start
beep()








###rcp26

####
#extract reef suit vals
#find omega arag for each point
#classify as suitable or not - save as .csv
###

setwd("oa/rcp26")

reefs <- st_read("reefs.shp")
oa <- brick("oa/co2sys_input/rcp26/OArcp26_modelmedian.tif")
rcp26suitcsvout <- paste("rcp26_",2000+(1:20)*5,"_suit.csv",sep="")

start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

reefext <- rgis::fast_extract(sf=reefs,ras=oa,funct=function(x){max(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
for (i in 1:20){
  x <- paste0("oa_OA_rcp26_modelmedian.",i)
  dfi <- subset(dfext,select=c("pointid","grid_code", "iscoral","ID",x))
  names(dfi) <- c("pointid","grid_code", "iscoral","ID","oa")
  dfi$suit <- with(dfi, ifelse(oa >= 3.3, 1, 0)) #if omega arag is > 3.3, it is suitable(1). otherwise not (0)
  write.csv(dfi,rcp26suitcsvout[i], row.names=FALSE)
}

stopCluster(cl)
finish <- Sys.time()
finish-start
beep()





###rcp45

####
#extract reef suit vals
#find omega arag for each point
#classify as suitable or not - save as .csv
###

setwd("oa/rcp45")

reefs <- st_read("reefs.shp")
oa <- brick("oa/co2sys_input/rcp45/OArcp45_modelmedian.tif")
rcp45suitcsvout <- paste("rcp45_",2000+(1:20)*5,"_suit.csv",sep="")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

reefext <- rgis::fast_extract(sf=reefs,ras=oa,funct=function(x){max(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
for (i in 1:20){
  x <- paste0("oa_OA_rcp45_modelmedian.",i)
  dfi <- subset(dfext,select=c("pointid","grid_code", "iscoral","ID",x))
  names(dfi) <- c("pointid","grid_code", "iscoral","ID","oa")
  dfi$suit <- with(dfi, ifelse(oa >= 3.3, 1, 0)) #if omega arag is > 3.3, it is suitable(1). otherwise not (0)
  write.csv(dfi,rcp45suitcsvout[i], row.names=FALSE)
}

stopCluster(cl)
finish <- Sys.time()
finish-start






###rcp85

####
#extract reef suit vals
#find omega arag for each point
#classify as suitable or not - save as .csv
###

setwd("oa/rcp85")

reefs <- st_read("reefs.shp")
oa <- brick("oa/co2sys_input/rcp85/OArcp85_modelmedian.tif")
rcp85suitcsvout <- paste("rcp85_",2000+(1:20)*5,"_suit.csv",sep="")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

reefext <- rgis::fast_extract(sf=reefs,ras=oa,funct=function(x){max(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
for (i in 1:20){
  x <- paste0("oa_OA_rcp85_modelmedian.",i)
  dfi <- subset(dfext,select=c("pointid","grid_code", "iscoral","ID",x))
  names(dfi) <- c("pointid","grid_code", "iscoral","ID","oa")
  dfi$suit <- with(dfi, ifelse(oa >= 3.3, 1, 0)) #if omega arag is > 3.3, it is suitable(1). otherwise not (0)
  write.csv(dfi,rcp85suitcsvout[i], row.names=FALSE)
}

stopCluster(cl)
finish <- Sys.time()
finish-start


#get count of suitability
csvfiles <- list.files("oa/rcp85",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count85<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  b <- a[!(a$iscoral==0),]
  c <- count(b, vars=suit)
  y <- add_column(c,year=paste0(i*5+2000),.before=T)
  count85 <- rbind(y,count85)
} 
beep()
