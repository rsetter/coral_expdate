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


###RCP26

####
#extract reef suit vals
#find dhw for each point
#classify as suitable or not - save as .csv
###

setwd("dhw/rcp26")


reefs <- st_read("reefs.shp")
m <- list.files("dhw/rcp26/", pattern="median",full.names=T)
dhw <- raster()
for(i in 2:length(m)){
  r <- raster(m[i])
  dhw <- stack(dhw,r)
}
rcp26suitcsvout <- paste("rcp26_",2000+(1:20)*5,"_suit.csv",sep="")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

reefext <- rgis::fast_extract(sf=reefs,ras=dhw,funct=function(x){max(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
count<-data.frame()
foreach::foreach(i=1:20, .packages=c("raster")) %dopar% {
  x <- paste0("dhw_modelmedian_bc_5y",2000+i*5)
  dfi <- subset(dfext,select=c("pointid","grid_code", "iscoral","ID",x))
  names(dfi) <- c("pointid","grid_code", "iscoral","ID","dhw")
  dfi$suit <- with(dfi, ifelse(dhw <= 8,1, 0)) #suitable number of dhw is <8
  write.csv(dfi,rcp26suitcsvout[i], row.names=FALSE)
  
  dfi$dhw_coral <- dfi$suit * dfi$iscoral
  c <- count(dfi, vars="dhw_coral")
  y <- add_column(c,year=paste0(i*5+2000),.before=T)
  count <- rbind(y,count)
}

stopCluster(cl)
finish <- Sys.time()
finish-start
#22 min





###RCP45

####
#extract reef suit vals
#find dhw for each point
#classify as suitable or not - save as .csv
###

setwd("dhw/rcp45")


m <- list.files("dhw/rcp45/", pattern="median",full.names=T)
dhw <- raster()
for(i in 2:length(m)){
  r <- raster(m[i])
  dhw <- stack(dhw,r)
}
rcp45suitcsvout <- paste("rcp45_",2000+(1:20)*5,"_suit.csv",sep="")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

reefext <- rgis::fast_extract(sf=reefs,ras=dhw,funct=function(x){max(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
foreach::foreach(i=1:20, .packages=c("raster","plyr","tibble")) %dopar% {
  x <- paste0("dhw_modelmedian_bc_5y",2000+i*5)
  dfi <- subset(dfext,select=c("pointid","grid_code", "iscoral","ID",x))
  names(dfi) <- c("pointid","grid_code", "iscoral","ID","dhw")
  dfi$suit05[dfi$dhw <=8] <- 1
  dfi$suit05[dfi$dhw > 8] <- 0 #suitable number of dhw is <8
  write.csv(dfi,rcp45suitcsvout[i], row.names=FALSE)
  
  dfi$dhw_coral <- dfi$suit * dfi$iscoral
  c <- count(dfi, vars="dhw_coral")
  y <- add_column(c,year=paste0(i*5+2000),.before=T)
  count <- rbind(y,count)
}

stopCluster(cl)
finish <- Sys.time()
finish-start




###RCP85

####
#extract reef suit vals
#find dhw for each point
#classify as suitable or not - save as .csv
###

setwd("dhw/rcp85")


m <- list.files("dhw/rcp85/", pattern="median",full.names=T)
dhw <- raster()
for(i in 2:length(m)){
  r <- raster(m[i])
  dhw <- stack(dhw,r)
}
rcp85suitcsvout <- paste("rcp85_",2000+(1:20)*5,"_suit.csv",sep="")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

reefext <- rgis::fast_extract(sf=reefs,ras=dhw,funct=function(x){max(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
for (i in 1:20){
  x <- paste0("dhw_modelmedian_bc_5y",2000+i*5)
  dfi <- subset(dfext,select=c("pointid","grid_code", "iscoral","ID",x))
  names(dfi) <- c("pointid","grid_code", "iscoral","ID","dhw")
  dfi$suit05[dfi$dhw <=8] <- 1
  dfi$suit05[dfi$dhw > 8] <- 0 #suitable number of dhw is <8
  write.csv(dfi,rcp85suitcsvout[i], row.names=FALSE)
}

stopCluster(cl)
finish <- Sys.time()
finish-start



###historic oisst

####
#extract reef suit vals
#find dhw for each point
#classify as suitable or not - save as .csv
###

setwd("dhw/historic/empirical")


m <- list.files("dhw/emp/", pattern=".tif",full.names=T)
dhw <- raster()
for(i in 1:length(m)){
  r <- raster(m[i])
  dhw <- stack(dhw,r)
}
histsuitcsvout <- paste("hist_",1980+(1:7)*5,"_suit.csv",sep="")



start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

count<-data.frame()

reefext <- rgis::fast_extract(sf=reefs,ras=dhw,funct=function(x){max(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
foreach::foreach(i=1:7, .packages=c("raster","plyr")) %dopar% {
  x <- paste0("dhw_emp_oisst_5y",1980+i*5)
  dfi <- subset(dfext,select=c("pointid","grid_code", "iscoral","ID",x))
  names(dfi) <- c("pointid","grid_code", "iscoral","ID","dhw")
  dfi$suit <- with(dfi, ifelse(dhw <= 8,1, 0)) #suitable number of dhw is <8
  write.csv(dfi,histsuitcsvout[i], row.names=FALSE)
  
  dfi$dhw_coral <- dfi$suit * dfi$iscoral
  c <- count(dfi, vars="dhw_coral")
  y <- add_column(c,year=paste0(i*5+2000),.before=T)
  count <- rbind(y,count)
}

stopCluster(cl)
finish <- Sys.time()
finish-start






#historic model

####
#extract reef suit vals
#find dhw for each point
#classify as suitable or not - save as .csv
###

setwd("dhw/historic")


m <- list.files("dhw/hist/", pattern="median",full.names=T)
dhw <- raster()
for(i in 2:length(m)){
  r <- raster(m[i])
  dhw <- stack(dhw,r)
}
histsuitcsvout <- paste("hist_",1845+(1:16)*10,"_suit.csv",sep="")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

reefext <- rgis::fast_extract(sf=reefs,ras=dhw,funct=function(x){max(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
foreach::foreach(i=1:16, .packages=c("raster","plyr")) %dopar% {
  x <- paste0("dhw_modelmedian_bc_10y",1845+i*10,sep="")
  dfi <- subset(dfext,select=c("pointid","grid_code", "iscoral","ID",x))
  names(dfi) <- c("pointid","grid_code", "iscoral","ID","dhw")
  dfi$suit <- with(dfi, ifelse(dhw <= 8,1, 0)) #suitable number of dhw is <8
  write.csv(dfi,histsuitcsvout[i], row.names=FALSE)
  

}

stopCluster(cl)
finish <- Sys.time()
finish-start

