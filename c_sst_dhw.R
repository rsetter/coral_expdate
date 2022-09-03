library(raster)
library(ncdf4)
library(doParallel)
library(velox)
library(sf)
library(foreach)


# Utility functions to track how long functions are taking
# Additionally sets process cores and enables parallelism
startTimer = function() {
  start = Sys.time()
  cores = detectCores()-1
  cl = makeCluster(cores, output="")
  registerDoParallel(cl)
  return(start)
}

# Utility function which takes a startTime and returns function duration
stopTimer = function(start) {
  stopCluster(cl)
  finish = Sys.time()
  return(finish - start)
}

# Utility function to process models - monthly bias
processModelsmonth = function(models, biasPath, dhwPath, sstPath, max=95, year=2005) {
  for(i in 1:length(models)) {

    # Get all 12 bias rasters for each month
    for(j in 1:12){
      object <- paste("bias.", j, sep = "")
      raster <- raster(paste0(biasPath, modelnames[i], "_bias", j, ".tif", sep=""))
      assign(a,r)
    }

    biasmonths <- lapply(paste0('bias.',1:12),get)

    # Get model into comparable raster format
    i_brickC <- calc(models[[i]], fun=function(x){x - 273.15}) #change from K to C
    i_brickCd <- disaggregate(i_brickC,fact=2)
    i_brick1 <- crop(i_brickCd,extent(180,359.5,-90,90))
    i_brick2 <- crop(i_brickCd,extent(-0.5,180,-90,90))
    extent(i_brick1) <- c(-180,-0.5,-90,90)
    i_brickCc <- merge(i_brick1,i_brick2)

    foreach::foreach(l=1:max, .packages=c("raster","foreach")) %dopar% {
      l1 <- 12*(l-1)+1
      l2 <- 12*l

      SSTbc <- foreach::foreach(k=l1:l2, .packages=c("raster"), .combine=raster::stack) %dopar% {
          i_brickf <- focal(i_brickCc[[k]],w=matrix(1,nrow=3,ncol=3),fun=mean,na.rm=T,NAonly=T)
          raster::overlay(i_brickf,biasmonths[[monthlabel[k]]],fun=function(x,y,na.rm=T){return((x+y))})
      }

      DHWx <- foreach::foreach(h=1:12, .packages=c("raster"), .combine=raster::stack) %dopar% {
          DHW <- raster::overlay(SSTbc[[h]],mmm,fun=function(x,y,na.rm=T){return((x-y)*4.34)})
          reclassify(DHW, c(-Inf,1,0))
      }

      filename.dhw <- paste0(dhwPath, modelnames[i], "/", modelnames[i], "_", year+l, ".tif", sep="")
      filename.sst <- paste0(sstPath, modelnames[i], "/", modelnames[i], "_bc_", year+l, ".tif", sep="")

      calc(DHWx,fun=sum,filename=filename.dhw)
      calc(SSTbc,fun=mean,filename=filename.sst)

      # Collect garbage - save the planet
      gc()
    }
  }
}

# Utility function to process models - overall bias
processModels = function(models, modelnames, biasPath, dhwPath, sstPath, maxyr, yearstart,MMM) {
  for(i in 1:length(models)) {
    
    # Get model bias raster
    bias <- raster(paste0(biasPath,modelnames[i],"_bias.tif"))
    
    # Get model into comparable raster format
    i_brickC <- calc(models[[i]], fun=function(x){x - 273.15}) #change from K to C
    i_brickCd <- disaggregate(i_brickC,fact=2)
    i_brick1 <- crop(i_brickCd,extent(180,359.5,-90,90))
    i_brick2 <- crop(i_brickCd,extent(-0.5,180,-90,90))
    extent(i_brick1) <- c(-180,-0.5,-90,90)
    i_brickCc <- merge(i_brick1,i_brick2)
    
    #apply bias
    SSTbc <- overlay(i_brickCc,bias,fun=sum)
    
    foreach::foreach(l=1:maxyr, .packages=c("raster","foreach")) %dopar% {
      l1 <- 12*(l-1)+1
      l2 <- 12*l
      
      SSTyr <- stack(SSTbc[[l1:l2]])
      
      DHWx <- foreach::foreach(h=l1:l2, .packages=c("raster"), .combine=raster::stack) %dopar% {
        sstf <- focal(SSTbc[[h]],w=matrix(1,nrow=3,ncol=3),fun=mean,na.rm=T,NAonly=T)
        DHW <- raster::overlay(sstf,MMM,fun=function(x,y,na.rm=T){return((x-y)*4.34)})
        reclassify(DHW, c(-Inf,1,0))
      }
      
      filename.dhw <- paste0(dhwPath, modelnames[i], "/", modelnames[i], "_", yearstart+l, ".tif", sep="")
      filename.sst <- paste0(sstPath, modelnames[i], "/", modelnames[i], "_bc_", yearstart+l, ".tif", sep="")
      
      calc(DHWx,fun=sum,filename=filename.dhw,overwrite=T)
      calc(SSTyr,fun=mean,filename=filename.sst, overwrite=T)
      
      # Collect garbage - save the planet
      gc()
    }
  }
}





#### SST BIAS CALCULATION



#open oisst. only keep 1982-2005
#use oisst monthly
oisst <- brick("oisstmasked.nc")
oisst <- rotate(oisst)
oisst<- dropLayer(oisst,c(1,290:450))

#break it up by months. get mean monthly values
jans <- c(1+12*(0:23))

for(j in 1:12){
  months <- stack(oisst[[jans[1:24]+j-1]])
  calc(months,fun=mean,na.rm=T,filename=paste0("bias/oisst_mean_",j,".tif"))
}

#calculate mean sst from entire 1982-2005 timespan
meanoisst <- calc(oisst,fun=mean,na.rm=T,filename=paste0("bias/oisst_mean.tif"))
meanoisstd <- disaggregate(meanoisst,fact=2)
meanoisstdf <- focal(meanoisstd,w=matrix(1,nrow=3,ncol=3),fun=mean,na.rm=T,NAonly=T)

#calculate MMM - maximum monthly mean from oisst dataset (using 1982-1992)
oisstmmm <- dropLayer(oisst,c(133:288))
jans <- c(1+12*(0:10))
maxmm <- raster()
for(j in 1:12){
  months <- stack(oisstmmm[[jans[1:10]+j-1]])
  newmmm <- calc(months,fun=max,na.rm=T)
  maxmm <- stack(maxmm,newmmm)
}
mmm <- calc(maxmm,fun=max,na.rm=T)
mmmf <- focal(mmm,w=matrix(1,nrow=3,ncol=3),fun=mean,na.rm=T,NAonly=T,filename="oisst_mmm8292f.tif")
mmmfd <- disaggregate(mmmf,fact=2)



#open historic models.

setwd("sst/hist")

parent.folder <- "SST/CMIP5/historical"
# Get list of file names matching the .nc pattern
m <- list.files(parent.folder, pattern=".nc",full.names=T)
# Remove files from given indices
n <- m[-c(1,6,7,12,19,20,21,23,24)]
# Create mh.[N] objects assigned to corresponding raster
for(i in 1:length(n)){
  object <- paste("mh.", i, sep = "")
  r <- brick(n[i])
  assign(object,r)
}

# adjust rasters so have same temporal range
# clean up datasets
e <- raster(nrow=180,ncol=360,ext=extent(mh.1),crs=crs(mh.1))
e[]<-NA
# too many layers. remove the last 60 layers
mh.5<-dropLayer(mh.5,c(1873:1932))
# too many layers. remove the last 60 layers
mh.7<-dropLayer(mh.7,c(1873:1932))
e120<-stack(replicate(120,e))
# add 120
mh.9<-stack(e120,mh.9)
# add 120
mh.10<-stack(e120,mh.10)
e119<-stack(replicate(119,e))
mh.11<-stack(e119,mh.11)
mh.15[mh.15<270] <- NA
mh.5[mh.5<270] <- NA

# Get the mh.1 to mh:N models
models <- lapply(paste0('mh.',1:length(n)),get)
modelnames <- c("CanESM2","CMCC-CESM","CMCC-CM","CMCC-CMS","GISS-E2-H-CC","GISS-E2-H","GISS-E2-R-CC","GISS-E2-R", 
                "HadGEM2-AO","HadGEM2-CC","HadGEM2-ES","inmcm4","MIROC-ESM-CHEM","MIROC-ESM","MRI-CGCM3")

#calculate mean overall value for 1982-2005
#calculate bias (difference between model and oisst)
startTime <- startTimer()
for(i in 1:length(models)){
  model <- models[[i]]
  model[model < 270] <- NA
  meanmodel <- calc(model,fun=mean,na.rm=T)
  modelC <- calc(meanmodel, fun=function(x){x - 273.15}) #change from K to C
  modelCd <- disaggregate(modelC,fact=2)
  model1 <- crop(modelCd,extent(180,359.5,-90,90))
  model2 <- crop(modelCd,extent(-0.5,180,-90,90))
  extent(model1) <- c(-180,-0.5,-90,90)
  modelm <- merge(model1,model2)
  modelmf <- focal(modelm,w=matrix(1,nrow=3,ncol=3),fun=mean,na.rm=T,NAonly=T)
  bias <- raster::overlay(meanoisstdf,modelmf,fun=function(x,y){return((x-y))},filename=paste0("sst/bias/",modelnames[i],"_bias.tif"))
}
stopTimer(startTime)


#adjust bias to models.
#calculate DHW
#calculate mean year SST. save
bias.folder <- "sst/bias/"
memory.limit(size=50000000) 

startTime <- startTimer()
processModels(models, modelnames, biasPath=bias.folder, dhwPath="dhw/hist/", sstPath="sst/hist/", maxyr=156, yearstart=1849,MMM=mmmfd)
stopTimer(startTime)
#30 min per model
#total time = 6.78 hours


#calculate model mean and model median
startTime <- startTimer()
foreach::foreach(i=1:156, .packages=c("raster")) %dopar% {
  dhwb <- raster()
  sstb <- raster()
  for(j in 1:length(modelnames)){
    dhwr <- raster(paste0("dhw/hist/",modelnames[j],"/",modelnames[j],"_",1849+i,".tif",sep=""))
    dhwb <- stack(dhwb,dhwr)
    sstr <- raster(paste0("sst/hist/",modelnames[j],"/",modelnames[j],"_bc_",1849+i,".tif",sep=""))
    sstb <- stack(sstb,sstr)
  }
  calc(dhwb,fun=mean,na.rm=T,filename=paste0("dhw/hist/modelmean/modelmean_bc_",1849+i,".tif",sep=""),overwrite=T)
  calc(dhwb,fun=median,na.rm=T,filename=paste0("dhw/hist/modelmedian/modelmedian_bc_",1849+i,".tif",sep=""),overwrite=T)
  calc(sstb,fun=mean,na.rm=T,filename=paste0("sst/hist/modelmean/modelmean_bc_",1849+i,".tif",sep=""),overwrite=T)
  calc(sstb,fun=median,na.rm=T,filename=paste0("sst/hist/modelmedian/modelmedian_bc_",1849+i,".tif",sep=""),overwrite=T)
}
stopTimer(startTime)
#total time = 5 min

#calculate 10 year mean
startTime <- startTimer()
foreach::foreach(i=1:16, .packages=c("raster")) %dopar% {
  one = i*10-5
  ten = i*10+5
  dhwmb <- raster()
  dhwdb <- raster()
  sstmb <- raster()
  sstdb <- raster()
  for(j in one:ten){
    dhwmr <- raster(paste0("dhw/hist/modelmean/modelmean_bc_",1845+j,".tif",sep=""))
    dhwmb <- stack(dhwmb,dhwmr)
    dhwdr <- raster(paste0("dhw/hist/modelmedian/modelmedian_bc_",1845+j,".tif",sep=""))
    dhwdb <- stack(dhwdr,dhwdb)
    sstmr <- raster(paste0("sst/hist/modelmean/modelmean_bc_",1845+j,".tif",sep=""))
    sstmb <- stack(sstmr,sstmb)
    sstdr <- raster(paste0("sst/hist/modelmedian/modelmedian_bc_",1845+j,".tif",sep=""))
    sstdb <- stack(sstdr,sstdb)
  }
  calc(dhwmb,fun=mean,filename=paste0("dhw/hist/modelmean_bc_10y",1845+i*10,".tif"),overwrite=T)
  calc(dhwdb,fun=mean,filename=paste0("dhw/hist/modelmedian_bc_10y",1845+i*10,".tif"),overwrite=T)
  calc(sstmb,fun=mean,filename=paste0("sst/hist/modelmean_bc_10y",1845+i*10,".tif"),overwrite=T)
  calc(sstdb,fun=mean,filename=paste0("sst/hist/modelmedian_bc_10y",1845+i*10,".tif"),overwrite=T)
}
stopTimer(startTime)
#total time = 20 sec





##rcp85

#open rcp85 models.

setwd("sst/rcp85")

parent.folder <- "SST/CMIP5/RCP85/grid"
m <- list.files(parent.folder, pattern=".nc",full.names=T)
n <- m[-c(1,2,7,8,9,10,11,12,20,21,22,25,27)]
for(i in 1:length(n)){
  a <- paste("m85.", i, sep = "")
  r <- brick(n[i])
  assign(a,r)
}

#adjust rasters so have same temporal range
#clean up datasets
m85.10 <- dropLayer(m85.10, c(1141))

models <- lapply(paste0('m85.',1:length(n)),get)
modelnames <- c("CanESM2","CMCC-CESM","CMCC-CM","CMCC-CMS","GISS-E2-H-CC","GISS-E2-H","GISS-E2-R-CC","GISS-E2-R",
                "HadGEM2-AO","HadGEM2-CC","inmcm4","MIROC-ESM-CHEM","MIROC-ESM","MRI-CGCM3")





#adjust bias to models.
#calculate DHW
#calculate mean year SST. save
bias.folder <- "sst/bias/"
memory.limit(size=50000000)

startTime <- startTimer()
processModels(models, modelnames, biasPath=bias.folder, dhwPath="dhw/rcp85/", sstPath="sst/rcp85/", maxyr=95, yearstart=2005,MMM=mmmfd)
stopTimer(startTime)
#total time = 3.7 hr

#calculate model mean and model median
startTime <- startTimer()
foreach::foreach(i=1:95, .packages=c("raster")) %dopar% {
  dhwb <- raster()
  sstb <- raster()
  for(j in 1:length(modelnames)){
    dhwr <- raster(paste0("dhw/rcp85/",modelnames[j],"/",modelnames[j],"_",2005+i,".tif",sep=""))
    dhwb <- stack(dhwb,dhwr)
    sstr <- raster(paste0("sst/rcp85/",modelnames[j],"/",modelnames[j],"_bc_",2005+i,".tif",sep=""))
    sstb <- stack(sstb,sstr)
  }
  calc(dhwb,fun=mean,na.rm=T,filename=paste0("dhw/rcp85/modelmean/modelmean_bc_",2005+i,".tif",sep=""),overwrite=T)
  calc(dhwb,fun=median,na.rm=T,filename=paste0("dhw/rcp85/modelmedian/modelmedian_bc_",2005+i,".tif",sep=""),overwrite=T)
  calc(sstb,fun=mean,na.rm=T,filename=paste0("sst/rcp85/modelmean/modelmean_bc_",2005+i,".tif",sep=""),overwrite=T)
  calc(sstb,fun=median,na.rm=T,filename=paste0("sst/rcp85/modelmedian/modelmedian_bc_",2005+i,".tif",sep=""),overwrite=T)
}
stopTimer(startTime)
#total time = 3 min

#calculate 10 year mean
startTime <- startTimer()
foreach::foreach(i=1:20, .packages=c("raster")) %dopar% {
  one = i*5-2
  five = i*5+2
  dhwmb <- raster()
  dhwdb <- raster()
  sstmb <- raster()
  sstdb <- raster()
  for(j in one:five){
    dhwmr <- raster(paste0("dhw/rcp85/modelmean/modelmean_bc_",2000+j,".tif",sep=""))
    dhwmb <- stack(dhwmb,dhwmr)
    dhwdr <- raster(paste0("dhw/rcp85/modelmedian/modelmedian_bc_",2000+j,".tif",sep=""))
    dhwdb <- stack(dhwdr,dhwdb)
    sstmr <- raster(paste0("sst/rcp85/modelmean/modelmean_bc_",2000+j,".tif",sep=""))
    sstmb <- stack(sstmr,sstmb)
    sstdr <- raster(paste0("sst/rcp85/modelmedian/modelmedian_bc_",2000+j,".tif",sep=""))
    sstdb <- stack(sstdr,sstdb)
  }
  calc(dhwmb,fun=mean,filename=paste0("dhw/rcp85/modelmean_bc_5y",2000+i*5,".tif"),overwrite=T)
  calc(dhwdb,fun=mean,filename=paste0("dhw/rcp85/modelmedian_bc_5y",2000+i*5,".tif"),overwrite=T)
  calc(sstmb,fun=mean,filename=paste0("sst/rcp85/modelmean_bc_5y",2000+i*5,".tif"),overwrite=T)
  calc(sstdb,fun=mean,filename=paste0("sst/rcp85/modelmedian_bc_5y",2000+i*5,".tif"),overwrite=T)
}
stopTimer(startTime)
#total time = 1.98 min


##rcp45

#open rcp45 models.

setwd("sst/rcp45")

parent.folder <- "SST/CMIP5/RCP45/grid"
m <- list.files(parent.folder, pattern=".nc",full.names=T)
n <- m[-c(1,2,6:12,21:23,26)]
for(i in 1:length(n)){
  a <- paste("m45.", i, sep = "")
  r <- brick(n[i])
  assign(a,r)
}

#adjust rasters so have same temporal range
#clean up datasets

models <- lapply(paste0('m45.',1:length(n)),get)
modelnames <- c("CanESM2","CMCC-CM","CMCC-CMS","GISS-E2-H-CC","GISS-E2-H","GISS-E2-R-CC","GISS-E2-R",
                "HadGEM2-AO","HadGEM2-CC","HadGEM2-ES","inmcm4","MIROC-ESM-CHEM","MIROC-ESM","MRI-CGCM3")






#adjust bias to models.
#calculate DHW
#calculate mean year SST. save
bias.folder <- "sst/bias/"
memory.limit(size=50000000) 

startTime <- startTimer()
processModels(models, modelnames, biasPath=bias.folder, dhwPath="dhw/rcp45/", sstPath="sst/rcp45/", maxyr=95, yearstart=2005,MMM=mmmfd)
stopTimer(startTime)
#total time = 2.98 hr



#calculate model mean and model median
startTime <- startTimer()
foreach::foreach(i=1:95, .packages=c("raster")) %dopar% {
  dhwb <- raster()
  sstb <- raster()
  for(j in 1:length(modelnames)){
    dhwr <- raster(paste0("dhw/rcp45/",modelnames[j],"/",modelnames[j],"_",2005+i,".tif",sep=""))
    dhwb <- stack(dhwb,dhwr)
    sstr <- raster(paste0("sst/rcp45/",modelnames[j],"/",modelnames[j],"_bc_",2005+i,".tif",sep=""))
    sstb <- stack(sstb,sstr)
  }
  calc(dhwb,fun=mean,na.rm=T,filename=paste0("dhw/rcp45/modelmean/modelmean_bc_",2005+i,".tif",sep=""),overwrite=T)
  calc(dhwb,fun=median,na.rm=T,filename=paste0("dhw/rcp45/modelmedian/modelmedian_bc_",2005+i,".tif",sep=""),overwrite=T)
  calc(sstb,fun=mean,na.rm=T,filename=paste0("sst/rcp45/modelmean/modelmean_bc_",2005+i,".tif",sep=""),overwrite=T)
  calc(sstb,fun=median,na.rm=T,filename=paste0("sst/rcp45/modelmedian/modelmedian_bc_",2005+i,".tif",sep=""),overwrite=T)
}
stopTimer(startTime)
#total time = 1.895363 hours

#calculate 5 year mean
startTime <- startTimer()
foreach::foreach(i=1:20, .packages=c("raster")) %dopar% {
  one <- i*5-2
  five <- i*5+2
  dhwmb <- raster()
  dhwdb <- raster()
  sstmb <- raster()
  sstdb <- raster()
  for(j in one:five){
    dhwmr <- raster(paste0("dhw/rcp45/modelmean/modelmean_bc_",2000+j,".tif",sep=""))
    dhwmb <- stack(dhwmb,dhwmr)
    dhwdr <- raster(paste0("dhw/rcp45/modelmedian/modelmedian_bc_",2000+j,".tif",sep=""))
    dhwdb <- stack(dhwdr,dhwdb)
    sstmr <- raster(paste0("sst/rcp45/modelmean/modelmean_bc_",2000+j,".tif",sep=""))
    sstmb <- stack(sstmr,sstmb)
    sstdr <- raster(paste0("sst/rcp45/modelmedian/modelmedian_bc_",2000+j,".tif",sep=""))
    sstdb <- stack(sstdr,sstdb)
  }
  calc(dhwmb,fun=mean,filename=paste0("dhw/rcp45/modelmean_bc_5y",2000+i*5,".tif"),overwrite=T)
  calc(dhwdb,fun=mean,filename=paste0("dhw/rcp45/modelmedian_bc_5y",2000+i*5,".tif"),overwrite=T)
  calc(sstmb,fun=mean,filename=paste0("sst/rcp45/modelmean_bc_5y",2000+i*5,".tif"),overwrite=T)
  calc(sstdb,fun=mean,filename=paste0("sst/rcp45/modelmedian_bc_5y",2000+i*5,".tif"),overwrite=T)
}
stopTimer(startTime)
#total time = 2 min




##rcp26

#open rcp26 models.

setwd("sst/rcp26")

parent.folder <- "SST/CMIP5/RCP26/grid"
m <- list.files(parent.folder, pattern=".nc",full.names=T)
n <- m[-c(2:6,10,11,12,15)]
for(i in 1:length(n)){
  a <- paste("m26.", i, sep = "")
  r <- brick(n[i])
  assign(a,r)
}


models <- lapply(paste0('m26.',1:length(n)),get)
modelnames <- c("CanESM2","GISS-E2-H","GISS-E2-R","HadGEM2-AO","MIROC-ESM-CHEM","MIROC-ESM","MRI-CGCM3")


#adjust bias to models.
#calculate DHW
#calculate mean year SST. save
bias.folder <- "sst/bias/"
memory.limit(size=50000000) 

startTime <- startTimer()
processModels(models, modelnames, biasPath=bias.folder, dhwPath="dhw/rcp26/", sstPath="sst/rcp26/", maxyr=95, yearstart=2005,MMM=mmmfd)
stopTimer(startTime)


#calculate model mean and model median
startTime <- startTimer()
foreach::foreach(i=1:95, .packages=c("raster")) %dopar% {
  dhwb <- raster()
  sstb <- raster()
  for(j in 1:length(modelnames)){
    dhwr <- raster(paste0("dhw/rcp26/",modelnames[j],"/",modelnames[j],"_",2005+i,".tif",sep=""))
    dhwb <- stack(dhwb,dhwr)
    sstr <- raster(paste0("sst/rcp26/",modelnames[j],"/",modelnames[j],"_bc_",2005+i,".tif",sep=""))
    sstb <- stack(sstb,sstr)
  }
  calc(dhwb,fun=mean,na.rm=T,filename=paste0("dhw/rcp26/modelmean/modelmean_bc_",2005+i,".tif",sep=""),overwrite=T)
  calc(dhwb,fun=median,na.rm=T,filename=paste0("dhw/rcp26/modelmedian/modelmedian_bc_",2005+i,".tif",sep=""),overwrite=T)
  calc(sstb,fun=mean,na.rm=T,filename=paste0("sst/rcp26/modelmean/modelmean_bc_",2005+i,".tif",sep=""),overwrite=T)
  calc(sstb,fun=median,na.rm=T,filename=paste0("sst/rcp26/modelmedian/modelmedian_bc_",2005+i,".tif",sep=""),overwrite=T)
}
stopTimer(startTime)
#total time = 2 hr

#calculate 5 year mean
startTime <- startTimer()
foreach::foreach(i=1:20, .packages=c("raster")) %dopar% {
  one = i*5-2
  five = i*5+2
  dhwmb <- raster()
  dhwdb <- raster()
  sstmb <- raster()
  sstdb <- raster()
  for(j in one:five){
    dhwmr <- raster(paste0("dhw/rcp26/modelmean/modelmean_bc_",2000+j,".tif",sep=""))
    dhwmb <- stack(dhwmb,dhwmr)
    dhwdr <- raster(paste0("dhw/rcp26/modelmedian/modelmedian_bc_",2000+j,".tif",sep=""))
    dhwdb <- stack(dhwdr,dhwdb)
    sstmr <- raster(paste0("sst/rcp26/modelmean/modelmean_bc_",2000+j,".tif",sep=""))
    sstmb <- stack(sstmr,sstmb)
    sstdr <- raster(paste0("sst/rcp26/modelmedian/modelmedian_bc_",2000+j,".tif",sep=""))
    sstdb <- stack(sstdr,sstdb)
  }
  calc(dhwmb,fun=mean,filename=paste0("dhw/rcp26/modelmean_bc_5y",2000+i*5,".tif"),overwrite=T)
  calc(dhwdb,fun=mean,filename=paste0("dhw/rcp26/modelmedian_bc_5y",2000+i*5,".tif"),overwrite=T)
  calc(sstmb,fun=mean,filename=paste0("sst/rcp26/modelmean_bc_5y",2000+i*5,".tif"),overwrite=T)
  calc(sstdb,fun=mean,filename=paste0("sst/rcp26/modelmedian_bc_5y",2000+i*5,".tif"),overwrite=T)
}
stopTimer(startTime)
#total time = 2 min




####
###empirical

setwd("dhw/emp")


#use oisst monthly
oisst <- brick("SST/OISST/oisstmasked.nc")
oisst <- rotate(oisst)
oisst <- dropLayer(oisst,c(1,446:450))

#calculate accumulated DHW/yr. save
monthlabel <- c(rep(1:12,37))
#memory free = 113524476
memory.limit(size=50000000)

startTime <- startTimer()
foreach::foreach(l=1:37, .packages=c("raster","foreach")) %dopar% {
  l1 <- 12*(l-1)+1
  l2 <- 12*l
  DHWx <- foreach::foreach(h=l1:l2, .packages=c("raster"), .combine=raster::stack) %dopar% {
    DHW <- raster::overlay(oisst[[h]],mmmf,fun=function(x,y,na.rm=T){return((x-y)*4.34)})
    reclassify(DHW, c(-Inf,1,0))
  }
  calc(DHWx,fun=sum,filename=paste0("dhw/emp/annual/emp_dhw_",1981+l,".tif",sep=""))
  gc()
}

stopTimer(startTime)
#6.4 min


#calculate 5 year mean
startTime <- startTimer()

foreach::foreach(i=1:7, .packages=c("raster")) %dopar% {
  one = i*5-2
  five = i*5+2
  dhwmb <- raster()
  for(j in one:five){
    dhwmr <- raster(paste0("dhw/emp/annual/emp_dhw_",1980+j,".tif",sep=""))
    dhwmb <- stack(dhwmb,dhwmr)
  }
  calc(dhwmb,fun=mean,filename=paste0("dhw/emp/emp_oisst_5y",1980+i*5,".tif",sep=""))
}

stopTimer(startTime)
