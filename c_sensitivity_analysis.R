library(raster)
library(gdalUtils)
library(gdal)
library(foreach)
library(doParallel)
library(terra) 
library(sf)
library(data.table)



setwd("sensitivity/")

#open coral points shapefile. keep just coral points, not rocky reefs
reefs <- st_read("reefs.shp")
reefcoral <- reefs[reefs$iscoral == 1, ]


#function: open rcp45 raster for variable. 
#extract values for each point, store that in a csv 
variablevaldhw = function(wdvar,datawd) {
  setwd(wdvar)
  
  m <- list.files(datawd, pattern="median",full.names=T)
  var <- raster()
  for(i in 2:length(m)){
    r <- raster(m[i])
    var <- stack(var,r)
  }
  
  start <- Sys.time()
  cores<- detectCores()-1
  cl <- makeCluster(cores, output="") 
  registerDoParallel(cl)
  
  reefext <- rgis::fast_extract(sf=reefcoral,ras=var,funct=function(x){max(x,na.rm=T)},parallel=T) #3.5 min
  
  dfext <- as.data.frame(reefext) #0.3 sec
  names(dfext) <- c("pointid","grid_code", "iscoral","ID","X2005","X2010","X2015","X2020","X2025","X2030","X2035","X2040","X2045","X2050","X2055","X2060","X2065",
                    "X2070","X2075","X2080","X2085","X2090","X2095","X2100","geometry")
  dfval <- dfext[c("ID","X2005","X2010","X2015","X2020","X2025","X2030","X2035","X2040","X2045","X2050","X2055","X2060","X2065",
                   "X2070","X2075","X2080","X2085","X2090","X2095","X2100")]
  
  write.csv(dfval,"dhwrcp45values.csv") #9 min
  
  
  stopCluster(cl)
  finish <- Sys.time()
  finish-start
} #15 min

variablevaloa = function(wdvar,datawd) {
  setwd(wdvar)
  
  m <- list.files(datawd, pattern="median",full.names=T)
  var <- brick(m[1])
  
  start <- Sys.time()
  cores<- detectCores()-1
  cl <- makeCluster(cores, output="") 
  registerDoParallel(cl)
  
  reefext <- rgis::fast_extract(sf=reefcoral,ras=var,funct=function(x){max(x,na.rm=T)},parallel=T) #3.5 min
  
  dfext <- as.data.frame(reefext) #0.3 sec
  names(dfext) <- c("pointid","grid_code", "iscoral","ID","X2005","X2010","X2015","X2020","X2025","X2030","X2035","X2040","X2045","X2050","X2055","X2060","X2065",
                    "X2070","X2075","X2080","X2085","X2090","X2095","X2100","geometry")
  dfval <- dfext[c("ID","X2005","X2010","X2015","X2020","X2025","X2030","X2035","X2040","X2045","X2050","X2055","X2060","X2065",
                   "X2070","X2075","X2080","X2085","X2090","X2095","X2100")]
  
  write.csv(dfval,"oarcp45values.csv") #
  
  
  stopCluster(cl)
  finish <- Sys.time()
  finish-start
} #13 min

variablevalstorm = function(wdvar,datawd) {
  setwd(wdvar)
  
  var <- brick("storms/rcp45/rcp45_returnyr.tif")
  
  start <- Sys.time()
  cores<- detectCores()-1
  cl <- makeCluster(cores, output="") 
  registerDoParallel(cl)
  
  reefext <- rgis::fast_extract(sf=reefcoral,ras=var,funct=function(x){min(x,na.rm=T)},parallel=T) #3.5 min
  
  dfext <- as.data.frame(reefext) #0.3 sec
  names(dfext) <- c("pointid","grid_code", "iscoral","ID","X2005","X2010","X2015","X2020","X2025","X2030","X2035","X2040","X2045","X2050","X2055","X2060","X2065",
                    "X2070","X2075","X2080","X2085","X2090","X2095","X2100","geometry")
  dfval <- dfext[c("ID","X2005","X2010","X2015","X2020","X2025","X2030","X2035","X2040","X2045","X2050","X2055","X2060","X2065",
                   "X2070","X2075","X2080","X2085","X2090","X2095","X2100")]
  
  write.csv(dfval,"stormrcp45values.csv") #
  
  
  stopCluster(cl)
  finish <- Sys.time()
  finish-start
} #8 min

variablevalpop = function(wdvar,datawd) {
  setwd(wdvar)
  
  polyfiles <- list.files("reef_polygon",pattern="\\.shp$",full.names=T)
  nfiles <- length(polyfiles)
  var <- read.csv("land_hurtt/polygonIDmatch.csv",header=T)
  names(var) <- c("ObjectID","iscoral","ID","pointid")
  ssp2files <- list.files(datawd,pattern="_dens", full.names=T)
  nssp2files <- length(ssp2files)
  csvout <- paste("ssp2_",2000+(1:nssp2files)*5,"_val.csv",sep="")
  
  start <- Sys.time()
  cores<- detectCores()-1
  cl <- makeCluster(cores, output="") 
  registerDoParallel(cl)
  for (i in 1:nssp2files){
    pop <- raster(ssp2files[i])
    df.extract1 <- foreach::foreach(j=1:5, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
      poly <- st_read(polyfiles[j])
      reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
    }
    gc()
    df.extract2 <- foreach::foreach(j=6:10, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
      poly <- st_read(polyfiles[j])
      reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
    }
    gc()
    df.extract3 <- foreach::foreach(j=11:15, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
      poly <- st_read(polyfiles[j])
      reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
    }
    gc()
    df.extract4 <- foreach::foreach(j=16:20, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
      poly <- st_read(polyfiles[j])
      reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
    }
    gc()
    df.extract5 <- foreach::foreach(j=20:25, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
      poly <- st_read(polyfiles[j])
      reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
    }
    gc()
    df.extract6 <- foreach::foreach(j=26:30, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
      poly <- st_read(polyfiles[j])
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
    names(dfext) <- c("pointid","zone",paste0(i*5+2000))
    dfext <- subset(dfext,select=-c(zone))
    result <- transform(dfext,pointid=as.numeric(pointid))
    b<-merge(var,result,by="pointid",all.x=T,all.y=F)
    c <- b[!duplicated(b$ID),]
    d <- c[base::order(c$ID),]
    dfval <- d[!(d$iscoral==0),]
    write.csv(dfval,csvout[i], row.names=FALSE)
    
    gc()
  }
  
  stopCluster(cl)
  finish <- Sys.time()
  finish-start
  
} #~1hr/year. 20 hours total. may crash, so adjust i range for whatever years remain 

variablevalland = function(wdvar,datawd) {
  setwd(wdvar)
  
  polyfiles <- list.files("reef_polygon",pattern="\\.shp$",full.names=T)
  nfiles <- length(polyfiles)
  var <- read.csv("land_hurtt/polygonIDmatch.csv",header=T)
  names(var) <- c("ObjectID","iscoral","ID","pointid")
  tu10 <- brick("land_hurtt/rcp45/rcp45_unsuitland10.tif")
  csvout <- paste("ssp2_",2000+(1:nlayers(tu10))*5,"_val.csv",sep="")
  
  start <- Sys.time()
  cores<- detectCores()-1
  cl <- makeCluster(cores, output="") 
  registerDoParallel(cl)
  for (i in 1:nlayers(tu10)){
    land <- tu10[[i]]
    df.extract1 <- foreach::foreach(j=1:5, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
      poly <- st_read(polyfiles[j])
      reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
    }
    gc()
    df.extract2 <- foreach::foreach(j=6:10, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
      poly <- st_read(polyfiles[j])
      reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
    }
    gc()
    df.extract3 <- foreach::foreach(j=11:15, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
      poly <- st_read(polyfiles[j])
      reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
    }
    gc()
    df.extract4 <- foreach::foreach(j=16:20, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
      poly <- st_read(polyfiles[j])
      reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
    }
    gc()
    df.extract5 <- foreach::foreach(j=20:25, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
      poly <- st_read(polyfiles[j])
      reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
    }
    gc()
    df.extract6 <- foreach::foreach(j=26:30, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
      poly <- st_read(polyfiles[j])
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
    names(dfext) <- c("pointid","zone",paste0(i*5+2000))
    dfext <- subset(dfext,select=-c(zone))
    result <- transform(dfext,pointid=as.numeric(pointid))
    b<-merge(var,result,by="pointid",all.x=T,all.y=F)
    c <- b[!duplicated(b$ID),]
    d <- c[base::order(c$ID),]
    dfval <- d[!(d$iscoral==0),]
    write.csv(dfval,csvout[i], row.names=FALSE)
    
    gc()
  }

  stopCluster(cl)
  finish <- Sys.time()
  finish-start
} 
# ~20 min per year, 6.5 hours total (run on HPC - 1 node, 24 cores, 256G RAM) <- on HPC


variablevaldhw("dhw/rcp45","dhw/rcp45")
variablevaloa("oa/rcp45","oa/co2sys_input/rcp45")
variablevalstorm("storm/rcp45","storms/rcp45")
variablevalpop("pop/ssp2","pop_hr/ssp2")
variablevalland("land/ssp2","land_hurtt/rcp45")


#combine pop values to one csv
lf <- list.files("pop/ssp2", full.names=T)
popvaldf <- read.csv(lf[1],header=T)
popvaldf <- subset(popvaldf,select=-c(pointid, ObjectID, iscoral))
for(i in 2:length(lf)){
  t <- read.csv(lf[i],header=T)
  t <- subset(t,select=-c(pointid, ObjectID, iscoral))
  popvaldf <- merge(popvaldf, t, by="ID", all.x=T, all.y=F)
}
write.csv(popvaldf,"pop/ssp2/poprcp45values.csv") 

#combine land values to one csv
lf <- list.files("land/ssp2", full.names=T)
landvaldf <- read.csv(lf[1],header=T)
landvaldf <- subset(landvaldf,select=-c(pointid, ObjectID, iscoral))
for(i in 2:length(lf)){
  t <- read.csv(lf[i],header=T)
  t <- subset(t,select=-c(pointid, ObjectID, iscoral))
  landvaldf <- merge(landvaldf, t, by="ID", all.x=T, all.y=F)
}
write.csv(landvaldf,"land/ssp2/landrcp45values.csv") 






#calculate median expiration date for new threshold

#function: given variable and new threshold value
newthreshold <- function(threshold,stressor,newthreshstressordir,stressordir1,stressordir2,stressordir3,stressordir4){
  start <- Sys.time()
  
  #load stressor values of interest
  lf <- list.files(newthreshstressordir, pattern="rcp45values", full.names=T)
  newstressval <- fread(lf[1],header=T)
  newstressSuit <- subset(newstressval,select=c(ID))
  namecol <- names(newstressval)
  
  #apply threshold: make a dataframe of suitability and nonsuitability (1's and 0's)
  for(i in 3:ncol(newstressval)){
    
    if (stressor == "dhw") {
      newstressSuit$x <- ifelse(newstressval[,..i] <= threshold,1, 0)
      names(newstressSuit)[names(newstressSuit) == 'x'] <- namecol[i]
    } else if (stressor == "oa") {
      newstressSuit$x <- ifelse(newstressval[,..i] >= threshold, 1, 0)
      names(newstressSuit)[names(newstressSuit) == 'x'] <- namecol[i]
    } else if (stressor == "storms") {
      newstressSuit$x <- ifelse(newstressval[,..i] <= threshold,0, 1) 
      newstressSuit$x <- with(newstressSuit, ifelse(is.na(x),1,x)) 
      names(newstressSuit)[names(newstressSuit) == 'x'] <- namecol[i]
    } else if (stressor == "pop"){
      newstressval[newstressval[[i]]==0] <- NA
      newstressval$logpop <- log10(newstressval[,..i])
      newstressSuit$x <- with(newstressval, ifelse(logpop <threshold, 1, 0))
      newstressSuit$x <- with(newstressSuit, ifelse(is.na(x),1,x))
      names(newstressSuit)[names(newstressSuit) == 'x'] <- namecol[i]
    } else if (stressor == "land"){
      newstressSuit$x <- ifelse(newstressval[,..i]<threshold,1,0)
      names(newstressSuit)[names(newstressSuit) == 'x'] <- namecol[i]
    } else {
      print("stressor name misspelled")
    }
  }
  
  #adjust expiration date so only when it expires forever
  cores<- detectCores()-1
  cl <- makeCluster(cores, output="")
  registerDoParallel(cl)
  for (j in (ncol(newstressSuit)-1):2){
    column <- namecol[j]
    future <- namecol[j+1]
    newstressSuit[[column]][newstressSuit[[future]]>newstressSuit[[column]]] <- 1
  } 
  
  stopCluster(cl)
  
  #find expiration year for each site under new stressor threhsold
  newstressSuit$year2005 <- ifelse(newstressSuit$X2005 == 0, 2005,1)
  newstressSuit$year2010 <- ifelse(newstressSuit$year2005 == 1, ifelse(newstressSuit$X2010 ==0,2010,1),newstressSuit$year2005)
  newstressSuit$year2015 <- ifelse(newstressSuit$year2010 == 1, ifelse(newstressSuit$X2015 ==0,2015,1),newstressSuit$year2010)
  newstressSuit$year2020 <- ifelse(newstressSuit$year2015 == 1, ifelse(newstressSuit$X2020 ==0,2020,1),newstressSuit$year2015)
  newstressSuit$year2025 <- ifelse(newstressSuit$year2020 == 1, ifelse(newstressSuit$X2025==0,2025,1),newstressSuit$year2020)
  newstressSuit$year2030 <- ifelse(newstressSuit$year2025 == 1, ifelse(newstressSuit$X2030==0,2030,1),newstressSuit$year2025)
  newstressSuit$year2035 <- ifelse(newstressSuit$year2030 == 1, ifelse(newstressSuit$X2035==0,2035,1),newstressSuit$year2030)
  newstressSuit$year2040 <- ifelse(newstressSuit$year2035 == 1, ifelse(newstressSuit$X2040==0,2040,1),newstressSuit$year2035)
  newstressSuit$year2045 <- ifelse(newstressSuit$year2040 == 1, ifelse(newstressSuit$X2045==0,2045,1),newstressSuit$year2040)
  newstressSuit$year2050 <- ifelse(newstressSuit$year2045 == 1, ifelse(newstressSuit$X2050==0,2050,1),newstressSuit$year2045)
  newstressSuit$year2055 <- ifelse(newstressSuit$year2050 == 1, ifelse(newstressSuit$X2055==0,2055,1),newstressSuit$year2050)
  newstressSuit$year2060 <- ifelse(newstressSuit$year2055 == 1, ifelse(newstressSuit$X2060==0,2060,1),newstressSuit$year2055)
  newstressSuit$year2065 <- ifelse(newstressSuit$year2060 == 1, ifelse(newstressSuit$X2065==0,2065,1),newstressSuit$year2060)
  newstressSuit$year2070 <- ifelse(newstressSuit$year2065 == 1, ifelse(newstressSuit$X2070==0,2070,1),newstressSuit$year2065)
  newstressSuit$year2075 <- ifelse(newstressSuit$year2070 == 1, ifelse(newstressSuit$X2075==0,2075,1),newstressSuit$year2070)
  newstressSuit$year2080 <- ifelse(newstressSuit$year2075 == 1, ifelse(newstressSuit$X2080==0,2080,1),newstressSuit$year2075)
  newstressSuit$year2085 <- ifelse(newstressSuit$year2080 == 1, ifelse(newstressSuit$X2085==0,2085,1),newstressSuit$year2080)
  newstressSuit$year2090 <- ifelse(newstressSuit$year2085 == 1, ifelse(newstressSuit$X2090==0,2090,1),newstressSuit$year2085)
  newstressSuit$year2095 <- ifelse(newstressSuit$year2090 == 1, ifelse(newstressSuit$X2095==0,2095,1),newstressSuit$year2090)
  newstressSuit$year2100 <- ifelse(newstressSuit$year2095 == 1, ifelse(newstressSuit$X2100==0,2100,1),newstressSuit$year2095)
  newstressSuit$expyr <- ifelse(newstressSuit$year2100 == 1, 2105,newstressSuit$year2100)
  
  newstressSuit$expyr[newstressSuit$expyr == 0]<- NA
  
  
  
  #compare exp year with the rest of the variables: earliest year = overall exp year
  #import csv for each stressor
  str1files <- list.files(stressordir1,pattern="final_yrexp",full.names=T)
  str2files <- list.files(stressordir2,pattern="final_yrexp",full.names=T)
  str3files <- list.files(stressordir3,pattern="final_yrexp",full.names=T)
  str4files <- list.files(stressordir4,pattern="final_yrexp",full.names=T)
  
  strexp <- fread(file=str1files[1],header=T, select=c("ID","expyr"))
  str2exp <- fread(file=str2files[1], header=T, select=c("ID","expyr"))
  strexp <-merge(strexp,str2exp,by="ID")
  str3exp <- fread(file=str3files[1],header=T, select=c("ID","expyr"))
  strexp <-merge(strexp,str3exp,by="ID")
  str4exp <- fread(file=str4files[1],header=T, select=c("ID","expyr"))
  strexp <-merge(strexp,str4exp,by="ID")
  names(strexp) <- c("ID","strexp1","strexp2","strexp3","strexp4")
  newstressexp <- newstressSuit[, .(ID, expyr)]
  strexp <- merge(strexp,newstressexp,by="ID")
  
  #calculate overall expiration (minimum expiration year = overall expiration for that site)
  strexp <- strexp[, .(strexp1,strexp2,strexp3,strexp4,expyr)]
  strexp$overallexp <- apply(strexp, 1, FUN = min, na.rm = TRUE)
  
  #calculate median expiration year
  medexp <- median(strexp$overallexp)
  print(medexp)
  
  finish <- Sys.time()
  finish-start
}


###adjust threshold values as needed

#dhw #7.2 7.6 8 8.4 8.8
newthreshold(threshold = 8.8, stressor = "dhw",newthreshstressordir= "dhw/rcp45",stressordir1="land_hurtt/rcp45",
             stressordir2="oa/rcp45",stressordir3="storms/rcp45",stressordir4="pop_hr/ssp2")
#oa #3.63 3.465 3.3 3.135 2.97
newthreshold(threshold = 2.97, stressor = "oa",newthreshstressordir= "oa/rcp45",stressordir1="dhw/rcp45",
             stressordir2="land_hurtt/rcp45",stressordir3="storms/rcp45",stressordir4="pop_hr/ssp2")
#storms 4.5 4.75 5 5.25 5.5
newthreshold(threshold = 5.5, stressor = "storms",newthreshstressordir= "storm/rcp45",stressordir1="dhw/rcp45",
             stressordir2="oa/rcp45",stressordir3="land_hurtt/rcp45",stressordir4="pop_hr/ssp2")
#pop 1.8 1.9 2 2.1 2.2
newthreshold(threshold = 2.2, stressor = "pop",newthreshstressordir= "pop/ssp2",stressordir1="dhw/rcp45",
             stressordir2="oa/rcp45",stressordir3="storms/rcp45",stressordir4="land_hurtt/rcp45")
#land 0.45 0.475 0.5 0.525 0.55
newthreshold(threshold = 0.55, stressor = "land",newthreshstressordir= "land/ssp2",stressordir1="dhw/rcp45",
             stressordir2="oa/rcp45",stressordir3="storms/rcp45",stressordir4="pop_hr/ssp2")


#create plot: median expiration v threshold value with all variables
library(ggplot2)
library(tidyverse)

expthdfmed <- data.frame(threshold_change = c(-0.1,-0.05,0,0.05,0.1),
                         dhw=c(2020,2020,2020,2020,2025),
                         oa=c(2010,2020,2020,2025,2025), 
                         storms=c(2020,2020,2020,2020,2020),
                         pop=c(2005,2015,2020,2020,2025),
                         land=c(2020,2020,2020,2020,2020))

expthmed <- expthdfmed %>%
  select(threshold_change,dhw,oa,storms,pop,land) %>%
  gather(key = "Variable", value = "ExpirationYear", -threshold_change)

thmed_plot <- ggplot(expthmed, aes(x=threshold_change,y=ExpirationYear,group=Variable)) + 
  geom_line(aes(color = Variable,alpha=Variable),size=1,position = position_jitter(width = 0, h=0.2)) + #,position = position_dodge(width = 0.01)
  scale_alpha_manual(values=c(0.7,0.7,0.7,0.7,0.7)) 

thmed_plot + theme_bw() + theme(
  plot.title = element_text(face = "bold", size = 12),
  legend.background = element_rect(fill = "white", size = 4, colour = "white"),
  axis.ticks = element_line(colour = "grey70", size = 0.2),
  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
  strip.text.x = element_text(size=10, angle=0),
  strip.text.y = element_text(size=12, face="bold"),
  panel.grid.major = element_line(colour = "grey70", size = 0.2),
  panel.grid.minor = element_blank()
)+
  xlab("Threshold percent change")+ylab("Year")+ #ggtitle("Median Remaining Expiration Year")+
  scale_color_manual(name="Stressor",labels=c("DHW","Land use","OA","Population","Storms"),
                     values=c("#66C2A5", "#FC8D62", "#E78AC3", "#A6D854", "#FFD92F", "#4E555A"))+
  guides(alpha=FALSE)+
  scale_x_continuous(labels = scales::percent_format(1)) 

