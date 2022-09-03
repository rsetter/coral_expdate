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
library(fasterize)
library(viridisLite)
library(gridExtra)
library(rasterVis)
library(sp)
library(maps)
library(maptools)
library(data.table)



reefs <- st_read("reefs.shp")


#function to calc year that sites become unsuit 
expyr = function(wdhist,wdrcp) {
  setwd(wdrcp)
  csvfiles <- list.files(wdhist,pattern="suit.csv$",full.names=T)
  ncsvfiles <- length(csvfiles)
  acsv <- read.csv(file=csvfiles[1],row.names=NULL,header=T)
  keep <- acsv[c("ID","iscoral","suit")]
  keep1 <- keep[!(keep$iscoral==0),]
  keep2 <- keep1[c("ID","suit")]
  colnames(keep2) <- c("ID","X1855")
  suitcoralall <- distinct(keep2)
  
  start <- Sys.time()
  cores<- detectCores()-1
  cl <- makeCluster(cores, output="")
  registerDoParallel(cl)
  for (i in 2:ncsvfiles){
    a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
    keep <- a[c("ID","iscoral","suit")]
    keep1 <- keep[!(keep$iscoral==0),]
    keep2 <- keep1[c("ID","suit")]
    colnames(keep2) <- c("ID",paste0("X",i*10+1845))
    keep2 <- distinct(keep2)
    suitcoralall <- left_join(suitcoralall,keep2,by="ID")
  } 
  stopCluster(cl)
  finish <- Sys.time()
  finish-start
  write.csv(suitcoralall,"hist_allyr.csv", row.names=FALSE)
  
  #get csv files. keep suitcoral column. rename column the correct year. merge together into one dataframe
  csvfiles <- list.files(wdrcp,pattern="suit.csv$",full.names=T)
  ncsvfiles <- length(csvfiles)
  
  start <- Sys.time()
  cores<- detectCores()-1
  cl <- makeCluster(cores, output="")
  registerDoParallel(cl)
  for (i in 1:ncsvfiles){
    a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
    keep <- a[c("ID","iscoral","suit")]
    keep1 <- keep[!(keep$iscoral==0),]
    keep2 <- keep1[c("ID","suit")]
    colnames(keep2) <- c("ID",paste0("X",i*5+2000))
    keep2 <- distinct(keep2)
    suitcoralall <- left_join(suitcoralall,keep2,by="ID")
  }
  stopCluster(cl)
  finish <- Sys.time()
  finish-start
  write.csv(suitcoralall,"rcp85_allyr.csv", row.names=FALSE)
  
  #find year site becomes unsuitable
  
  suitcoralall$year1855 <- ifelse(suitcoralall$X1855 == 0, 1855,1)
  suitcoralall$year1865 <- ifelse(suitcoralall$year1855 == 1, ifelse(suitcoralall$X1865 ==0,1865,1),suitcoralall$year1855)
  suitcoralall$year1875 <- ifelse(suitcoralall$year1865 == 1, ifelse(suitcoralall$X1875 ==0,1875,1),suitcoralall$year1865)
  suitcoralall$year1885 <- ifelse(suitcoralall$year1875 == 1, ifelse(suitcoralall$X1885 ==0,1885,1),suitcoralall$year1875)
  suitcoralall$year1895 <- ifelse(suitcoralall$year1855 == 1, ifelse(suitcoralall$X1895 ==0,1895,1),suitcoralall$year1885)
  suitcoralall$year1905 <- ifelse(suitcoralall$year1895 == 1, ifelse(suitcoralall$X1905 ==0,1905,1),suitcoralall$year1895)
  suitcoralall$year1915 <- ifelse(suitcoralall$year1905 == 1, ifelse(suitcoralall$X1915 ==0,1915,1),suitcoralall$year1905)
  suitcoralall$year1925 <- ifelse(suitcoralall$year1915 == 1, ifelse(suitcoralall$X1925 ==0,1925,1),suitcoralall$year1915)
  suitcoralall$year1935 <- ifelse(suitcoralall$year1925 == 1, ifelse(suitcoralall$X1935 ==0,1935,1),suitcoralall$year1925)
  suitcoralall$year1945 <- ifelse(suitcoralall$year1935 == 1, ifelse(suitcoralall$X1945 ==0,1945,1),suitcoralall$year1935)
  suitcoralall$year1955 <- ifelse(suitcoralall$year1945 == 1, ifelse(suitcoralall$X1955 ==0,1955,1),suitcoralall$year1945)
  suitcoralall$year1965 <- ifelse(suitcoralall$year1955 == 1, ifelse(suitcoralall$X1965 ==0,1965,1),suitcoralall$year1955)
  suitcoralall$year1975 <- ifelse(suitcoralall$year1965 == 1, ifelse(suitcoralall$X1975 ==0,1975,1),suitcoralall$year1965)
  suitcoralall$year1985 <- ifelse(suitcoralall$year1975 == 1, ifelse(suitcoralall$X1985 ==0,1985,1),suitcoralall$year1975)
  suitcoralall$year1995 <- ifelse(suitcoralall$year1985 == 1, ifelse(suitcoralall$X1995 ==0,1995,1),suitcoralall$year1985)
  suitcoralall$year2005 <- ifelse(suitcoralall$year1995 == 1, ifelse(suitcoralall$X2005.y ==0,2005,1),suitcoralall$year1995)
  suitcoralall$year2010 <- ifelse(suitcoralall$year2005 == 1, ifelse(suitcoralall$X2010 ==0,2010,1),suitcoralall$year2005)
  suitcoralall$year2015 <- ifelse(suitcoralall$year2010 == 1, ifelse(suitcoralall$X2015 ==0,2015,1),suitcoralall$year2010)
  suitcoralall$year2020 <- ifelse(suitcoralall$year2015 == 1, ifelse(suitcoralall$X2020 ==0,2020,1),suitcoralall$year2015)
  suitcoralall$year2025 <- ifelse(suitcoralall$year2020 == 1, ifelse(suitcoralall$X2025==0,2025,1),suitcoralall$year2020)
  suitcoralall$year2030 <- ifelse(suitcoralall$year2025 == 1, ifelse(suitcoralall$X2030==0,2030,1),suitcoralall$year2025)
  suitcoralall$year2035 <- ifelse(suitcoralall$year2030 == 1, ifelse(suitcoralall$X2035==0,2035,1),suitcoralall$year2030)
  suitcoralall$year2040 <- ifelse(suitcoralall$year2035 == 1, ifelse(suitcoralall$X2040==0,2040,1),suitcoralall$year2035)
  suitcoralall$year2045 <- ifelse(suitcoralall$year2040 == 1, ifelse(suitcoralall$X2045==0,2045,1),suitcoralall$year2040)
  suitcoralall$year2050 <- ifelse(suitcoralall$year2045 == 1, ifelse(suitcoralall$X2050==0,2050,1),suitcoralall$year2045)
  suitcoralall$year2055 <- ifelse(suitcoralall$year2050 == 1, ifelse(suitcoralall$X2055==0,2055,1),suitcoralall$year2050)
  suitcoralall$year2060 <- ifelse(suitcoralall$year2055 == 1, ifelse(suitcoralall$X2060==0,2060,1),suitcoralall$year2055)
  suitcoralall$year2065 <- ifelse(suitcoralall$year2060 == 1, ifelse(suitcoralall$X2065==0,2065,1),suitcoralall$year2060)
  suitcoralall$year2070 <- ifelse(suitcoralall$year2065 == 1, ifelse(suitcoralall$X2070==0,2070,1),suitcoralall$year2065)
  suitcoralall$year2075 <- ifelse(suitcoralall$year2070 == 1, ifelse(suitcoralall$X2075==0,2075,1),suitcoralall$year2070)
  suitcoralall$year2080 <- ifelse(suitcoralall$year2075 == 1, ifelse(suitcoralall$X2080==0,2080,1),suitcoralall$year2075)
  suitcoralall$year2085 <- ifelse(suitcoralall$year2080 == 1, ifelse(suitcoralall$X2085==0,2085,1),suitcoralall$year2080)
  suitcoralall$year2090 <- ifelse(suitcoralall$year2085 == 1, ifelse(suitcoralall$X2090==0,2090,1),suitcoralall$year2085)
  suitcoralall$year2095 <- ifelse(suitcoralall$year2090 == 1, ifelse(suitcoralall$X2095==0,2095,1),suitcoralall$year2090)
  suitcoralall$year2100 <- ifelse(suitcoralall$year2095 == 1, ifelse(suitcoralall$X2100==0,2100,1),suitcoralall$year2095)
  suitcoralall$expyr <- ifelse(suitcoralall$year2100 == 1, 2105,suitcoralall$year2100)
  
  
  
  suitcoralall$expyr[suitcoralall$expyr == 0]<- NA
  
  write.csv(suitcoralall,"rcp85_yrexp.csv", row.names=FALSE)
  
}


#specifically for storms
expyrstorm = function(wdhist,wdrcp) {
  setwd(wdrcp)
  csvfiles <- list.files(wdhist,pattern="suit.csv$",full.names=T)
  ncsvfiles <- length(csvfiles)
  acsv <- read.csv(file=csvfiles[1],row.names=NULL,header=T)
  keep <- acsv[c("ID","iscoral","suit")]
  keep1 <- keep[!(keep$iscoral==0),]
  keep2 <- keep1[c("ID","suit")]
  colnames(keep2) <- c("ID","X1950")
  suitcoralall <- distinct(keep2)
  
  start <- Sys.time()
  cores<- detectCores()-1
  cl <- makeCluster(cores, output="")
  registerDoParallel(cl)
  for (i in 2:ncsvfiles){
    a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
    keep <- a[c("ID","iscoral","suit")]
    keep1 <- keep[!(keep$iscoral==0),]
    keep2 <- keep1[c("ID","suit")]
    colnames(keep2) <- c("ID",paste0("X",i*5+1945))
    keep2 <- distinct(keep2)
    suitcoralall <- left_join(suitcoralall,keep2,by="ID")
  } 
  stopCluster(cl)
  finish <- Sys.time()
  finish-start
  write.csv(suitcoralall,"hist_allyr.csv", row.names=FALSE)
  
  #get csv files. keep suitcoral column. rename column the correct year. merge together into one dataframe
  csvfiles <- list.files(wdrcp,pattern="suit.csv$",full.names=T)
  ncsvfiles <- length(csvfiles)
  
  start <- Sys.time()
  cores<- detectCores()-1
  cl <- makeCluster(cores, output="")
  registerDoParallel(cl)
  for (i in 1:ncsvfiles){
    a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
    keep <- a[c("ID","iscoral","suit")]
    keep1 <- keep[!(keep$iscoral==0),]
    keep2 <- keep1[c("ID","suit")]
    colnames(keep2) <- c("ID",paste0("X",i*5+2000))
    keep2 <- distinct(keep2)
    suitcoralall <- left_join(suitcoralall,keep2,by="ID")
  }
  stopCluster(cl)
  finish <- Sys.time()
  finish-start
  write.csv(suitcoralall,"rcp85_allyr.csv", row.names=FALSE)
  
  #find year site becomes unsuitable
  
  suitcoralall$year1950 <- ifelse(suitcoralall$X1950 == 0, 1950,1)
  suitcoralall$year1955 <- ifelse(suitcoralall$year1950 == 1, ifelse(suitcoralall$X1955 ==0,1955,1),suitcoralall$year1950)
  suitcoralall$year1960 <- ifelse(suitcoralall$year1955 == 1, ifelse(suitcoralall$X1960 ==0,1960,1),suitcoralall$year1955)
  suitcoralall$year1965 <- ifelse(suitcoralall$year1960 == 1, ifelse(suitcoralall$X1965 ==0,1965,1),suitcoralall$year1960)
  suitcoralall$year1970 <- ifelse(suitcoralall$year1965 == 1, ifelse(suitcoralall$X1970 ==0,1970,1),suitcoralall$year1965)
  suitcoralall$year1975 <- ifelse(suitcoralall$year1970 == 1, ifelse(suitcoralall$X1975 ==0,1975,1),suitcoralall$year1970)
  suitcoralall$year1980 <- ifelse(suitcoralall$year1975 == 1, ifelse(suitcoralall$X1980 ==0,1980,1),suitcoralall$year1975)
  suitcoralall$year1985 <- ifelse(suitcoralall$year1980 == 1, ifelse(suitcoralall$X1985 ==0,1985,1),suitcoralall$year1980)
  suitcoralall$year1990 <- ifelse(suitcoralall$year1985 == 1, ifelse(suitcoralall$X1990 ==0,1990,1),suitcoralall$year1985)
  suitcoralall$year1995 <- ifelse(suitcoralall$year1990 == 1, ifelse(suitcoralall$X1995 ==0,1995,1),suitcoralall$year1990)
  suitcoralall$year2000 <- ifelse(suitcoralall$year1995 == 1, ifelse(suitcoralall$X2000 ==0,2000,1),suitcoralall$year1995)
  suitcoralall$year2005 <- ifelse(suitcoralall$year2000 == 1, ifelse(suitcoralall$X2005.y ==0,2005,1),suitcoralall$year2000)
  suitcoralall$year2010 <- ifelse(suitcoralall$year2005 == 1, ifelse(suitcoralall$X2010 ==0,2010,1),suitcoralall$year2005)
  suitcoralall$year2015 <- ifelse(suitcoralall$year2010 == 1, ifelse(suitcoralall$X2015 ==0,2015,1),suitcoralall$year2010)
  suitcoralall$year2020 <- ifelse(suitcoralall$year2015 == 1, ifelse(suitcoralall$X2020 ==0,2020,1),suitcoralall$year2015)
  suitcoralall$year2025 <- ifelse(suitcoralall$year2020 == 1, ifelse(suitcoralall$X2025==0,2025,1),suitcoralall$year2020)
  suitcoralall$year2030 <- ifelse(suitcoralall$year2025 == 1, ifelse(suitcoralall$X2030==0,2030,1),suitcoralall$year2025)
  suitcoralall$year2035 <- ifelse(suitcoralall$year2030 == 1, ifelse(suitcoralall$X2035==0,2035,1),suitcoralall$year2030)
  suitcoralall$year2040 <- ifelse(suitcoralall$year2035 == 1, ifelse(suitcoralall$X2040==0,2040,1),suitcoralall$year2035)
  suitcoralall$year2045 <- ifelse(suitcoralall$year2040 == 1, ifelse(suitcoralall$X2045==0,2045,1),suitcoralall$year2040)
  suitcoralall$year2050 <- ifelse(suitcoralall$year2045 == 1, ifelse(suitcoralall$X2050==0,2050,1),suitcoralall$year2045)
  suitcoralall$year2055 <- ifelse(suitcoralall$year2050 == 1, ifelse(suitcoralall$X2055==0,2055,1),suitcoralall$year2050)
  suitcoralall$year2060 <- ifelse(suitcoralall$year2055 == 1, ifelse(suitcoralall$X2060==0,2060,1),suitcoralall$year2055)
  suitcoralall$year2065 <- ifelse(suitcoralall$year2060 == 1, ifelse(suitcoralall$X2065==0,2065,1),suitcoralall$year2060)
  suitcoralall$year2070 <- ifelse(suitcoralall$year2065 == 1, ifelse(suitcoralall$X2070==0,2070,1),suitcoralall$year2065)
  suitcoralall$year2075 <- ifelse(suitcoralall$year2070 == 1, ifelse(suitcoralall$X2075==0,2075,1),suitcoralall$year2070)
  suitcoralall$year2080 <- ifelse(suitcoralall$year2075 == 1, ifelse(suitcoralall$X2080==0,2080,1),suitcoralall$year2075)
  suitcoralall$year2085 <- ifelse(suitcoralall$year2080 == 1, ifelse(suitcoralall$X2085==0,2085,1),suitcoralall$year2080)
  suitcoralall$year2090 <- ifelse(suitcoralall$year2085 == 1, ifelse(suitcoralall$X2090==0,2090,1),suitcoralall$year2085)
  suitcoralall$year2095 <- ifelse(suitcoralall$year2090 == 1, ifelse(suitcoralall$X2095==0,2095,1),suitcoralall$year2090)
  suitcoralall$year2100 <- ifelse(suitcoralall$year2095 == 1, ifelse(suitcoralall$X2100==0,2100,1),suitcoralall$year2095)
  suitcoralall$expyr <- ifelse(suitcoralall$year2100 == 1, 2105,suitcoralall$year2100)
  
  
  
  suitcoralall$expyr[suitcoralall$expyr == 0]<- NA
  
  write.csv(suitcoralall,"rcp85_yrexp.csv", row.names=FALSE)
  
}


#function to make unsuitability shapefile 
expirationyr = function(wdrcp) {
  
  setwd(wdrcp)
  csvfiles <- list.files(wdrcp,pattern="yrexp",full.names=T)
  suitcoralall <- fread(csvfiles[1],header=T,select=c(1,73))
  
  
  #merge with shapefile and save
  suitcoralyr <- merge(reefs,suitcoralall,by="ID")
  st_write(suitcoralyr, dsn = "exp_yr", driver = "ESRI Shapefile")
  
}

#function to recalc expiration, such that unsuitable reefs are only those that are unsuitable forever
expyrnew <- function(wdrcp) {
  setwd(wdrcp)
  csvfiles <- list.files(wdrcp,pattern="allyr.csv$",full.names=T)
  
  #open csv. this one has all years from 1855-2100, only coral points
  suitcoralall <-fread(file=csvfiles[2],header=T)

  #modify unsuitability date: remove patterns where reef goes unsuitable to not to unsuitable
  start <- Sys.time()
  cores<- detectCores()-1
  cl <- makeCluster(cores, output="")
  registerDoParallel(cl)
  cols <- colnames(suitcoralall)
    for (j in (ncol(suitcoralall)-1):2){
      column <- cols[j]
      future <- cols[j+1]
      suitcoralall[[column]][suitcoralall[[future]]>suitcoralall[[column]]] <- 1
    } 
  
  stopCluster(cl)
  finish <- Sys.time()
  finish-start #total time=
  write.csv(suitcoralall,"final_allyr.csv", row.names=FALSE)
  
  #find year site becomes unsuitable
  
  suitcoralall$year1855 <- ifelse(suitcoralall$X1855 == 0, 1855,1)
  suitcoralall$year1865 <- ifelse(suitcoralall$year1855 == 1, ifelse(suitcoralall$X1865 ==0,1865,1),suitcoralall$year1855)
  suitcoralall$year1875 <- ifelse(suitcoralall$year1865 == 1, ifelse(suitcoralall$X1875 ==0,1875,1),suitcoralall$year1865)
  suitcoralall$year1885 <- ifelse(suitcoralall$year1875 == 1, ifelse(suitcoralall$X1885 ==0,1885,1),suitcoralall$year1875)
  suitcoralall$year1895 <- ifelse(suitcoralall$year1855 == 1, ifelse(suitcoralall$X1895 ==0,1895,1),suitcoralall$year1885)
  suitcoralall$year1905 <- ifelse(suitcoralall$year1895 == 1, ifelse(suitcoralall$X1905 ==0,1905,1),suitcoralall$year1895)
  suitcoralall$year1915 <- ifelse(suitcoralall$year1905 == 1, ifelse(suitcoralall$X1915 ==0,1915,1),suitcoralall$year1905)
  suitcoralall$year1925 <- ifelse(suitcoralall$year1915 == 1, ifelse(suitcoralall$X1925 ==0,1925,1),suitcoralall$year1915)
  suitcoralall$year1935 <- ifelse(suitcoralall$year1925 == 1, ifelse(suitcoralall$X1935 ==0,1935,1),suitcoralall$year1925)
  suitcoralall$year1945 <- ifelse(suitcoralall$year1935 == 1, ifelse(suitcoralall$X1945 ==0,1945,1),suitcoralall$year1935)
  suitcoralall$year1955 <- ifelse(suitcoralall$year1945 == 1, ifelse(suitcoralall$X1955 ==0,1955,1),suitcoralall$year1945)
  suitcoralall$year1965 <- ifelse(suitcoralall$year1955 == 1, ifelse(suitcoralall$X1965 ==0,1965,1),suitcoralall$year1955)
  suitcoralall$year1975 <- ifelse(suitcoralall$year1965 == 1, ifelse(suitcoralall$X1975 ==0,1975,1),suitcoralall$year1965)
  suitcoralall$year1985 <- ifelse(suitcoralall$year1975 == 1, ifelse(suitcoralall$X1985 ==0,1985,1),suitcoralall$year1975)
  suitcoralall$year1995 <- ifelse(suitcoralall$year1985 == 1, ifelse(suitcoralall$X1995 ==0,1995,1),suitcoralall$year1985)
  suitcoralall$year2005h <- ifelse(suitcoralall$year1995 == 1, ifelse(suitcoralall$X2005.x ==0,2005,1),suitcoralall$year1995)
  suitcoralall$year2005 <- ifelse(suitcoralall$year1995 == 1, ifelse(suitcoralall$X2005.y ==0,2005,1),suitcoralall$year1995)
  suitcoralall$year2010 <- ifelse(suitcoralall$year2005 == 1, ifelse(suitcoralall$X2010 ==0,2010,1),suitcoralall$year2005)
  suitcoralall$year2015 <- ifelse(suitcoralall$year2010 == 1, ifelse(suitcoralall$X2015 ==0,2015,1),suitcoralall$year2010)
  suitcoralall$year2020 <- ifelse(suitcoralall$year2015 == 1, ifelse(suitcoralall$X2020 ==0,2020,1),suitcoralall$year2015)
  suitcoralall$year2025 <- ifelse(suitcoralall$year2020 == 1, ifelse(suitcoralall$X2025==0,2025,1),suitcoralall$year2020)
  suitcoralall$year2030 <- ifelse(suitcoralall$year2025 == 1, ifelse(suitcoralall$X2030==0,2030,1),suitcoralall$year2025)
  suitcoralall$year2035 <- ifelse(suitcoralall$year2030 == 1, ifelse(suitcoralall$X2035==0,2035,1),suitcoralall$year2030)
  suitcoralall$year2040 <- ifelse(suitcoralall$year2035 == 1, ifelse(suitcoralall$X2040==0,2040,1),suitcoralall$year2035)
  suitcoralall$year2045 <- ifelse(suitcoralall$year2040 == 1, ifelse(suitcoralall$X2045==0,2045,1),suitcoralall$year2040)
  suitcoralall$year2050 <- ifelse(suitcoralall$year2045 == 1, ifelse(suitcoralall$X2050==0,2050,1),suitcoralall$year2045)
  suitcoralall$year2055 <- ifelse(suitcoralall$year2050 == 1, ifelse(suitcoralall$X2055==0,2055,1),suitcoralall$year2050)
  suitcoralall$year2060 <- ifelse(suitcoralall$year2055 == 1, ifelse(suitcoralall$X2060==0,2060,1),suitcoralall$year2055)
  suitcoralall$year2065 <- ifelse(suitcoralall$year2060 == 1, ifelse(suitcoralall$X2065==0,2065,1),suitcoralall$year2060)
  suitcoralall$year2070 <- ifelse(suitcoralall$year2065 == 1, ifelse(suitcoralall$X2070==0,2070,1),suitcoralall$year2065)
  suitcoralall$year2075 <- ifelse(suitcoralall$year2070 == 1, ifelse(suitcoralall$X2075==0,2075,1),suitcoralall$year2070)
  suitcoralall$year2080 <- ifelse(suitcoralall$year2075 == 1, ifelse(suitcoralall$X2080==0,2080,1),suitcoralall$year2075)
  suitcoralall$year2085 <- ifelse(suitcoralall$year2080 == 1, ifelse(suitcoralall$X2085==0,2085,1),suitcoralall$year2080)
  suitcoralall$year2090 <- ifelse(suitcoralall$year2085 == 1, ifelse(suitcoralall$X2090==0,2090,1),suitcoralall$year2085)
  suitcoralall$year2095 <- ifelse(suitcoralall$year2090 == 1, ifelse(suitcoralall$X2095==0,2095,1),suitcoralall$year2090)
  suitcoralall$year2100 <- ifelse(suitcoralall$year2095 == 1, ifelse(suitcoralall$X2100==0,2100,1),suitcoralall$year2095)
  suitcoralall$expyr <- ifelse(suitcoralall$year2100 == 1, 2105,suitcoralall$year2100)

  
  suitcoralall$expyr[suitcoralall$expyr == 0]<- NA
  
  write.csv(suitcoralall,"final_yrexp.csv", row.names=FALSE)
} #timer: ~45min


#function to recalc overall expiration date using new stressor expiration dates
expyroverallnew  <- function(wdrcp,wddhwrcp, wdlandrcp, wdoarcp,wdpoprcp,wdstormrcp) {
  setwd(wdrcp)
  #load dhw
  dhwfiles <- list.files(wddhwrcp,pattern="final",full.names=T)
  #load land
  landfiles <- list.files(wdlandrcp,pattern="final",full.names=T)
  #load oa
  oafiles <- list.files(wdoarcp,pattern="final",full.names=T)
  #load pop
  popfiles <- list.files(wdpoprcp,pattern="final",full.names=T)
  #load storms
  stormfiles <- list.files(wdstormrcp,pattern="final",full.names=T)
  
  
  
  #multiply columns together. save csv
  start <- Sys.time()
  cores<- detectCores()-3
  cl <- makeCluster(cores, output="") 
  registerDoParallel(cl)
  cols <- c( "X1855","X1865","X1875","X1885","X1895","X1905","X1915","X1925","X1935","X1945","X1955","X1965","X1975","X1985","X1995",
             "X2005.x","X2005.y","X2010","X2015","X2020","X2025","X2030","X2035","X2040","X2045","X2050","X2055","X2060","X2065","X2070","X2075",
             "X2080","X2085","X2090","X2095","X2100")
  overall <- fread(dhwfiles[1], header=T, select=c(1))
  
  for (i in 2:37) {
    column <- cols[i-1]
    
    #load each file
    dhw <- fread(file=dhwfiles[1],header=T, select=c(1,i))
    dhw[is.na(dhw)] <- 1
    dhw <- dhw[order(ID)]
    land <- fread(file=landfiles[1], header=T, select=c(1,i))
    land[is.na(land)] <- 1
    land <- land[order(ID)]
    oa <- fread(file=oafiles[1],header=T, select=c(1,i))
    oa[is.na(oa)] <- 1
    oa <- oa[order(ID)]
    pop <- fread(file=popfiles[1],header=T, select=c(1,i))
    pop[is.na(pop)] <- 1
    pop <- pop[order(ID)]
    storms <- fread(file=stormfiles[1],header=T, select=c(1,i))
    storms[is.na(storms)] <- 1
    storms <- storms[order(ID)]
    
    #join together overall suit
    overall[[column]] <- dhw[[column]] * land[[column]] * oa[[column]] * pop[[column]] * storms[[column]]
    
  }
  
  write.csv(overall,"final_allyr.csv", row.names=FALSE)
  
  stopCluster(cl)
  finish <- Sys.time()
  finish-start
  
  #keep a copy for expyr calculation below
  suitcoralall <- overall
  
  #find the count for number of stressors in 2100
  overall[["nstressors2100"]] <- dhw[["X2100"]] + land[["X2100"]] + oa[["X2100"]] + pop[["X2100"]] + storms[["X2100"]]
  stress <- subset(overall, select=c(ID,X2100,nstressors2100))
  write.csv(stress,"final_stressors2100.csv", row.names=FALSE)
  
  
  #find year site becomes unsuitable
  suitcoralall$year1855 <- ifelse(suitcoralall$X1855 == 0, 1855,1)
  suitcoralall$year1865 <- ifelse(suitcoralall$year1855 == 1, ifelse(suitcoralall$X1865 ==0,1865,1),suitcoralall$year1855)
  suitcoralall$year1875 <- ifelse(suitcoralall$year1865 == 1, ifelse(suitcoralall$X1875 ==0,1875,1),suitcoralall$year1865)
  suitcoralall$year1885 <- ifelse(suitcoralall$year1875 == 1, ifelse(suitcoralall$X1885 ==0,1885,1),suitcoralall$year1875)
  suitcoralall$year1895 <- ifelse(suitcoralall$year1855 == 1, ifelse(suitcoralall$X1895 ==0,1895,1),suitcoralall$year1885)
  suitcoralall$year1905 <- ifelse(suitcoralall$year1895 == 1, ifelse(suitcoralall$X1905 ==0,1905,1),suitcoralall$year1895)
  suitcoralall$year1915 <- ifelse(suitcoralall$year1905 == 1, ifelse(suitcoralall$X1915 ==0,1915,1),suitcoralall$year1905)
  suitcoralall$year1925 <- ifelse(suitcoralall$year1915 == 1, ifelse(suitcoralall$X1925 ==0,1925,1),suitcoralall$year1915)
  suitcoralall$year1935 <- ifelse(suitcoralall$year1925 == 1, ifelse(suitcoralall$X1935 ==0,1935,1),suitcoralall$year1925)
  suitcoralall$year1945 <- ifelse(suitcoralall$year1935 == 1, ifelse(suitcoralall$X1945 ==0,1945,1),suitcoralall$year1935)
  suitcoralall$year1955 <- ifelse(suitcoralall$year1945 == 1, ifelse(suitcoralall$X1955 ==0,1955,1),suitcoralall$year1945)
  suitcoralall$year1965 <- ifelse(suitcoralall$year1955 == 1, ifelse(suitcoralall$X1965 ==0,1965,1),suitcoralall$year1955)
  suitcoralall$year1975 <- ifelse(suitcoralall$year1965 == 1, ifelse(suitcoralall$X1975 ==0,1975,1),suitcoralall$year1965)
  suitcoralall$year1985 <- ifelse(suitcoralall$year1975 == 1, ifelse(suitcoralall$X1985 ==0,1985,1),suitcoralall$year1975)
  suitcoralall$year1995 <- ifelse(suitcoralall$year1985 == 1, ifelse(suitcoralall$X1995 ==0,1995,1),suitcoralall$year1985)
  suitcoralall$year2005 <- ifelse(suitcoralall$year1995 == 1, ifelse(suitcoralall$X2005.y ==0,2005,1),suitcoralall$year1995)
  suitcoralall$year2010 <- ifelse(suitcoralall$year2005 == 1, ifelse(suitcoralall$X2010 ==0,2010,1),suitcoralall$year2005)
  suitcoralall$year2015 <- ifelse(suitcoralall$year2010 == 1, ifelse(suitcoralall$X2015 ==0,2015,1),suitcoralall$year2010)
  suitcoralall$year2020 <- ifelse(suitcoralall$year2015 == 1, ifelse(suitcoralall$X2020 ==0,2020,1),suitcoralall$year2015)
  suitcoralall$year2025 <- ifelse(suitcoralall$year2020 == 1, ifelse(suitcoralall$X2025==0,2025,1),suitcoralall$year2020)
  suitcoralall$year2030 <- ifelse(suitcoralall$year2025 == 1, ifelse(suitcoralall$X2030==0,2030,1),suitcoralall$year2025)
  suitcoralall$year2035 <- ifelse(suitcoralall$year2030 == 1, ifelse(suitcoralall$X2035==0,2035,1),suitcoralall$year2030)
  suitcoralall$year2040 <- ifelse(suitcoralall$year2035 == 1, ifelse(suitcoralall$X2040==0,2040,1),suitcoralall$year2035)
  suitcoralall$year2045 <- ifelse(suitcoralall$year2040 == 1, ifelse(suitcoralall$X2045==0,2045,1),suitcoralall$year2040)
  suitcoralall$year2050 <- ifelse(suitcoralall$year2045 == 1, ifelse(suitcoralall$X2050==0,2050,1),suitcoralall$year2045)
  suitcoralall$year2055 <- ifelse(suitcoralall$year2050 == 1, ifelse(suitcoralall$X2055==0,2055,1),suitcoralall$year2050)
  suitcoralall$year2060 <- ifelse(suitcoralall$year2055 == 1, ifelse(suitcoralall$X2060==0,2060,1),suitcoralall$year2055)
  suitcoralall$year2065 <- ifelse(suitcoralall$year2060 == 1, ifelse(suitcoralall$X2065==0,2065,1),suitcoralall$year2060)
  suitcoralall$year2070 <- ifelse(suitcoralall$year2065 == 1, ifelse(suitcoralall$X2070==0,2070,1),suitcoralall$year2065)
  suitcoralall$year2075 <- ifelse(suitcoralall$year2070 == 1, ifelse(suitcoralall$X2075==0,2075,1),suitcoralall$year2070)
  suitcoralall$year2080 <- ifelse(suitcoralall$year2075 == 1, ifelse(suitcoralall$X2080==0,2080,1),suitcoralall$year2075)
  suitcoralall$year2085 <- ifelse(suitcoralall$year2080 == 1, ifelse(suitcoralall$X2085==0,2085,1),suitcoralall$year2080)
  suitcoralall$year2090 <- ifelse(suitcoralall$year2085 == 1, ifelse(suitcoralall$X2090==0,2090,1),suitcoralall$year2085)
  suitcoralall$year2095 <- ifelse(suitcoralall$year2090 == 1, ifelse(suitcoralall$X2095==0,2095,1),suitcoralall$year2090)
  suitcoralall$year2100 <- ifelse(suitcoralall$year2095 == 1, ifelse(suitcoralall$X2100==0,2100,1),suitcoralall$year2095)
  suitcoralall$expyr <- ifelse(suitcoralall$year2100 == 1, 2105,suitcoralall$year2100)
  
  suitcoralall$expyr[suitcoralall$expyr == 0]<- NA
  
  write.csv(suitcoralall,"final_yrexp.csv", row.names=FALSE)
}

#function to count the percentage of suitable sites
suitpercent  <- function(wdrcp) {
  setwd(wdrcp)
  
  #get count of suitability
  csvfiles <- list.files(wdrcp,pattern="final_",full.names=T)
  dt <- fread(csvfiles[1])
  countcoral<-data.table()
  years <- c( "X1855","X1865","X1875","X1885","X1895","X1905","X1915","X1925","X1935","X1945","X1955","X1965","X1975","X1985","X1995",
              "X2005.x","X2005.y","X2010","X2015","X2020","X2025","X2030","X2035","X2040","X2045","X2050","X2055","X2060","X2065","X2070","X2075",
              "X2080","X2085","X2090","X2095","X2100")
  for (i in 1:36){
    c <- count(dt, vars=years[i])
    colnames(c) <- c("exp", "freq")
    y <- add_column(c,year=years[i],.before=T)
    countcoral <- rbind(y,countcoral)
    
  } 
  
  write.csv(countcoral,"exp_count.csv", row.names=FALSE)
  
  
}

#function to make stressor count shapefile 
expiredstressors = function(wdrcp) {
  
  setwd(wdrcp)
  csvfiles <- list.files(wdrcp,pattern="stressor",full.names=T)
  suitcoralall <- fread(csvfiles[1],header=T)
  
  
  #merge with shapefile and save
  suitcoralyr <- merge(reefs,suitcoralall,by="ID")
  st_write(suitcoralyr, dsn = "exp_stressors", driver = "ESRI Shapefile")
  
}

expstressorpercent = function(wdrcp){
  
  setwd(wdrcp)

  csvfiles <- list.files(wdrcp,pattern="stressor",full.names=T)
  dt <- fread(csvfiles[2])
  dt$nexp2100 <- 5-dt$nstressors2100

    c <- count(dt, vars="nexp2100")

  write.csv(c,"stressor2100_count.csv", row.names=FALSE)
  
}
  

#rcp85
expyr(wdhist="dhw/historic",wdrcp="dhw/rcp85") 
expyr(wdhist="land_hurtt/historic",wdrcp="land_hurtt/rcp85")
expyr(wdhist="oa/hist",wdrcp="oa/rcp85")
expyr(wdhist="pop_hr/hist_hyde",wdrcp="pop_hr/ssp5")
expyrstorm(wdhist="storms/historic",wdrcp="storms/rcp85")

#rcp45
expyr(wdhist="dhw/historic",wdrcp="dhw/rcp45") 
expyr(wdhist="land_hurtt/historic",wdrcp="land_hurtt/rcp45")
expyr(wdhist="oa/hist",wdrcp="oa/rcp45")
expyr(wdhist="pop_hr/hist_hyde",wdrcp="pop_hr/ssp2")
expyrstorm(wdhist="storms/historic",wdrcp="storms/rcp45")

#rcp26
expyr(wdhist="dhw/historic",wdrcp="dhw/rcp26")
expyr(wdhist="land_hurtt/historic",wdrcp="land_hurtt/rcp26")
expyr(wdhist="oa/hist",wdrcp="oa/rcp26")
expyr(wdhist="pop_hr/hist_hyde",wdrcp="pop_hr/ssp1")
expyrstorm(wdhist="storms/historic",wdrcp="storms/rcp26")


#adjust expiration year calculation

#rcp85
expyrnew(wdrcp="dhw/rcp85") 
expyrnew(wdrcp="land_hurtt/rcp85")
expyrnew(wdrcp="oa/rcp85")
expyrnew(wdrcp="pop_hr/ssp5")
expyrnew(wdrcp="storms/rcp85")

#rcp45
expyrnew(wdrcp="dhw/rcp45") 
expyrnew(wdrcp="land_hurtt/rcp45")
expyrnew(wdrcp="oa/rcp45")
expyrnew(wdrcp="pop_hr/ssp2")
expyrnew(wdrcp="storms/rcp45")

#rcp26
expyrnew(wdrcp="dhw/rcp26")
expyrnew(wdrcp="land_hurtt/rcp26")
expyrnew(wdrcp="oa/rcp26")
expyrnew(wdrcp="pop_hr/ssp1")
expyrnew(wdrcp="storms/rcp26")


#adjust storms dataset so columns match the other stressors
stormfiles <- list.files("storms/rcp85",pattern="final",full.names=T)
stormfiles <- list.files("storms/rcp45",pattern="final",full.names=T)
stormfiles <- list.files("storms/rcp26",pattern="final",full.names=T)

storms <- fread(file=stormfiles[1],header=T)
storms[is.na(storms)] <- 1
storms <- storms[order(ID)]
#add empty rows for 1855-1950 for storms
storms[ , c("X1855","X1865", "X1875","X1885","X1895", "X1905", "X1915", "X1925", "X1935", "X1945")] <- 1
#remove unneeded columns
storms <- subset(storms, select=-c(X1950,X1960,X1970,X1980,X1990,X2000))
#reorder
setcolorder(storms,c("ID", "X1855","X1865","X1875","X1885","X1895","X1905","X1915","X1925","X1935","X1945","X1955","X1965","X1975","X1985","X1995",
                     "X2005.x","X2005.y","X2010","X2015","X2020","X2025","X2030","X2035","X2040","X2045","X2050","X2055","X2060","X2065","X2070","X2075",
                     "X2080","X2085","X2090","X2095","X2100"))
write.csv(storms,"storms/rcp26/final_aallyr.csv", row.names=FALSE)


##overall
expyroverallnew("suit/rcp85-ssp5","dhw/rcp85","land_hurtt/rcp85","oa/rcp85","pop_hr/ssp5","storms/rcp85")
expyroverallnew("suit/rcp45-ssp2","dhw/rcp45","land_hurtt/rcp45","oa/rcp45","pop_hr/ssp2","storms/rcp45")
expyroverallnew("suit/rcp26-ssp1","dhw/rcp26","land_hurtt/rcp26","oa/rcp26","pop_hr/ssp1","storms/rcp26")


#calculate percent expiration
#rcp85
suitpercent("suit/rcp85-ssp5")
suitpercent("dhw/rcp85")
suitpercent("land_hurtt/rcp85")
suitpercent("oa/rcp85")
suitpercent("pop_hr/ssp5")
suitpercent("storms/rcp85")

#rcp45
suitpercent("suit/rcp45-ssp2")
suitpercent("dhw/rcp45")
suitpercent("land_hurtt/rcp45")
suitpercent("oa/rcp45")
suitpercent("pop_hr/ssp2")
suitpercent("storms/rcp45")

#rcp26
suitpercent("suit/rcp26-ssp1")
suitpercent("dhw/rcp26")
suitpercent("land_hurtt/rcp26")
suitpercent("oa/rcp26")
suitpercent("pop_hr/ssp1")
suitpercent("storms/rcp26")

#make shapefiles for figures. then go to ArcPro to rasterize for figure

#rcp85
expirationyr(wdrcp="dhw/rcp85") 
expirationyr(wdrcp="oa/rcp85") 
expirationyr(wdrcp="sst/rcp85")
expirationyr(wdrcp="land_hurtt/rcp85")
expirationyr(wdrcp="pop_hr/ssp5")
expirationyr(wdrcp="storms/rcp85")

#rcp45
expirationyr(wdrcp="dhw/rcp45") 
expirationyr(wdrcp="oa/rcp45") 
expirationyr(wdrcp="sst/rcp45")
expirationyr(wdrcp="land_hurtt/rcp45")
expirationyr(wdrcp="pop_hr/ssp2")
expirationyr(wdrcp="storms/rcp45")

#rcp26
expirationyr(wdrcp="dhw/rcp26") 
expirationyr(wdrcp="oa/rcp26") 
expirationyr(wdrcp="sst/rcp26")
expirationyr(wdrcp="land_hurtt/rcp26")
expirationyr(wdrcp="pop_hr/ssp1")
expirationyr(wdrcp="storms/rcp26")

#overall
expirationyr(wdrcp="suit/rcp85-ssp5")
expirationyr(wdrcp="suit/rcp45-ssp2")
expirationyr(wdrcp="suit/rcp26-ssp1")

#stressor count
expiredstressors("suit/rcp85-ssp5")
expiredstressors("suit/rcp45-ssp2")
expiredstressors("suit/rcp26-ssp1")

#distribution of number of stressors
expstressorpercent("suit/rcp85-ssp5")
expstressorpercent("suit/rcp45-ssp2")
expstressorpercent("suit/rcp26-ssp1")


