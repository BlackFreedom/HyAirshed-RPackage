#' Extract FRP from Daily HYSPLIT Airshed
#' 
#' Extract the daily FRP from HYSPLIT-derived airsheds of selected point location(s)
#' @usage dailyHyFRP(loc, day_select, year_select, buf=10000,
#'    hy_alt, nHours, local_time, frp_mode="combined", ...)
#' @param loc character. Name of location, starting point, used in naming subdirectories and files
#' @param day_select numeric vector, integers, [1,365]. Days to be run, in julian days, non-leap years
#' @param year_select numeric vector, integers. Years to be run
#' @param buf numeric. Width of buffer, in meters, from original trajectory [default is 10000]
#' @param hy_alt character string or vector. Name(s) of hysplit heights to be run
#' @param nHours integer. Total hours in back trajectories.
#' @param local_time numeric. local time of HYSPLIT starting time, in hours
#' @param frp_mode character: "terra", "aqua", or "combined". Choose Terra-only, Aqua-only, or Terra+Aqua FRP. [default is "combined"]
#' @param ... see global arguments: rasDir, pointDir, ask_home
#' @keywords hysplit airshed daily frp modis
#' @export
#' @examples 
#' #extracts daily FRP in October and November, 2007-2013
#' #from a 10-km buffer (20 km in width) of 72-hour HYSPLIT
#' #back trajectories started from Delhi at 5:00pm (17:00) local time
#' dailyHyFRP(loc="delhi",day_select=274:334,year_select=2007:2013,buf=10000,hy_alt="500m",local_time=17)

dailyHyFRP <- function(loc, day_select, year_select, buf=10000,
                       hy_alt, nHours, local_time, frp_mode="combined",
                       rasDir="./rasters/", pointDir=paste0("./",loc,"/"),
                       ask_home=TRUE) {
  
  #open packages (install first if not in "Packages")
  library("raster");library("rgeos");library("maptools");
  library("mapproj");library("maps");library("rgdal");
  library("spdep");library("sp");library("fields")
  
  options(error=NULL)
  
  #folder locations
  Home <- getwd()
  
  if (ask_home=="TRUE") {
    DirAns <- readline(paste("Is the project home directory:",Home,"? (enter y/n) "))
    if (DirAns!="y") {stop("Set project home directory with setwd()")}
  }
  
  if (length(day_select[which(day_select>365 | day_select<1)])>0) {stop("day_select must be integers [1,365]")}
  
  #INPUT: hysplit gis & txt files
  #subfolders = iYear (e.g. 2007) -> month (e.g. oct)
  trajHome <- paste0(pointDir,"hysplit/")
  ifelse(!file.exists(trajHome),stop("missing hysplit coordinates folder"),FALSE)
  
  #INPUT: frp daily stack rasters, ocean masked out
  #subfolders = iYear (e.g. 2007)
  if (frp_mode=="terra") {nameAdd <- "_terra"}
  if (frp_mode=="aqua") {nameAdd <- "_aqua"}
  if (frp_mode=="combined") {nameAdd <- ""}
  
  frpHome <- paste0(rasDir,"frp_daily",nameAdd,"/")
  ifelse(!file.exists(frpHome),stop("missing frp directory"),FALSE)
  
  #OUTPUT: frp in daily airshed
  #subfolders = 500m, 1000m, 1500m
  outputHome <- paste0(pointDir,"frp_airshed",nameAdd,"/")
  
  #create output directory, folders
  ifelse(!file.exists(outputHome),dir.create(outputHome),FALSE)
  
  print("Okay, wait for function to finish.")
  
  for (iAltitude in 1:length(hy_alt)) {
    iAltitudeName <- hy_alt[iAltitude] #get altitude
    hr.s <- c(1,seq(1+local_time,nHours,24)); hr.e <- c(hr.s[-1],nHours+1) #start/end time of hysplit for each day to split 72-hour trajectories by day
    day.div <- length(hr.s) #number of passes to partition trajectory
    
    #select iYear(s)
    for(iYear in year_select) {
      
      #convert to julian days
      if ((iYear %% 4)!=0) { #non-leap year
        nDays <- 365
        st.jul = c(1,32,60,91,121,152,182,213,244,274,305,335)
        end.jul = c(st.jul[-1]-1,nDays)
        
        DaysInMon <- c(31,28,31,30,31,30,31,31,30,31,30,31)
        DaysInAll = c(1:31,1:28,1:31,1:30,1:31,1:30,1:31,1:31,1:30,1:31,1:30,1:31)
        day_select.adj <- day_select
      }
      if ((iYear %% 4)==0) { #leap year
        nDays <- 366
        st.jul = c(1,32,61,92,122,153,183,214,245,275,306,336)
        end.jul = c(st.jul[-1]-1,nDays)
        
        DaysInMon = c(31,29,31,30,31,30,31,31,30,31,30,31)
        DaysInAll = c(1:31,1:29,1:31,1:30,1:31,1:30,1:31,1:31,1:30,1:31,1:30,1:31)
        
        day_select.adj <- day_select
        day_select.adj <- c(day_select.adj[which(day_select.adj<28)],day_select.adj[which(day_select.adj>=28)]+1)
        if (day_select[1]==28) {day_select.adj <- c(28,day_select.adj)}
      }
      
      if (file.exists(paste0(outputHome,iAltitudeName,"/",loc,"_daily",iYear,".csv"))=="TRUE") {
        collectFRP <- read.table(paste0(outputHome,iAltitudeName,"/",loc,"_daily",iYear,".csv"),sep=",")
        if (length(collectFRP[,1])!=nDays) {collectFRP <- do.call(rbind,rep(list(rep(NA,4)),nDays))}
      } else {collectFRP <- do.call(rbind,rep(list(rep(NA,4)),nDays))}
      
      for (iDay in day_select.adj) {
        #open hysplit coordinates
        shape.lat <- read.table(paste0(trajHome,loc,"_HysplitLat_",iYear,"_",iAltitudeName,".csv"),sep=",")
        shape.lon <- read.table(paste0(trajHome,loc,"_HysplitLon_",iYear,"_",iAltitudeName,".csv"),sep=",")
        
        if (length(shape.lat[iDay,][!is.na(shape.lat[iDay,])])>0) {
          #subset coordinates of day in loop
          shape.coor <- cbind(as.numeric(shape.lat[iDay,]),as.numeric(shape.lon[iDay,]))
          shape.coor <- apply(shape.coor,2,na.omit)
          shape <- SpatialPoints(shape.coor,proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
          
          #reproject trajectory to Mercator to match frp
          shape <- spTransform(shape,CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
          
          for (iPass in 1:length(hr.s)) {
            
            if (iPass>1 & iPass<=4) {
              #if trajectory ends before iPass starts, set frp as 0
              if (length(shape) <= hr.e[iPass-1]) {
                point <- 0
                collectFRP[iDay,iPass] <- 0}
            }
            
            if (iPass==1) {
              #subset hysplit trajectories by day
              if (length(shape) < hr.e[iPass] & length(shape) > 0) {
                point <- shape[hr.s[iPass]:length(shape),]
              }
            }
            
            if (iPass>1 & iPass<=4) {  
              #subset hysplit trajectories by day
              if (length(shape) < hr.e[iPass] & length(shape) > hr.e[iPass-1]) {
                point <- shape[hr.s[iPass]:length(shape),]
              }
            }
            
            if (length(shape) >= hr.e[iPass]) {
              point <- shape[hr.s[iPass]:hr.e[iPass],]
            }
            
            if (class(point)=="SpatialPoints") {
              #extract frp from hysplit trajectory
              traj <- as(point,"SpatialLines")
              traj$id <- 1
              
              traj_buf <- gBuffer(traj,byid=T,width=buf)
              traj_buf <- spTransform(traj_buf,CRS("+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
              traj <- spTransform(traj,CRS("+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
              
              #open frp files
              #open frp raster and multiply by 0.1 factor
              
              if (iDay-(iPass-1)<=0) {
                frp <- raster(paste0(frpHome,iYear,"/","FRP",(iYear-1)*1000+iDay-(iPass-1)+nDays,".tif"))*0.1
              } else {frp <- raster(paste0(frpHome,iYear,"/","FRP",iYear*1000+iDay-(iPass-1),".tif"))*0.1}
              collectFRP[iDay,iPass] <- cellStats(mask(frp,traj_buf),sum) #extract frp from buffer
            }
          }
        }
      }
      #save frp tables that compiles daily airshed frp of all months, by altitude
      ifelse(!file.exists(paste0(outputHome,iAltitudeName)),dir.create(paste0(outputHome,iAltitudeName)),FALSE)
      
      collectFRP <- cbind(collectFRP,rowSums(collectFRP,na.rm=T))
      write.table(collectFRP,file=paste0(outputHome,iAltitudeName,"/",loc,"_daily",iYear,".csv"),sep=",",row.names=paste(rep(tolower(month.abb),DaysInMon),DaysInAll),col.names=c(c("Same Day",paste(1:(day.div-1),"Days Ago")),"Total"))
      
      cat(iYear,": year",iYear-(iYear[1]-1),"of",length(year_select),";",hy_alt[iAltitude], ": altitude",iAltitude,"of",length(hy_alt),as.character(Sys.time()),'\n')
    }
  }
  
  #collect total frp of all altitudes by iYear and save as table
  for (iYear in year_select) {
    
    altTableAll <- list()
    for(iAltitude in 1:length(hy_alt)) {
      #copy for each altitude
      altTable <- read.table(paste0(outputHome,hy_alt[iAltitude],"/",loc,"_daily",iYear,".csv"),header=T,sep=",")
      #collect each altitude
      altTableAll[[iAltitude]] <- altTable$Total
    }
    altTableAll <- do.call(cbind,altTableAll)
    
    #save file in a table
    write.table(altTableAll,paste0(outputHome,loc,"_daily_airshed_frp",iYear,".csv"),sep=",",col.names=hy_alt,row.names=rownames(altTable))
  }
  print("Complete!")
}
