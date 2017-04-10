#' Create HYSPLIT-derived Seasonal/Monthly and Composite Airshed
#' 
#' Create monthly/seasonal and composite airsheds from HYSPLIT back trajectories for a given point location
#' @usage airshed(loc, lon, lat, day_select=NA, month_select=NA, year_select,
#'    season_name="season",mode="seasonal", hy_alt, sPix=0.3,
#'    nPix=100, cutoff=10, composite=TRUE, ...)
#' @param loc character. Name of location, starting point, used in naming subdirectories and files
#' @param lon numeric. Longitude of starting point
#' @param lat numeric. Latitude of starting point
#' @param day_select numeric vector, integers, [1,365]. Days to be run, in Julian days, non-leap years, useful if mode="seasonal"
#' @param month_select numeric vector, integers, [1,12]. Months to be run, useful if mode="monthly"
#' @param year_select numeric vector, integers. Years to be run
#' @param season_name character. Name of season, used in the names of the output folder and file [default is "season"]
#' @param mode character: "monthly" or "seasonal" [default is "seasonal"]
#' @param hy_alt character string or vector. Name(s) of hysplit heights to be run
#' @param sPix numeric. Size of pixel, in degrees [default is 0.3]
#' @param nPix numeric. Number of pixels in horizontal and vertical direction [default is 100]
#' @param cutoff numeric. Threshold, in percent, for pixels to be considered part of seasonal airshed  [default to 10]
#' @param composite logical. Create composite airshed? [default is TRUE]
#' @param ... see global arguments: pointDir, ask_home
#' @keywords hysplit airshed seasonal monthly composite
#' @export
#' @examples 
#' #creates winter delhi airshed for 2007-2013 and composite winter airshed for these years
#' airshed(loc="delhi",lon=77.2090,lat=28.6139,month_select=10:11,year_select=2007:2013,season_name="win",hy_alt="500m")

airshed <- function(loc, lon, lat, day_select=NA, month_select=NA, year_select,
                    season_name="season",mode="seasonal", hy_alt, sPix=0.3,
                    nPix=100, cutoff=10, composite=TRUE,
                    pointDir=paste0("./",loc,"/"), ask_home=TRUE) {
  
  #Open packages
  library("raster");library("rgeos");library("maptools");
  library("mapproj");library("maps");library("rgdal");
  library("spdep");library("sp")
  
  options(error=NULL)
  
  #folder locations
  Home <- getwd()
  
  if (ask_home=="TRUE") {
    DirAns <- readline(paste("Is the project home directory:",Home,"? (enter y/n) "))
    if (DirAns!="y") {stop("Set project home directory with setwd()")}
  }
  
  if (mode=="seasonal") {
    if (is.na(day_select)==TRUE) {stop("set day_select")}
    if (length(day_select[which(day_select>365 | day_select<1)])>0) {stop("day_select must be integers [1,365]")}
  }
  if (mode=="monthly") {
    if (is.na(month_select)==TRUE) {stop("set month_select")}
    if (length(month_select[which(month_select>12 | month_select<1)])>0) {stop("month_select must be integers [1,12]")}
  }
  
  #INPUT: hysplit gis & txt files
  #subfolders = iYear (e.g. 2007) -> month (e.g. oct)
  trajHome <- paste0(pointDir,"hysplit/")
  ifelse(!file.exists(trajHome),stop("missing hysplit coordinates folder"),FALSE)
  
  if (mode=="seasonal") {
    outputHome <- paste0(pointDir,"airshed_",season_name,"/shapefiles/")
    ifelse(!file.exists(outputHome),dir.create(outputHome,recursive=T),FALSE)
  }
  
  if (mode=="monthly") {
    outputHome <- paste0(pointDir,"airshed_monthly/shapefiles/")
    ifelse(!file.exists(outputHome),dir.create(outputHome,recursive=T),FALSE)
  }
  
  print("Okay, wait for function to finish.")
  
  gt <- GridTopology(c(lon-(sPix*nPix/2),lat-(sPix*nPix/2)),c(sPix,sPix),c(nPix,nPix))
  grd <- SpatialGrid(gt,proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  spix <- as(grd,"SpatialPixels")
  spol <- as(spix,"SpatialPolygons")
  spol$coor <- coordinates(spol)
  
  for (iYear in year_select) {
    
    if (mode=="monthly") {
      if ((iYear %% 4)!=0) { #non-leap year
        st.jul = c(1,32,60,91,121,152,182,213,244,274,305,335)
        end.jul = c(st.jul[-1]-1,365)
      }
      
      if ((iYear %% 4)==0) { #leap year
        st.jul = c(1,32,61,92,122,153,183,214,245,275,306,336)
        end.jul = c(st.jul[-1]-1,366)
      }
      
      for (iMonth in month_select) {
        
        TrajName2 <- list()
        for (iAltitude in 1:length(hy_alt)) {
          iAltitudeName <- hy_alt[iAltitude]
          shape.lat <- read.table(paste0(trajHome,loc,"_HysplitLat_",iYear,"_",iAltitudeName,".csv"),sep=",")
          shape.lon <- read.table(paste0(trajHome,loc,"_HysplitLon_",iYear,"_",iAltitudeName,".csv"),sep=",")
          
          TrajName1 <- c()
          for(iDay in st.jul[iMonth]:end.jul[iMonth]) {
            if (length(shape.lat[iDay,][!is.na(shape.lat[iDay,])])>0) {
              TempName <- paste0("alt",iAltitude,"traj",iDay)
              TrajName1[iDay] <- TempName
              
              shape.coor <- cbind(as.numeric(shape.lat[iDay,]),as.numeric(shape.lon[iDay,]))
              shape.coor <- apply(shape.coor,2,na.omit)
              shape <- SpatialPoints(shape.coor,proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
              
              ras <- rasterize(as(shape,"SpatialLines"),raster(spix))
              
              assign(TempName,ras)
            }
          }
          TrajName1 <- TrajName1[complete.cases(TrajName1)]
          TrajName2[[iAltitude]] <- TrajName1
        }
        TrajNames <- do.call(c,TrajName2)
        airshedS <- sum(stack(mget(TrajNames)),na.rm=T)
        air_ext <- extract(airshedS,spix)
        spol$den <- air_ext/length(TrajNames)*100
        
        writeOGR(spol[which(spol$den>cutoff),],layer=paste0("airshed",iYear,"monthly"),dsn=outputHome,driver="ESRI Shapefile",overwrite=T)
      }
    }   
    
    if (mode=="seasonal") {
      
      #convert to julian days
      if ((iYear %% 4)!=0) { #non-leap year
        st.jul = c(1,32,60,91,121,152,182,213,244,274,305,335)
        end.jul = c(st.jul[-1]-1,365)
        DaysInMon = c(31,28,31,30,31,30,31,31,30,31,30,31)
        day_select.adj <- day_select
      }
      if ((iYear %% 4)==0) { #leap year
        st.jul = c(1,32,61,92,122,153,183,214,245,275,306,336)
        end.jul = c(st.jul[-1]-1,366)
        DaysInMon = c(31,29,31,30,31,30,31,31,30,31,30,31)
        
        day_select.adj <- day_select
        day_select.adj <- c(day_select.adj[which(day_select.adj<28)],day_select.adj[which(day_select.adj>=28)]+1)
        if (day_select[1]==28) {day_select.adj <- c(28,day_select.adj)}
      }
      
      TrajName2 <- list()
      for (iAltitude in 1:length(hy_alt)) {
        iAltitudeName <- hy_alt[iAltitude]
        shape.lat <- read.table(paste0(trajHome,loc,"_HysplitLat_",iYear,"_",iAltitudeName,".csv"),sep=",")
        shape.lon <- read.table(paste0(trajHome,loc,"_HysplitLon_",iYear,"_",iAltitudeName,".csv"),sep=",")
        
        TrajName1 <- c()
        for(iDay in day_select.adj) {
          if (length(shape.lat[iDay,][!is.na(shape.lat[iDay,])])>0) {
            TempName <- paste0("alt",iAltitude,"traj",iDay)
            TrajName1[iDay] <- TempName
            
            shape.coor <- cbind(as.numeric(shape.lat[iDay,]),as.numeric(shape.lon[iDay,]))
            shape.coor <- apply(shape.coor,2,na.omit)
            shape <- SpatialPoints(shape.coor,proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
            
            ras <- rasterize(as(shape,"SpatialLines"),raster(spix))
            
            assign(TempName,ras)
          }
        }
        TrajName1 <- TrajName1[complete.cases(TrajName1)]
        TrajName2[[iAltitude]] <- TrajName1
      }
      TrajNames <- do.call(c,TrajName2)
      airshedS <- sum(stack(mget(TrajNames)),na.rm=T)
      air_ext <- extract(airshedS,spix)
      spol$den <- air_ext/length(TrajNames)*100
      
      writeOGR(spol[which(spol$den>cutoff),],layer=paste0("airshed",iYear,season_name),dsn=outputHome,driver="ESRI Shapefile",overwrite=T)
      
      if (composite==TRUE) {
        airshedS <- airshedS/length(TrajNames)*100
        airshedYr <- airshedS
        airshedYr[airshedYr<=cutoff] <- NA
        airshedYr[airshedYr>cutoff] <- 1
        
        assign(paste0("airshed",iYear),airshedS)
        assign(paste0("airshed",iYear,"_Yr"),airshedYr)
      }
    }
    
    
    cat(iYear,": year",which(year_select==iYear),"of",length(year_select),as.character(Sys.time()),'\n')
  }
  
  #Composite
  if (mode=="seasonal") {
    if (length(year_select)>1) {
      if (composite==TRUE) {
        all_air <- paste0("airshed",year_select)
        airshedC <- mean(stack(mget(all_air)),na.rm=T)
        air_ext <- extract(airshedC,spix)
        spol$den <- air_ext
        
        writeOGR(spol[which(spol$den>cutoff),],layer=paste0("airshed_",season_name,"_composite"),dsn=outputHome,driver="ESRI Shapefile",overwrite=T)
        
        all_air <- paste0("airshed",year_select,"_Yr")
        airshedC <- sum(stack(mget(all_air)),na.rm=T)
        air_ext <- extract(airshedC,spix)
        spol$den <- air_ext
        
        writeOGR(spol[which(spol$den>0),],layer=paste0("airshed_",season_name,"_compositeYr"),dsn=outputHome,driver="ESRI Shapefile",overwrite=T)
      }
    }
  }
  
  print("Complete!")
}