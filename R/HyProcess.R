#' Process HYSPLIT files and save coordinates to table
#' 
#' Process HYSPLIT files and save coordinates to table
#' @usage HyProcess(loc, year_select, hy_alt, nHours, alt.adj=0, ...)
#' @param loc character. Name of location, starting point, used in naming subdirectories and files
#' @param year_select numeric vector. Years to be run
#' @param hy_alt character string or vector. Name(s) of hysplit heights to be run
#' @param nHours numeric. Total hours in back trajectories
#' @param alt.adj numeric. Adjustment if more than one location [default is 0]
#' @param ... see global arguments: hyDir, pointDir, ask_home
#' @keywords hysplit process
#' @export
#' @examples 
#' #Process all available HYSPLIT downloads and organizes the extracted coordinates of the back trajectories into tables 
#' HyProcess(loc="delhi",year_select=2007:2013,hy_alt="0.5km",nHours=72)

HyProcess <- function(loc, year_select, hy_alt, nHours, alt.adj=0, 
                      hyDir="./hysplit/", pointDir=paste0("./",loc,"/"),ask_home="TRUE") {
  
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
  
  #INPUT: hysplit gis & txt files
  #subfolders = iYear (e.g. 2007) -> month (e.g. oct)
  trajHome <- paste0(hyDir,"hysplit_gis/")
  ifelse(!file.exists(trajHome),stop("missing hysplit gis directory"),FALSE)
  
  #INPUT: hysplit pdf files
  #subfolders = iYear (e.g. 2007) -> month (e.g. oct)
  hypdfHome <- "./hysplit/hysplit_plot/"
  ifelse(!file.exists(trajHome),stop("missing hysplit pdf plots directory"),FALSE)
  
  #OUTPUT: frp in daily airshed
  #subfolders = 500m, 1000m, 1500m
  outputHome <- paste0(pointDir,"hysplit/")
  
  #create output directory, folders
  ifelse(!file.exists(pointDir),dir.create(pointDir),FALSE)
  ifelse(!file.exists(outputHome),dir.create(outputHome),FALSE)
  
  ifelse(!file.exists("hysplitDays.csv"),stop("hysplitDays.csv"),FALSE)
  availableDays <- read.table("hysplitDays.csv",sep=",",header=T)
  
  print("Okay, wait for function to finish.")
  
  for(iAltitude in 1:length(hy_alt)) {
    iAltitudeName <- hy_alt[iAltitude] #get altitude
    st.shp <- 1000*iAltitude+alt.adj; end.shp <- st.shp+1000 #select trajectories by altitude (useful if more than one altitudes or cities downloaded)
    
    #select iYear(s)
    for(iYear in year_select) {
      
      #convert to julian days
      if ((iYear %% 4)!=0) { #non-leap year
        st.jul = c(1,32,60,91,121,152,182,213,244,274,305,335)
        end.jul = c(st.jul[-1]-1,365)
        DaysInMon = c(31,28,31,30,31,30,31,31,30,31,30,31)
      }
      if ((iYear %% 4)==0) { #leap year
        st.jul = c(1,32,61,92,122,153,183,214,245,275,306,336)
        end.jul = c(st.jul[-1]-1,366)
        DaysInMon = c(31,29,31,30,31,30,31,31,30,31,30,31)
      }
      
      month.lower <- tolower(month.abb)
      month_folders <- dir(paste0(hypdfHome,iYear,"/"))
      
      month.id <- c()
      for (iFolder in 1:length(month_folders)) {
        month.id[iFolder] <- which(month.lower==month_folders[iFolder])
      }
      month_select <- sort(month.id)
      
      collectLat <- do.call(rbind,rep(list(rep(NA,nHours+1)),end.jul[12]));
      collectLon <- do.call(rbind,rep(list(rep(NA,nHours+1)),end.jul[12]));
      
      for (iMonth in 1:length(month_select)) {
        mon.actual <- month_select[iMonth] #month number
        mon.name <- month.lower[mon.actual] #month number, starting from 1
        iMday <- eval(parse(text=as.character(availableDays[which(availableDays$X==iYear),which(names(availableDays)==mon.name)]))) #days of month
        Jdays <- c(st.jul[mon.actual]:end.jul[mon.actual])[iMday] #convert days to Julian
        
        hy <- dir(paste0(hypdfHome,iYear,"/",mon.name),recursive=TRUE,pattern="*\\.pdf$",full.names = T) #pdf file
        details <- file.info(hy) #meta data of pdf file
        details <- details[with(details, order(as.POSIXct(mtime))), ] #order by time pdf file was modified
        pdffiles <- rownames(details) #get names of pdf files
        
        #get number of hysplit file to make a list of ordered hysplit .shp files
        temp <- as.character(matrix(unlist(strsplit(pdffiles,split="_")), ncol=3, byrow=TRUE)[,3])
        temp <- as.character(matrix(unlist(strsplit(temp,split=".pdf"))))
        gisfiles <- paste0("gis_",temp,"/GIS_traj01_",temp,".shp")
        
        for (iDay in 1:length(iMday)) {
          emptyPoints <- rep(NA,nHours+1)
          
          #open hysplit table
          shape <- readShapePoints(paste0(trajHome,iYear,"/",mon.name,"/",gisfiles[iDay]))
          #select subset trajectory within gis file
          shape <- shape[which(shape$id >= st.shp & shape$id < end.shp),]
          
          collectLat[Jdays[iDay],] <- replace(collectLat[Jdays[iDay],],1:length(shape),as.vector(coordinates(shape)[,1]))
          collectLon[Jdays[iDay],] <- replace(collectLon[Jdays[iDay],],1:length(shape),as.vector(coordinates(shape)[,2]))
        } 
      }
      
      write.table(collectLat,file=paste0(outputHome,loc,"_HysplitLat_",iYear,"_",iAltitudeName,".csv"),sep=",",col.names=paste0("point",1:(nHours+1)),row.names=paste0("coor",1:end.jul[12]))
      write.table(collectLon,file=paste0(outputHome,loc,"_HysplitLon_",iYear,"_",iAltitudeName,".csv"),sep=",",col.names=paste0("point",1:(nHours+1)),row.names=paste0("coor",1:end.jul[12]))
      
    }
    cat("altitude",iAltitude,'of',length(hy_alt),as.character(Sys.time()),'\n')
  }
  print("Complete!")
}
