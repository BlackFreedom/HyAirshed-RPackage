#' Extract FRP from Monthly/Seasonal HYSPLIT Airshed
#' 
#' Extract the monthly/seasonal FRP from HYSPLIT-derived airsheds of selected point location(s)
#' @usage modeHyFRP(loc, month_select=NA, year_select,
#'    season_name="season", mode="seasonal", frp_mode="combined", ...)
#' @param loc character. Name of location, starting point, used in naming subdirectories and files
#' @param month_select numeric vector. Months to be run, only useful when mode="monthly" [default is NA]
#' @param year_select numeric vector. Years to be run
#' @param season_name character. Name of season, used in the name of the output table when mode="seasonal" [default is "season"]
#' @param mode character: "monthly" or "seasonal". [default is "seasonal"]
#' @param frp_mode character: "terra", "aqua", or "combined". Choose Terra-only, Aqua-only, or Terra+Aqua FRP. [default is "combined"]
#' @param ... see global arguments: pointDir, ask_home
#' @keywords hysplit airshed monthly seasonal frp modis
#' @export
#' @examples
#' #extracts Oct-Nov FRP from Delhi winter airshed in 2007-2013, FRP already stacked by modeFRPstack()
#' modeHyFRP(loc="delhi", year_select=2007:2013, season_name="win")
#' 
#' #extracts monthly FRP from Delhi monthly airshed in 2007-2013, FRP already stacked by modeFRPstack()
#' modeHyFRP(loc="delhi", month_select=1:12, year_select=2007:2013, mode="monthly")


modeHyFRP <- function(loc, month_select=NA, year_select, season_name="season",
                      mode="seasonal", frp_mode="combined", pointDir=paste0("./",loc,"/"),
                      ask_home=TRUE) {

  options(error=NULL)
  
  #folder locations
  Home <- getwd()
  
  if (ask_home=="TRUE") {
    DirAns <- readline(paste("Is the project home directory:",Home,"? (enter y/n) "))
    if (DirAns!="y") {stop("Set project home directory with setwd()")}
  }
  
  if (mode=="monthly") {
    if (is.na(month_select)==TRUE) {
      stop("Set months (month_select argument)")
    }
  }
  
  if (frp_mode=="terra") {nameAdd <- "_terra"}
  if (frp_mode=="aqua") {nameAdd <- "_aqua"}
  if (frp_mode=="combined") {nameAdd <- ""}
  
  #INPUT: point location airshed directory
  pointDir <- paste0("./",loc,"/")
  ifelse(!file.exists(pointDir),stop("missing directory for point location"),FALSE)
  
  #OUTPUT: frp in airshed
  outputHome <- paste0(pointDir,"frp_airshed",nameAdd,"/")
  
  #create output folders
  ifelse(!file.exists(outputHome),dir.create(outputHome),FALSE)
  
  print("Okay, wait for function to finish.")
  
  if (mode=="monthly") {
    
    rasDir <- "./rasters/"
    ifelse(!file.exists(rasDir),dir.create(rasDir),FALSE)
    
    #INPUT: airshed
    airDir <- paste0(pointDir,"airshed_monthly/shapefiles/")
    ifelse(!file.exists(airDir),stop("missing airshed folder"),FALSE)
    
    #INPUT: frp monthly rasters
    frpDir <- paste0(rasDir,"frp_monthly",nameAdd,"/")
    
    #create output folders
    ifelse(!file.exists(outputHome),dir.create(outputHome),FALSE)
    
    collectYear <- list()
    for (iYear in year_select) {
      collectMonth <- rep(NA,12)
      for (iMonth in month_select) {
        
        frp <- raster(paste0(frpDir,iYear,"/FRP",iYear*100+iMonth,"monthly.tif"))
        
        if (file.exists(paste0(airDir,iYear,"/airshed",iYear*100+iMonth,".shp"))==TRUE) {
          airshed <- readShapePoly(paste0(airDir,iYear,"/airshed",iYear*100+iMonth),proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
          airshed <- spTransform(gUnaryUnion(airshed),crs(frp))
          
          frp_extract <- extract(frp,airshed,sum,na.rm=T)*0.1
          collectMonth[iMonth] <- frp_extract
        }
      }
      collectYear[[iYear-(year_select[1]-1)]] <- collectMonth
      cat(iYear,as.character(Sys.time()),'\n') #completion time by iYear in loop
    }
    collectTable <- do.call(rbind,collectYear)
    rownames(collectTable) <- year_select[1]:year_select[length(year_select)]
    colnames(collectTable) <- tolower(month.abb)
    
    #save file in a table
    write.table(collectTable,paste0(outputHome,loc,"_monthly_airshed_frp.csv"),sep=",")
  }
  
  if (mode=="seasonal") {
    
    #INPUT: airshed
    airDir <- paste0(pointDir,"airshed_",season_name,"/shapefiles/")
    ifelse(!file.exists(airDir),stop("missing airshed folder"),FALSE)
    
    #INPUT: frp seasonal rasters
    frpDir <- paste0(pointDir,"rasters/frp_",season_name,nameAdd,"/")
    ifelse(!file.exists(frpDir),stop("missing frp directory"),FALSE)
    
    collectSeason <- c()
    for (iYear in year_select) {
      frp <- raster(paste0(frpDir,"/FRP",iYear,season_name,".tif"))
      
      airshed <- readShapePoly(paste0(airDir,"airshed",iYear,season_name),proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
      airshed <- spTransform(gUnaryUnion(airshed),crs(frp))
      
      frp_extract <- extract(frp,airshed,sum,na.rm=T)*0.1
      collectSeason[iYear-(year_select[1]-1)] <- frp_extract
      cat(iYear,as.character(Sys.time()),'\n') #completion time by iYear in loop
    }
    
    collectTable <- as.data.frame(collectSeason)
    rownames(collectTable) <- year_select[1]:year_select[length(year_select)]
    colnames(collectTable) <- paste0("frp_",season_name)
    
    #save file in a table
    write.table(collectTable,paste0(outputHome,loc,"_",season_name,"_airshed_frp.csv"),sep=",")
  }
  
  print("Complete!")
}
