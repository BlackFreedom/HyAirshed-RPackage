#' Sum and Stack MODIS FRP to Daily
#'
#' Stack and sum Terra and Aqua daily FRP & mask ocean (for HYSPLIT). Input takes processed FRP .tif files the MODIS Reprojection Tool.
#' @usage ... dailyFRPstack(year_select, sDay, eDay,
#'    frp_mode="combined", ocean_mask=TRUE, ...)
#' @param year_select numeric vector. Years to be run
#' @param sDay numeric. First day of starting 8-day composite
#' @param eDay numeric. First day of ending 8-day composite
#' @param frp_mode character: "terra", "aqua", or "combined". Choose Terra-only, Aqua-only, or Terra+Aqua FRP. [default is "combined"]
#' @param ocean_mask logical. Mask ocean? [default is TRUE]
#' @param ... see global arguments: rawDir, rasDir, ask_home
#' @keywords frp modis daily stack
#' @export
#' @examples 
#' #renames and stacks FRP in the years 2007-2013 for all days and masks ocean
#' dailyFRPstack(year_select=2007:2013,sDay=1,eDay=361)
#' 
#' #renames terra FRP in the years 2007-2013 for all days and masks ocean
#' dailyFRPstack(year_select=2007:2013,sDay=1,eDay=361,mode="terra")

dailyFRPstack <- function(year_select, sDay, eDay, frp_mode="combined",
                          ocean_mask=TRUE, rawDir="./raw_files/", 
                          rasDir="./rasters/", ask_home=TRUE) {
  
  options(error=NULL)
  
  #folder locations
  Home <- getwd()
  
  if (ask_home=="TRUE") {
    DirAns <- readline(paste("Is the project home directory:",Home,"? (enter y/n) "))
    if (DirAns!="y") {stop("Set project home directory with setwd()")}
  }
  
  print("Okay, wait for function to finish.")
  
  #open ocean_mask
  if (ocean_mask==TRUE) {ocean <- raster("rasters/ocean_mask.tif")} #ocean mask
  days <- seq(1,361,8) #start days of 8-day composites
  
  if(length(which(days==sDay))==0) {
    stop(paste("sDay is invalid, must be one of these:",paste(days,collapse = ",")))
  }
  if(length(which(days==eDay))==0) {
    stop(paste("eDay is invalid, must be one of these:",paste(days,collapse = ",")))
  }
  
  for (iYear in year_select) {
    
    #INPUT: daily frp .tif files processed by modis reprojection tool;
    #subfolders = iYear (e.g. 2007)
    rawHome <- paste0(rawDir,iYear,"/")
    ifelse(!file.exists(rawHome),stop("missing raw files/folder"),FALSE)
    
    #OUTPUT: daily frp rasters
    #subfolders = iYear (e.g. 2007)
    
    if (frp_mode=="terra") {nameAdd <- "_terra"}
    if (frp_mode=="aqua") {nameAdd <- "_aqua"}
    if (frp_mode=="combined") {nameAdd <- ""}
    
    frpHome <- paste0(rasDir,"frp_daily",nameAdd,"/")
    frpSub <- paste0(frpHome,iYear,"/")

    ifelse(!file.exists(rasDir),dir.create(rasDir),FALSE)
    ifelse(!file.exists(frpHome),dir.create(frpHome),FALSE)
    ifelse(!file.exists(frpSub),dir.create(frpSub),FALSE)
    
    for(iDay in which(days==sDay):which(days==eDay)) {
      
      #STACK
      if (days[iDay] < 361) {
        TComp = 8}
      
      if (days[iDay] == 361) {
        if(iYear%%4>0) {TComp = 5}
        if(iYear%%4==0) {TComp = 6}
      }
      
      for (iComp in 1:TComp) {
        terra <- raster(paste0(rawHome,"MOD14A1.A",iYear*1000+days[iDay],".MaxFRP.Number_of_Days_0",iComp,".tif")) #names of objects to stack
        aqua <- raster(paste0(rawHome,"MYD14A1.A",iYear*1000+days[iDay],".MaxFRP.Number_of_Days_0",iComp,".tif")) #names of objects to stack
        
        if (frp_mode=="terra") {frp <- terra}
        if (frp_mode=="aqua") {frp <- aqua}
        if (frp_mode=="combined") {frp <- sum(stack(c(aqua,terra)),na.rm=T)} #stack terra/aqua daily
        
        if (ocean_mask==TRUE) {stack <- frp*ocean} else {stack <- frp} #mask ocean if ocean_mask=TRUE
        
        writeRaster(stack,paste0(frpHome,"FRP",iYear*1000+days[iDay]+(iComp-1)),format="GTiff",overwrite=T)
        }
      
      cat(iYear,as.character(Sys.time()),'\n') #completion time by iYear in loop
    }
  }
  print("Complete!")
}