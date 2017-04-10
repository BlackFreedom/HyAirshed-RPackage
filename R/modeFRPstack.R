#' Sum and Stack MODIS Daily FRP to Monthly or Seasonal
#'
#' Stack and sum daily FRP to monthly and seasonal, up to a two-year period (for HYSPLIT). Input takes dailyFRPstack() output.
#' @usage modeFRPstack(loc=NA, day_select=NA, month_select=NA,
#'    year_select, mode="seasonal", frp_mode="combined",
#'    season_name="season", next_yr=NA, ...)
#' @param loc character. Name of location, starting point, used in naming subdirectories and files when mode="seasonal" [default is NA]
#' @param day_select numeric vector, integers, [1,365]. Days to be run, in Julian days, non-leap years, useful if mode="seasonal"
#' @param month_select numeric vector, integers, [1,12]. Months to be run, useful if mode=="monthly"
#' @param year_select numeric vector. Years to be run
#' @param mode character: "monthly" or "seasonal".  [default is "seasonal"]
#' @param frp_mode character: "terra", "aqua", or "combined". Choose Terra-only, Aqua-only, or Terra+Aqua FRP. [default is "combined"]
#' @param season_name character. Name of season, used in the names of the output subfolder and file [default is “season"]
#' @param next_yr numeric [1,365]. Days in next year, format is Julian days in non-leap year, useful for stacking days in winter months, only for "seasonal" mode [default is NA]
#' @param ... see global arguments: rasDir, pointDir, ask_home
#' @keywords frp modis monthly seasonal stack
#' @export
#' @examples
#' #stacks daily FRP to monthly for January-April from 2007-2013
#' modeFRPstack(year_select=2007:2013,month_select=1:4,mode="monthly")
#'  
#' #stacks daily FRP for the days 152-243 (leap:152-244) (JJA) in 2007-2013
#' modeFRPstack(year_select=2007:2013,mode="seasonal",day_select=152:243,season_name="summer")
#' 
#' #stacks daily FRP for the days 335-365 (leap:336-366) plus days 1-59 (leap:1-60) of the next year,
#' or DJF, in 2007-2013 (requires January FRP in 2014)
#' modeFRPstack(year_select=2007:2013,mode="seasonal",day_select=274:300,season_name="winter",next_yr=1:31)
#' 
#' #stacks daily FRP for Feburary 25-29, 2008
#' #in leap years, a day_select ending on day 59
#' #has day 60 added to the stack
#' modeFRPstack(year_select=2008,mode="seasonal",day_select=56:59,season_name="EndofFeb")

modeFRPstack <- function(loc=NA, day_select=NA, month_select=NA, year_select, mode="seasonal",
                         frp_mode="combined", season_name="season", next_yr=NA,
                         rasDir="./rasters/", pointDir=paste0("./",loc,"/"), ask_home=TRUE) {
  
  options(error=NULL)
  
  #folder locations
  Home <- getwd()
  
  if (ask_home=="TRUE") {
    DirAns <- readline(paste("Is the project home directory:",Home,"? (enter y/n) "))
    if (DirAns!="y") {stop("Set project home directory with setwd()")}
  }
  
  if (mode=="seasonal") {
    if (is.na(loc)==TRUE) {stop("set location (loc argument)")}
    if (is.na(day_select)==TRUE) {stop("set days (day_select argument)")}
    if (length(day_select[which(day_select>365 | day_select<1)])>0) {stop("day_select must be integers [1,365]")}
  }
  if (mode=="monthly") {
    if (is.na(month_select)==TRUE) {stop("set months (month_select argument)")}
    if (length(month_select[which(month_select>12 | month_select<1)])>0) {stop("month_select must be integers [1,12]")}
  }
  
  #INPUT: frp
  if (frp_mode=="terra") {nameAdd <- "_terra"}
  if (frp_mode=="aqua") {nameAdd <- "_aqua"}
  if (frp_mode=="combined") {nameAdd <- ""}
  
  ifelse(!file.exists(rasDir),stop("missing rasters folder with daily FRP"),FALSE)
  
  frpDir <- paste0(rasDir,"frp_daily",nameAdd,"/")
  ifelse(!file.exists(frpDir),dir.create(frpDir),FALSE)
  
  print("Okay, wait for function to finish.")
  
  for (iYear in year_select) {
    
    #INPUT: daily frp rasters
    #subfolders = iYear (e.g. 2007)
    frpHome <- paste0(frpDir,iYear,"/")
    ifelse(!file.exists(frpHome),stop("missing daily frp files/folder"),FALSE)
    
    if (mode=="monthly") {
      
      #OUTPUT: monthly frp rasters
      outfrpHome <- paste0(rasDir,"frp_monthly",nameAdd,"/")
      ifelse(!file.exists(outfrpHome),dir.create(outfrpHome),FALSE)
      
      outfrpSub <- paste0(outfrpHome,iYear,"/")
      ifelse(!file.exists(outfrpSub),dir.create(outfrpSub),FALSE)
      
      if ((iYear %% 4)!=0) { #non-leap year
        st.jul = c(1,32,60,91,121,152,182,213,244,274,305,335)
        end.jul = c(st.jul[-1]-1,365)
      }
      if ((iYear %% 4)==0) { #leap year
        st.jul = c(1,32,61,91,122,153,183,214,245,275,306,336)
        end.jul = c(st.jul[-1]-1,366)
      }
      
      for (iMonth in month_select) {
        #STACK
        yrjul = (iYear*1000+st.jul[iMonth]):(iYear*1000+end.jul[iMonth])
        frpnames <- paste0(frpHome,"FRP",yrjul,".tif")
        frp <- sum(stack(frpnames),na.rm=T) #stack daily to monthly
        writeRaster(frp,paste0(outfrpHome,"FRP",iYear*100+iMonth,"monthly"),format="GTiff",overwrite=T)
      }
    }
    
    if (mode=="seasonal") {
      
      if (next_yr!=FALSE) {
        #INPUT: Additional FRP
        frpHome_nextyr <- paste0(frpDir,iYear+1,"/")
        ifelse(!file.exists(frpHome_nextyr),stop("missing daily frp files/folder"),FALSE)
      }
      
      pointRasDir <- paste0(pointDir,"rasters/")
      
      ifelse(!file.exists(pointDir),dir.create(pointDir),FALSE)
      ifelse(!file.exists(pointRasDir),dir.create(pointRasDir),FALSE)
      
      #OUTPUT: monthly frp rasters
      outfrpDir <- paste0(pointRasDir,"frp_",season_name,nameAdd,"/")
      ifelse(!file.exists(outfrpDir),dir.create(outfrpDir),FALSE)
      
      st.jul = day_select[1] #start day
      end.jul = day_select[length(day_select)] #end day
      
      if ((iYear %% 4)==0) { #leap year compensate
        if (st.jul>28) {st.jul <- st.jul + 1}
        if (end.jul>=28) {end.jul <- end.jul + 1}
      }
      
      if (next_yr!=FALSE) { #next year days
        st.juln = next_yr[1]
        end.juln = next_yr[length(next_yr)]
        
        if ((iYear %% 4)==0) { #leap year
          if (st.juln>28) {st.juln <- st.juln + 1}
          if (end.juln>=28) {end.juln <- end.juln + 1}
        }
      }
      
      #STACK
      yrjul = (iYear*1000+st.jul):(iYear*1000+end.jul)
      frpnames <- paste0(frpHome,"FRP",yrjul,".tif")
      
      if (is.na(next_yr)==FALSE) {
        yrjuln = (iYear*1000+st.juln):(iYear*1000+end.juln)
        frpnames <- c(frpnames,paste0(frpHome_nextyr,"FRP",yrjuln,".tif"))
      }
      
      frp <- sum(stack(frpnames),na.rm=T) #stack daily to seasonal
      writeRaster(frp,paste0(outfrpDir,"FRP",iYear,season_name),format="GTiff",overwrite=T)
    }
    
    cat(iYear,": year",which(year_select==iYear),"of",length(year_select),as.character(Sys.time()),'\n')
  }
  print("Complete!")
}