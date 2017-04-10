#' Create an Ocean Mask
#'
#' Create an ocean mask from the countries.tif raster derived from the
#' Natural Earth shapefile ne_10m_admin_0_countries using QGIS.
#' @usage oceanMask(...)
#' @param ... see global arguments: ask_home
#' @keywords ocean mask
#' @export

oceanMask <- function(ask_home = TRUE) {
  
  options(error=NULL)
  
  Home <- getwd()

  if (ask_home=="TRUE") {
    DirAns <- readline(paste("Is the project home directory:",Home,"? (enter y/n) "))
    if (DirAns!="y") {stop("Set project home directory with setwd()")}
  }
  
  rasDir <- "./rasters/"
  ifelse(!file.exists(rasDir),dir.create(rasDir),FALSE)
  
  #INPUT
  #Choose a FRP .tif file processed by MODIS Reprojection Tool
  frpname <- dir("raw_files/",all.files = T,pattern="*\\.tif$",recursive=T)[1]
  print(frpname)
  
  frpAns <- readline(paste("This is the FRP file used, okay? (enter y/n) "))
  if (frpAns=="n") {frpname <- readline(paste("Enter path of FRP file (format it like the example) that you want to use. "))}
  
  frp <- raster(paste0("raw_files/",frpname))
  
  #Load countries.tif file, make sure it is in the rasters folder
  land <- raster("rasters/countries.tif")
  
  #Reproject to parameters of FRP file
  ocean_mask <- suppressWarnings(projectRaster(land,frp,method="ngb"))
  ocean_mask[ocean_mask==0] = NA #set ocean as NA
  ocean_mask[ocean_mask>0] = 1 #set land as 1
  
  #Save ocean mask
  writeRaster(ocean_mask,paste0(rasDir,"ocean_mask.tif"),format="GTiff",overwrite=T)
  print("Complete!")
}
