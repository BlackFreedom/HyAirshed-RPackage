#' Create Plots for HYSPLIT-derived Composite Seasonal Airshed
#' 
#' Create plot for composite seasonal airshed for a given point location
#' @usage plot.airshedC(loc, lon, lat, season_name, borderName="ne_10m_admin_0_countries",
#'    x.adj=FALSE, zoom=4.5, x.just=0, y.just=0, nKM=150, col.adj=c(0,0,0,0),
#'    scale.adj=c(0,0), caption=NA ,caption.adj=c(0,0),
#'    nYrs.adj=c(0,0), legend.adj=c(0,0), txt_size=0.75, ...)
#' @param loc character. Name of location, starting point, used in naming subdirectories and files
#' @param lon numeric. Longitude of starting point
#' @param lat numeric. Latitude of starting point
#' @param season_name character. Name of season, used in the names of the output folder and file
#' @param borderName character. Name of shapefile located in the "shapefiles" folder of project home directory to plot country/sub-country borders [default is "ne_10m_admin_0_countries"]
#' @param x.adj logical. Adjust x axis latitude to 0-360? [default is FALSE] (-180 to 180)
#' @param zoom numeric. Window of the plot, in degrees, from the lon, lat of the starting location in all 4 directions [default is 4.5]
#' @param x.just numeric. Adjust window laterally, positive is right and negative is left [default is 0]
#' @param y.just numeric. Adjust window vertically, positive is right and negative is left [default is 0]
#' @param nKM numeric. Length of scale bar, in kilometers [default is 150]
#' @param col.adj numeric vector, length=4. Adjust color bar position/size (format is c(x.min, x.max, y.min, y.max) adjustment) [default is c(0,0,0,0)]
#' @param scale.adj numeric vector, length=2. Adjust scale bar position (format is c(x, y) adjustment), [default is c(0,0)]
#' @param caption character or expression. Name of caption, [default is NA], which is no caption
#' @param caption.adj numeric vector, length=2. Adjust caption position (format is c(x, y) adjustment) [default is c(0,0)]
#' @param nYrs.adj numeric vector, length=2. Adjust number of years label position (format is c(x, y) adjustment) [default is c(0,0)]
#' @param legend.adj numeric vector, length=2. Adjust legend position (format is c(x, y) adjustment), [default is c(0,0)]
#' @param txt_size numeric. Adjust text size of plot labels, based on size of axis tick labels [default is 0.75]
#' @param ... see global arguments: pointDir, ask_home
#' @keywords hysplit airshed seasonal plot
#' @export
#' @examples 
#' #creates plot for composite winter delhi airshed for 2007-2013
#' plot.airshedC(loc="delhi",lon=77.2090,lat=28.6139,x.just=2,y.just=-2,caption="Delhi Winter Composite Airshed \nOct-Nov")

plot.airshedC <- function(loc, lon, lat, season_name, borderName="ne_10m_admin_0_countries", x.adj=FALSE,
                          zoom=4.5, x.just=0, y.just=0, nKM=150, col.adj=c(0,0,0,0),
                          scale.adj=c(0,0), caption=NA ,caption.adj=c(0,0), nYrs.adj=c(0,0),
                          legend.adj=c(0,0), txt_size=0.75, pointDir=paste0("./",loc,"/"),
                          ask_home=TRUE) {
  
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
  
  #INPUT: airsheds
  airHome <- paste0(pointDir,"/airshed_",season_name,"/")
  airHomeSub <- paste0(airHome,"shapefiles/")
  
  #OUTPUT: airshed plots
  outputHome <- paste0(airHome,"plots/")
  ifelse(!file.exists(outputHome),dir.create(outputHome,recursive=T),FALSE)
  
  border <- readShapePoly(paste0("./shapefiles/",borderName),proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  print("Okay, wait for function to finish.")
  
  mdeg <- function(reflat) {
    rlat = reflat*pi/180;
    mlat = 111412.84*cos(rlat)-93.5*cos(3*rlat)
    return(mlat)
  }
  
  x.labs <- seq(-180,180,2); x.labs_s = seq(-180,180,0.5)
  y.labs <- seq(-180,180,2); y.labs_s = seq(-180,180,0.5)
  
  x.labsN <- c()
  for(i in 1:length(x.labs)) {
    if (x.adj==FALSE) {
      if (x.labs[i]>0) {x.labsN[i] <- paste0(x.labs[i],"°E")}
      if (x.labs[i]<0) {x.labsN[i] <- paste0(abs(x.labs[i]),"°W")}
      if (x.labs[i]==0) {x.labsN[i] <- 0}
    }
    
    if (x.adj==TRUE) {
      if (x.labs[i]>180) {x.labsN[i] <- paste0(x.labs[i]-360,"°W")}
      if (x.labs[i]<180) {x.labsN[i] <- paste0(x.labs[i],"°E")}
      if (x.labs[i]==180) {x.labsN[i] <- 0}
    }
  }
  
  y.labsN <- c()
  for(i in 1:length(y.labs)) {
    if (y.labs[i]>0) {y.labsN[i] <- paste0(y.labs[i],"°N")}
    if (y.labs[i]<0) {y.labsN[i] <- paste0(abs(y.labs[i]),"°S")}
    if (y.labs[i]==0) {y.labsN[i] <- 0}
  }
  
  box <- c(lon-zoom-x.just,lon+zoom-x.just,lat-zoom-y.just,lat+zoom-y.just)
  city_coor <- SpatialPoints(as.matrix(cbind(lon,lat)),proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  
  quartz(w=5.8,h=4.8);par(mar=c(2.3,2,1,1.2))
  plot(border,xlim=c(box[1],box[2]),ylim=c(box[3],box[4]),lwd=0.7,border="gray50");box()
  axis(s=1,tck=0.02,labels=NA,at=x.labs);axis(s=1,tck=0.01,labels=NA,at=x.labs_s);
  axis(s=1,lwd=0,line=-0.9,cex.axis=txt_size,las=1,x.labsN,at=x.labs)
  axis(s=2,tck=0.02,labels=NA,at=y.labs);axis(s=2,tck=0.01,labels=NA,at=y.labs_s)
  axis(s=2,lwd=0,line=-0.7,cex.axis=txt_size,las=3,labels=y.labsN,at=y.labs)
  axis(s=3,tck=0.02,labels=NA,at=x.labs);axis(s=3,tck=0.01,labels=NA,at=x.labs_s);
  axis(s=4,tck=0.02,labels=NA,at=y.labs);axis(s=4,tck=0.01,labels=NA,at=y.labs_s)
  
  airshedC_Yr <- readShapePoly(paste0(airHomeSub,"/airshed_",season_name,"_compositeYr"),proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  airshedC <- readShapePoly(paste0(airHomeSub,"/airshed_",season_name,"_composite"),proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  polyC <- gUnaryUnion(airshedC)
  
  #colormap for pixels
  color <- c("transparent",rev(paste0("gray",seq(1,99,1))),1)
  pixcol <- c()
  for (iPix in 1:length(airshedC_Yr$den)) {
    pixcol[iPix] = color[round(airshedC_Yr$den[iPix]*100/max(airshedC_Yr$den,na.rm=T))]
  }
  
  plot(airshedC_Yr,add=T,col=pixcol,lwd=0.05)
  plot(polyC,add=T,border=1,lwd=2)
  
  plot(city_coor,col=2,add=T,cex=1,pch=4,lwd=2)
  x = 1:10; y = 1:10; z = outer(x,y,"+")
  obj = list(x=x,y=y,z=z)
  fac <- round(max(airshedC_Yr$den)/5)
  ipos <- c(0.134,0.164,0.14,0.5)
  image.plot(obj,add=T,legend.only=T,zlim=c(0,max(airshedC_Yr$den)),col=color[-1],axes=F,
             smallplot=ipos+col.adj,axis.args=list(at=seq(0,max(airshedC_Yr$den),fac),labels=NA,tck=-0.5))
  image.plot(obj,add=T,legend.only=T,zlim=c(0,max(airshedC_Yr$den)),col=color[-1],axes=F,
             smallplot=ipos+col.adj,axis.args=list(at=seq(0,max(airshedC_Yr$den),fac),line=-0.4,labels=seq(0,max(airshedC_Yr$den),fac),cex.axis=0.7,lwd=0),
             legend.args=list(text="Years (# of Total)",side=2,font=2,line=0.3,cex=0.75,adj=0.5))
  
  par(new=T);plot(NA,NA,xlim=c(0,10),ylim=c(0,10),axes=F,xlab="",ylab="")
  
  spos <- c(1.6,0.2)+scale.adj
  
  segments(spos[1],spos[2],spos[1]+nKM*1000/mdeg(spos[2]),spos[2],lwd=1.2)
  segments(spos[1],spos[2]-0.12,spos[1],spos[2]+0.12,lwd=1.2)
  segments(spos[1]+nKM*1000/mdeg(spos[2]),spos[2]-0.12,spos[1]+nKM*1000/mdeg(spos[2]),spos[2]+0.12,lwd=1.2)
  text(mean(c(spos[1],spos[1]+nKM*1000/mdeg(spos[2]))),spos[2]+0.35,paste(nKM,"km"),cex=txt_size)
  
  apos <- c(10.45,-1.2)+nYrs.adj
  text(apos[1],apos[2],paste0("Number of Years = ",max(airshedC_Yr$den)),xpd=NA,cex=0.6,font=3,adj=1)
  
  cpos <- c(0.05,9.9)+caption.adj
  if (is.na(caption)==FALSE) {
    text(cpos[1],cpos[2],caption,xpd=NA,cex=txt_size,font=2,adj=c(0,1))
  }
  
  quartz.save(paste0(outputHome,"airshed_",season_name,"_comp.pdf"),type="pdf")
  graphics.off()
  
  print("Complete!")
}