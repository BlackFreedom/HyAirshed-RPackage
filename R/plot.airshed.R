#' Create Plots for HYSPLIT-derived Seasonal Airshed
#' 
#' Create plots for seasonal airsheds with daily HYSPLIT back trajectories for a given point location
#' @usage plot.airshed(loc, lon, lat, day_select, year_select,
#'    tYears=year_select[1]:year_select[length(year_select)],
#'    season_name, hy_alt, alt.adj=0, sPix=0.25, nPix=100,
#'    borderName="ne_10m_admin_0_countries", x.adj=FALSE, zoom=4.5, x.just=0,
#'    y.just=0, nKM=150, traj.col=c("#a6cee3","#fb9a99","#cab2d6"),
#'    traj.lwd=0.25, col.adj=c(0,0,0,0), scale.adj=c(0,0),
#'    bar.adj=c(0,0,0,0), caption=NA, caption.adj=c(0,0),
#'    nAlt.adj=c(0,0), legend.lab=hy_alt,
#'    legend.adj=c(0,0), txt_size=0.75, ...)
#' @param loc character. Name of location, starting point, used in naming subdirectories and files
#' @param lon numeric. Longitude of starting point
#' @param lat numeric. Latitude of starting point
#' @param month_select numeric vector. Months to be run
#' @param year_select numeric vector. Years to be run
#' @param tYears numeric vector. Total range of years in study, ignoring missing years [default = year_select[1]:year_select[length(year_select)], which is the years from the first year to last year in year_select
#' @param season_name character. Name of season, used in the names of the output folder [default is “season"]
#' @param hy_alt character string or vector. Name(s) of hysplit heights to be run
#' @param alt.adj numeric. Adjustment if more than one location [default is 0]
#' @param sPix numeric. Size of pixel, in degrees [default is 0.3]
#' @param nPix numeric. Number of pixels in horizontal and vertical direction [default is 100]
#' @param borderName character. Name of shapefile located in the "shapefiles" folder of project home directory to plot country/sub-country borders [default is "ne_10m_admin_0_countries"]
#' @param x.adj logical. Adjust x axis latitude to 0-360? default is FALSE (-180 to 180)
#' @param zoom numeric. Window of the plot, in degrees, from the lon, lat of the starting location in all 4 directions [default is 4.5]
#' @param x.just numeric. Adjust window laterally, positive is right and negative is left [default is 0]
#' @param y.just numeric. Adjust window vertically, positive is right and negative is left [default is 0]
#' @param traj.col character string or vector in html or R color name. Colors of trajectories [default is c("#a6cee3","#fb9a99","#cab2d6")]
#' @param traj.lwd numeric. Width of trajectory lines [default is 0.3]
#' @param nKM numeric. Length of scale bar, in kilometers [default is 150]
#' @param col.adj numeric vector, length=4. Adjust color bar position/size (format is c(x.min, x.max, y.min, y.max) adjustment) [default is c(0,0,0,0)]
#' @param scale.adj numeric vector, length=2. Adjust scale bar position (format is c(x, y) adjustment), [default is c(0,0)]
#' @param bar.adj numeric vector, length=4. Adjust timestamp bar position/size (format is c(x.min, x.max, y.min, y.max) adjustment), [default is c(0,0,0,0)]
#' @param caption character or expression. Name of caption, [default is NA], which is no caption
#' @param caption.adj numeric vector, length=2. Adjust caption position (format is c(x, y) adjustment) [default is c(0,0)]
#' @param nAlt.adj numeric vector, length=2. Adjust number of trajectories label position (format is c(x, y) adjustment) [default is c(0,0)]
#' @param legend.lab character vector. Labels for back trajectories [default is hy_alt]
#' @param legend.adj numeric vector, length=2. Adjust legend position (format is c(x, y) adjustment), [default is c(0,0)]
#' @param txt_size numeric. Adjust text size of plot labels, based on size of axis tick labels [default is 0.75]
#' @param ... see global arguments: pointDir, ask_home
#' @keywords hysplit airshed seasonal plot
#' @export
#' @examples 
#' #creates plots for winter delhi airshed for 2007-2013 and composite winter airshed for these years
#' plot.airshed(loc="delhi",lon=77.2090,lat=28.6139,day_select=274:334,year_select=2007:2013,season_name="win",hy_alt="0.5km",x.just=1.5,y.just=-1.5,caption="Delhi Winter Airshed \n(OCT-NOV)",legend.lab="0.5 km")

plot.airshed <- function(loc, lon, lat, day_select, year_select,
                         tYears=year_select[1]:year_select[length(year_select)],
                         season_name, hy_alt, alt.adj=0, sPix=0.25, nPix=100,
                         borderName="ne_10m_admin_0_countries", x.adj=FALSE, zoom=4.5, x.just=0,
                         y.just=0, nKM=150, traj.col=c("#a6cee3","#fb9a99","#cab2d6"),
                         traj.lwd=0.25, col.adj=c(0,0,0,0), scale.adj=c(0,0),
                         bar.adj=c(0,0,0,0), caption=NA, caption.adj=c(0,0),
                         nAlt.adj=c(0,0), legend.lab=hy_alt,
                         legend.adj=c(0,0), txt_size=0.75,
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
  
  #INPUT: hysplit gis & txt files
  #subfolders = iYear (e.g. 2007) -> month (e.g. oct)
  trajHome <- paste0(pointDir,"hysplit/")
  ifelse(!file.exists(trajHome),stop("missing hysplit coordinates folder"),FALSE)
  
  #INPUT: hysplit pdf files
  #subfolders = iYear (e.g. 2007) -> month (e.g. oct)
  hypdfHome <- "./hysplit/hysplit_plot/"
  ifelse(!file.exists(trajHome),stop("missing hysplit pdf plots directory"),FALSE)
  
  #INPUT: airsheds
  airHome <- paste0("./",loc,"/airshed_",season_name,"/")
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
  
  gt <- GridTopology(c(lon-(sPix*nPix/2),lat-(sPix*nPix/2)),c(sPix,sPix),c(nPix,nPix))
  grd <- SpatialGrid(gt,proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  spix <- as(grd,"SpatialPixels")
  spol <- as(spix,"SpatialPolygons")
  spol$coor <- coordinates(spol)
  
  for (iYear in year_select) {
    
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
    
    quartz(w=5.8,h=5);par(mar=c(2.2,2,2.3,1.2))
    plot(border,xlim=c(box[1],box[2]),ylim=c(box[3],box[4]),lwd=0.7,border="gray50");box()
    axis(s=1,tck=0.02,labels=NA,at=x.labs);axis(s=1,tck=0.01,labels=NA,at=x.labs_s);
    axis(s=1,lwd=0,line=-0.9,cex.axis=txt_size,las=1,x.labsN,at=x.labs)
    axis(s=2,tck=0.02,labels=NA,at=y.labs);axis(s=2,tck=0.01,labels=NA,at=y.labs_s)
    axis(s=2,lwd=0,line=-0.7,cex.axis=txt_size,las=3,labels=y.labsN,at=y.labs)
    axis(s=3,tck=0.02,labels=NA,at=x.labs);axis(s=3,tck=0.01,labels=NA,at=x.labs_s);
    axis(s=4,tck=0.02,labels=NA,at=y.labs);axis(s=4,tck=0.01,labels=NA,at=y.labs_s)
    
    city_coor <- SpatialPoints(as.matrix(cbind(lon,lat)),proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
    
    arishedS <- readShapePoly(paste0(airHomeSub,"airshed",iYear,season_name),proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
    polyS <- gUnaryUnion(arishedS)
    
    #colormap for pixels
    color <- c("transparent",rev(paste0("gray",seq(1,99,1))),1)
    pixcol <- c()
    for (iPix in 1:length(arishedS)) {
      pixcol[iPix] = color[round(arishedS$den[iPix]*100/max(arishedS$den,na.rm=T))]
    }
    
    plot(polyS,add=T,border=1,lwd=2)
    plot(arishedS,add=T,col=pixcol,lwd=0.05)
    
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
    
    for (iAltitude in 1:length(hy_alt)) {
      iAltitudeName <- hy_alt[iAltitude]
      shape.lat <- read.table(paste0(trajHome,loc,"_HysplitLat_",iYear,"_",iAltitudeName,".csv"),sep=",")
      shape.lon <- read.table(paste0(trajHome,loc,"_HysplitLon_",iYear,"_",iAltitudeName,".csv"),sep=",")
      
      nTraj <- shape.lat[,1][day_select.adj]
      nTraj <- length(nTraj[!is.na(nTraj)])
      
      for(iDay in day_select.adj) {
        if (length(shape.lat[iDay,][!is.na(shape.lat[iDay,])])>0) {
          shape.coor <- cbind(as.numeric(shape.lat[iDay,]),as.numeric(shape.lon[iDay,]))
          shape.coor <- apply(shape.coor,2,na.omit)
          shape <- SpatialPoints(shape.coor,proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
          traj <- as(shape,"SpatialLines")
          
          plot(shape,cex=0.1,lwd=traj.lwd,add=T,col=traj.col[iAltitude])
          plot(traj,lwd=0.3,add=T,col=traj.col[iAltitude])
        }
      }
    }
    
    plot(city_coor,col=2,add=T,cex=1,pch=4,lwd=2)
    x <- 1:10; y <- 1:10; z <- outer(x,y,"+")
    obj <- list(x=x,y=y,z=z)
    ipos <- c(0.134,0.164,0.12,0.42)
    
    image.plot(obj,add=T,legend.only=T,zlim=c(0,100),col=color[-1],axes=F,
               smallplot=ipos+col.adj,axis.args=list(at=seq(0,100,20),labels=NA,tck=-0.5))
    image.plot(obj,add=T,legend.only=T,zlim=c(0,100),col=color[-1],axes=F,
               smallplot=ipos+col.adj,axis.args=list(at=seq(0,100,20),line=-0.4,labels=seq(0,100,20),cex.axis=0.7,lwd=0),
               legend.args=list(text="Trajectories (% of Total)",side=2,font=2,line=0.3,cex=0.75,adj=0.5))
    
    par(new=T);plot(NA,NA,xlim=c(0,10),ylim=c(0,10),axes=F,xlab="",ylab="")
    
    spos <- c(1.5,0.1)+scale.adj
    
    segments(spos[1],spos[2],spos[1]+nKM*1000/mdeg(spos[2]),spos[2],lwd=1.2)
    segments(spos[1],spos[2]-0.12,spos[1],spos[2]+0.12,lwd=1.2)
    segments(spos[1]+nKM*1000/mdeg(spos[2]),spos[2]-0.12,spos[1]+nKM*1000/mdeg(spos[2]),spos[2]+0.12,lwd=1.2)
    text(mean(c(spos[1],spos[1]+nKM*1000/mdeg(spos[2]))),spos[2]+0.35,paste(nKM,"km"),cex=txt_size)
    
    apos <- c(10.45,-1.27)+nAlt.adj
    text(apos[1],apos[2],paste0("Number of Trajectories Per Altitude = ",sum(nTraj)),xpd=NA,cex=0.6,font=3,adj=1)
    
    bar <- c(-0.4,10.4,10.55,11.35)+bar.adj
    nYears <- length(tYears)
    rect(bar[1],bar[3],bar[2],bar[4],xpd=T)
    rect(bar[1]+(bar[2]-bar[1])/length(tYears)*(iYear-min(tYears)),bar[3],bar[1]+(bar[2]-bar[1])/length(tYears)*(iYear-(min(tYears)-1)),bar[4],xpd=NA,col=1)
    text(bar[1]+(bar[2]-bar[1])/length(tYears)*(iYear-(min(tYears)-0.5)),mean(bar[3],bar[4])+0.5*(bar[4]-bar[3]),iYear,col="white",xpd=NA,cex=txt_size+0.05,adj=0.5)
    
    cpos <- c(0.05,9.9)+caption.adj
    if (is.na(caption)==FALSE) {
      text(cpos[1],cpos[2],caption,xpd=NA,cex=txt_size,font=2,adj=c(0,1))
    }
    
    lpos <- c(8.15,9.3)+legend.adj
    legend(lpos[1],lpos[2]-(length(hy_alt)-1)*0.2,legend=paste(legend.lab,""),lty=1,col=traj.col,cex=txt_size-0.05,title=expression(paste(bold("Altitudes"))),title.adj=0.25,yjust=0.5,bg="white")
    
    quartz.save(paste0(outputHome,"airshed",iYear,season_name,".pdf"),type="pdf")
    graphics.off()
    
    cat(iYear,": year",iYear-(year_select[1]-1),"of",length(year_select),as.character(Sys.time()),'\n')
  }
  
  print("Complete!")
}