% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dailyFRPstack.R
\name{dailyFRPstack}
\alias{dailyFRPstack}
\title{Sum and Stack MODIS FRP to Daily}
\usage{
... dailyFRPstack(year_select, sDay, eDay,
   frp_mode="combined", ocean_mask=TRUE, ...)
}
\arguments{
\item{year_select}{numeric vector. Years to be run}

\item{sDay}{numeric. First day of starting 8-day composite}

\item{eDay}{numeric. First day of ending 8-day composite}

\item{frp_mode}{character: "terra", "aqua", or "combined". Choose Terra-only, Aqua-only, or Terra+Aqua FRP. [default is "combined"]}

\item{ocean_mask}{logical. Mask ocean? [default is TRUE]}

\item{...}{see global arguments: rawDir, rasDir, ask_home}
}
\description{
Stack and sum Terra and Aqua daily FRP & mask ocean (for HYSPLIT). Input takes processed FRP .tif files the MODIS Reprojection Tool.
}
\examples{
#renames and stacks FRP in the years 2007-2013 for all days and masks ocean
dailyFRPstack(year_select=2007:2013,sDay=1,eDay=361)

#renames terra FRP in the years 2007-2013 for all days and masks ocean
dailyFRPstack(year_select=2007:2013,sDay=1,eDay=361,mode="terra")
}
\keyword{daily}
\keyword{frp}
\keyword{modis}
\keyword{stack}

