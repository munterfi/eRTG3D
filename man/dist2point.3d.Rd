% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrapper3D.R
\name{dist2point.3d}
\alias{dist2point.3d}
\title{Distance of each track point to a given point}
\usage{
dist2point.3d(track, point, groundDistance = FALSE)
}
\arguments{
\item{track}{a list containing data.frames with x,y,z coordinates or a data.frame}

\item{point}{a vector with x, y or x, y, z coordinates}

\item{groundDistance}{logical: calculate only ground distance in x-y plane?}
}
\value{
Returns the distance of each track point to the point.
}
\description{
Distance of each track point to a given point
}
\examples{
dist2point.3d(niclas, c(0, 0, 0))
}
