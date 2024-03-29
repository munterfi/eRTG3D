% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/glidingSoaring3D.R
\name{n.sim.glidingSoaring.3d}
\alias{n.sim.glidingSoaring.3d}
\title{Simulates multiple 'gliding & soaring' tracks with a given number of gliding steps}
\usage{
n.sim.glidingSoaring.3d(
  n.sim = 1,
  parallel = FALSE,
  MODE,
  dGliding,
  dSoaring,
  qGliding,
  start = c(0, 0, 0),
  end = start,
  a0,
  g0,
  error = TRUE,
  smoothTransition = TRUE,
  glideRatio = 20,
  DEM = NULL,
  BG = NULL,
  verbose = FALSE
)
}
\arguments{
\item{n.sim}{number of simulations to produce}

\item{parallel}{logical: run computations in parallel (n-1 cores)? Or numeric: the number of nodes (maximum: n - 1 cores)}

\item{MODE}{raster layer containing the number/index of the mode, which should be used at each location}

\item{dGliding}{density object returned by the \link[eRTG3D]{get.densities.3d} function for gliding mode}

\item{dSoaring}{density object returned by the \link[eRTG3D]{get.densities.3d} function for soaring mode}

\item{qGliding}{the Q probabilites for the steps in gliding mode (\link[eRTG3D]{qProb.3d})}

\item{start}{numeric vector of length 3 with the coordinates of the start point}

\item{end}{numeric vector of length 3 with the coordinates of the end point}

\item{a0}{initial incoming heading in radian}

\item{g0}{initial incoming gradient/polar angle in radian}

\item{error}{logical: add random noise to the turn angle, lift angle and step length to account for errors measurements?}

\item{smoothTransition}{logical: should the transitions between soaring and the following gliding sections be smoothed? Recommended to avoid dead ends}

\item{glideRatio}{ratio between vertical and horizontal movement, by default set to 15 meters forward movement per meter vertical movement}

\item{DEM}{raster layer containing a digital elevation model, covering the area between start and end point}

\item{BG}{a background raster layer that can be used to inform the choice of steps}

\item{verbose}{logical: print current mode used?}
}
\value{
A list containing 'soaring-gliding' trajectories or \code{NULL}s if dead ends have been encountered.
}
\description{
Creates conditional empirical random walks in gliding mode, between a start and end point.
The walk is performed on a MODE layer and, if provided, additionally on a background and digital elevation layer.
The gliding is simulated with \link[eRTG3D]{sim.cond.3d} and soaring with \link[eRTG3D]{sim.uncond.3d},
therefore soaring is not restricted towards the target and can happen completly free as long as there are good thermal conditions.
It is important to extract for every mode in the MODE raster layer a corresponding densities object with \link[eRTG3D]{get.densities.3d}
and pass them to the function.
}
\note{
The MODE raster layer must be in the following structure: Gliding pixels have the value 1 and soaring pixel the values 2. \code{NA}'s are not allowed in the raster.
}
\examples{
print("tbd.")
}
