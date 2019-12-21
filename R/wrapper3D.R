#' Extract tldCube and autodifferences functions from a consistent track
#'
#' Get densities creates a list consisting of the 3 dimensional
#' probability distribution cube for turning angle, lift angle and step length (\link[eRTG3D]{turnLiftStepHist})
#' as well as the uni-dimensional distributions of the differences
#' of the turning angles, lift angles and step lengths with a lag of 1 to maintain
#' minimal level of autocorrelation in each of the terms.
#'
#' @section Note:
#' The time between the acquisition of fix points  of the track must be constant,
#' otherwise this leads to distorted statistic distributions,
#' which increases the probability of dead ends. In this case please check
#' \link[eRTG3D]{track.split.3d} and \link[eRTG3D]{get.section.densities.3d}
#'
#' @param track a data.frame with 3 columns containing the x,y,z coordinates
#' @param gradientDensity logical: Should a distribution of the gradient angle be extracted and later used in the simulations?
#' @param heightDistEllipsoid logical: Should a distribution of the flight height over ellipsoid be extracted and later used in the sim.cond.3d()?
#' @param DEM a raster containing a digital elevation model, covering the same extent as the track
#' @param maxBin numeric scalar, maximum number of bins per dimension of the tld-cube (\link[eRTG3D]{turnLiftStepHist})
#'
#' @return A list containing the tldCube and the autodifferences functions (and additionally the height distribution function)
#' @export
#'
#' @examples
#' get.track.densities.3d(niclas, heightDist = TRUE)
get.track.densities.3d <- function(track, gradientDensity = TRUE, heightDistEllipsoid = TRUE, DEM = NULL, maxBin = 25)
{
  .is.df.xyz(track)
  track <- track.properties.3d(track)
  turnAngle <- track$t[2:nrow(track)]; liftAngle <- track$l[2:nrow(track)]; stepLength <- track$d[2:nrow(track)]
  deltaTurn <- diff(turnAngle); deltaLift <- diff(liftAngle); deltaStep <- diff(stepLength)
  if (gradientDensity) {gradientAngle <- track$g} else {gradientAngle <- NULL}
  if (heightDistEllipsoid) {heightEllipsoid <- track$z} else {heightEllipsoid <- NULL}
  if (!is.null(DEM)) {
    .check.extent(DEM = DEM, track = track)
    heightTopo <- track$z - raster::extract(DEM, track[,1:2])
  } else {heightTopo <- NULL}
  return(get.densities.3d(turnAngle = turnAngle, liftAngle = liftAngle, stepLength = stepLength,
                          deltaLift = deltaLift, deltaTurn = deltaTurn, deltaStep = deltaStep, gradientAngle = gradientAngle,
                          heightEllipsoid = heightEllipsoid, heightTopo = heightTopo, maxBin = maxBin))
}

#' Extract tldCube and autodifferences functions from track sections
#'
#' Creates a list consisting of the 3 dimensional
#' probability distribution cube for turning angle, lift angle and step length (\link[eRTG3D]{turnLiftStepHist})
#' as well as the uni-dimensional distributions of the differences
#' of the turning angles, lift angles and step lengths with a lag of 1 to maintain
#' minimal level of autocorrelation in each of the terms.
#'
#' @param trackSections list of track sections got by the \link[eRTG3D]{track.split.3d} function
#' @param gradientDensity logical: Should a distribution of the gradient angle be extracted and later used in the simulations?
#' @param heightDistEllipsoid logical: Should a distribution of the flight height over ellipsoid be extracted and later used in the sim.cond.3d()?
#' @param DEM a raster containing a digital elevation model, covering the same extent as the track sections
#' @param maxBin numeric scalar, maximum number of bins per dimension of the tld-cube (\link[eRTG3D]{turnLiftStepHist})
#'
#' @return A list containing the tldCube and the autodifferences functions (and additionally the height distribution function)
#' @export
#'
#' @examples
#' get.section.densities.3d(list(niclas[1:10, ], niclas[11:nrow(niclas), ]))
get.section.densities.3d <- function(trackSections, gradientDensity = TRUE, heightDistEllipsoid = TRUE, DEM = NULL, maxBin = 25)
{
  trackSections <- lapply(X=trackSections, FUN= function(X) track.properties.3d(X)[2:nrow(X), ])
  deltaTurn <- Reduce(c, lapply(X = trackSections, FUN = function(X) diff(X$t)))
  deltaLift <- Reduce(c, lapply(X = trackSections, FUN = function(X) diff(X$l)))
  deltaStep <- Reduce(c, lapply(X = trackSections, FUN = function(X) diff(X$d)))
  trackSections <- do.call(rbind, trackSections)
  turnAngle <- trackSections$t; liftAngle <- trackSections$l; stepLength <- trackSections$d
  if (gradientDensity) {gradientAngle <- trackSections$g} else {gradientAngle <- NULL}
  if (heightDistEllipsoid) {heightEllipsoid <- trackSections$z} else {heightEllipsoid <- NULL}
  if (!is.null(DEM)) {
    .check.extent(DEM = DEM, track = trackSections)
    heightTopo <- trackSections$z - raster::extract(DEM, trackSections[,1:2])
  } else {heightTopo <- NULL}
  return(get.densities.3d(turnAngle = turnAngle, liftAngle = liftAngle, stepLength = stepLength,
                          deltaLift = deltaLift, deltaTurn = deltaTurn, deltaStep = deltaStep, gradientAngle = gradientAngle,
                          heightEllipsoid = heightEllipsoid, heightTopo = heightTopo, maxBin = maxBin))
}

#' This function splits the by outliers in the time lag.
#'
#' The length of timeLag must be the the track's length minus 1 and represents
#' the time passed between the fix point acquisition
#'
#' @param track track data.frame with x, y and z coordinates
#' @param timeLag a numeric vector with the time passed between the fix point acquisition
#' @param lag NULL or a manually chosen lag
#' @param tolerance NULL or a manually chosen tolerance
#'
#' @return A list containing the splitted tracks.
#' @export
#'
#' @examples
#' track.split.3d(niclas, timeLag=rep(1, nrow(niclas)-1) + rnorm(nrow(niclas)-1, mean = 0, sd = 0.25))
track.split.3d <- function(track, timeLag, lag = NULL, tolerance = NULL)
{
  .is.df.xyz(track)
  if (any(is.na(timeLag))) stop("TimeLag is not allowed to contain NAs.")
  if (is.null(lag)) {
    m <- mean(timeLag)
  } else {m <- lag}
  if (is.null(tolerance)) {
    tolerance <- 0.5 * stats::sd(timeLag)
  }
  splitRows <- which(abs(m-timeLag) > tolerance)
  trackSections <- split(track, cumsum(1:nrow(track) %in% (splitRows+2))) # + 1 if the cut should be one step before
  nSplits <- length(splitRows); nChange <- round(sum(timeLag[splitRows]/m-1));
  message(paste("  |Mean time lag: ", round(m,2) , ", tolerance: ", round(tolerance,2),
                ", number of splits: ", nSplits, ", proposed change in steps: ", nChange, sep=""))
  return(trackSections)
}

#' Crops the DEM to the extent of the track with a buffer
#'
#' @param DEM a raster containing a digital elevation model, covering the extent as the track
#' @param track data.frame with x,y,z coordinates of the original track
#' @param buffer buffer with, by default set to 100
#'
#' @return A the cropped digital elevation model as a raster layer.
#' @export
#'
#' @examples
#' dem2track.extent(dem, niclas)
dem2track.extent <- function(DEM, track, buffer=100)
{
  .is.df.xyz(track = track)
  .check.extent(DEM = DEM, track = track)
  return(raster::crop(DEM, raster::extent(min(track$x)-buffer, max(track$x)+buffer, min(track$y)-buffer, max(track$y)+buffer)))
}

#' Extent of track(s)
#'
#' @param track a list containing data.frames with x,y,z coordinates or a data.frame
#' @param zAxis logical: return also the extent of the Z axis?
#'
#' @return
#' Returns an extent object of the raster package in the 2–D case and a vector in the 3–D case.
#' @export
#'
#' @examples
#' track.extent(niclas, zAxis = TRUE)
track.extent <- function(track, zAxis = FALSE){
  if (!is.list(track) || !is.list(track)) stop("Track input has to be of type list or data.frame.")
  if (is.list(track) && is.data.frame(track)) {track <- list(track)}
  extents <- do.call("rbind", lapply(X = track, FUN = function(track){
    c(floor(min(track$x)), floor(max(track$x))+1,
      floor(min(track$y)), floor(max(track$y))+1,
      floor(min(track$z)), floor(max(track$z))+1)
  }))
  minX <- min(extents[,1]); maxX <- max(extents[,2]);
  minY <- min(extents[,3]); maxY <- max(extents[,4]);
  minZ <- min(extents[,5]); maxZ <- max(extents[,6]);
  if (zAxis) {return(rbind(xmin=minX, xmax=maxX, ymin=minY, ymax=maxY, zmin=minZ, zmax=maxZ))}
  return(raster::extent(minX, maxX, minY, maxY))
}

#' Tests if the object is of type 'data.frame' and has x, y, z coordinates without NA values
#'
#' @param track any object to test
#'
#' @return A logical: TRUE if the track is the object needed, FALSE otherwise.
#' @export
#'
#' @examples
#' .is.df.xyz(niclas)
#' @noRd
.is.df.xyz <- function(track)
{
  if(!class(track)=="data.frame") stop("Input not of type 'data.frame'.")
  if(!any(colnames(track)[1:3]==c("x","y","z"))) stop("Colnames of first three cols not 'x, y, z'.")
  if(any(is.na(track[, 1:3]))) stop("Track 'data.frame' contains NA values.")
}

#' Checks if the track lies inside the digital elevation model.
#'
#' @param DEM a 'RasterLayer' containing a digital elevation model
#' @param track a data.frame with 3 columns containing the x,y,z coordinates
#'
#' @return A logical: TRUE if the track lies inside the DEM, FALSE otherwise.
#' @export
#'
#' @examples
#' .check.extent(dem, niclas)
#' @noRd
.check.extent <- function(DEM, track)
{
  if(!class(DEM)=="RasterLayer") stop("'DEM' is not of type 'RasterLayer'")
  e <- raster::extent(DEM)
  if(!(min(track[,1]) >= e[1] && max(track[,1]) <= e[2] &&
       min(track[,2]) >= e[3] && max(track[,2]) <= e[4])) stop("The track is not inside the area of the digital elevation model.")
}

#' Reproduce a track with the eRTG3D
#'
#' Simulates n tracks with the geometrical properties of the original track,
#' between the same start and end point.
#'
#' @param track data.frame with x,y,z coordinates of the original track
#' @param n.sim number of simulations that should be done
#' @param error logical: add error term to movement in simulation?
#' @param DEM a raster containing a digital elevation model, covering the same extent as the track
#' @param BG a raster influencing the probabilities.
#' @param parallel logical: run computations in parallel (n-1 cores)? Or numeric: the number of nodes (maximum: n - 1 cores)
#' @param plot2d logical: plot tracks on 2-D plane?
#' @param plot3d logical: plot tracks in 3-D?
#' @param maxBin numeric scalar, maximum number of bins per dimension of the tld-cube (\link[eRTG3D]{turnLiftStepHist})
#' @param gradientDensity logical: Should a distribution of the gradient angle be extracted and used in the simulations (\link[eRTG3D]{get.densities.3d})?
#' @param filterDeadEnds logical: Remove tracks that ended in a dead end?
#'
#' @return A list or data.frame containing the simulated track(s) (CERW).
#' @export
#'
#' @examples
#' reproduce.track.3d(niclas[1:10, ])
reproduce.track.3d <- function(track, n.sim = 1, parallel = FALSE, error = TRUE, DEM = NULL, BG = NULL, filterDeadEnds = TRUE, plot2d = FALSE, plot3d = FALSE, maxBin = 25, gradientDensity = TRUE)
{
  .is.df.xyz(track = track)
  track <- track.properties.3d(track)
  n.locs <- nrow(track)
  if (n.locs > 1500) message("  |Note: The track is very long (>1500 steps), the system may run out of memory.")
  turnAngle <- track$t[2:nrow(track)]; liftAngle <- track$l[2:nrow(track)]; stepLength <- track$d[2:nrow(track)]
  deltaTurn <- diff(turnAngle); deltaLift <- diff(liftAngle); deltaStep <- diff(stepLength)
  heightEllipsoid <- track$z
  if (gradientDensity) {gradientAngle <- track$g} else {gradientAngle <- NULL}
  if (!is.null(DEM)) {
    .check.extent(DEM = DEM, track = track)
    heightTopo <- track$z - raster::extract(DEM, track[,1:2])
  } else {heightTopo <- NULL}
  D <- get.densities.3d(liftAngle = liftAngle, turnAngle = turnAngle, stepLength = stepLength,
                        deltaLift = deltaLift, deltaTurn = deltaTurn, deltaStep = deltaStep, gradientAngle = gradientAngle,
                        heightEllipsoid = heightEllipsoid, heightTopo = heightTopo, maxBin = maxBin)
  uerw <- sim.uncond.3d(n.locs*1500, start = c(track$x[1],track$y[1],track$z[1]),
                        a0 = track$a[1], g0 = track$g[1], densities = D, error = error)
  Q <- qProb.3d(uerw, n.locs, parallel = parallel, maxBin = maxBin)
  cerwList <- suppressWarnings(n.sim.cond.3d(n.sim = n.sim, n.locs <- n.locs, start=c(track$x[1],track$y[1],track$z[1]), end=c(track$x[n.locs],track$y[n.locs],track$z[n.locs]),
                                             a0 = track$a[1], g0 = track$g[1], densities=D, qProbs=Q, error = error, parallel = parallel, DEM = DEM, BG = BG))
  if(filterDeadEnds){cerwList <- filter.dead.ends(cerwList)}
  if(plot2d){print(plot2d(origTrack = track, simTrack = cerwList, DEM = DEM))}
  if(plot3d){plot3d(origTrack = track, simTrack = cerwList, DEM = DEM)}
  return(cerwList)
}

#' Remove dead ends
#'
#' Function to filter out tracks that have found a dead end
#'
#' @param cerwList list of data.frames and NULL entries
#'
#' @return A list that is only containing valid tracks.
#' @export
#'
#' @examples
#' filter.dead.ends(list(niclas, niclas))
filter.dead.ends <- function(cerwList)
{
  if(is.null(cerwList)) {warning("No track made it to the end point."); return(NULL)}
  l1 <- length(cerwList)
  cerwList <- cerwList[!unlist(lapply(cerwList, is.null))]
  l2 <- length(cerwList)
  if (l2 == 0) {warning("No track made it to the end point."); return(NULL)}
  if (l1!=l2) {message(paste("  |Dead end tracks removed (n = ", (l1-l2), ", proportion: ", 1-round(l2/l1, 2), ")", sep=""))}
  return(cerwList)
}

#' Moving median in one dimension
#'
#' Applies a twosided moving median window on a vector,
#' where the window paramter is the total size of the window.
#' The value in the window middle is the index where the median of the window is written.
#' Therefore the window size has to be an uneven number.
#' The border region of the vetor is filled with a one-sided median.
#' There might be border effects.
#'
#' @param data numeric vector
#' @param window uneven number for the size of the moving window
#'
#' @return A numeric vector.
#' @export
#'
#' @examples
#' movingMedian(sequence(1:10), window = 5)
movingMedian <- function(data, window){
  if(!(window %% 2 == 0)) {window <- floor(window/2)} else{stop("Window must be an uneven number.")}
  total <- length(data)
  result <- vector(length = total)
  for(i in (window+1):(total-window)){
    result[i] <- stats::median(data[(i-window):(i+window)])
  }
  result[1:window] <- stats::median(data[1:window])
  result[(total-window):total] <- stats::median(data[(total-window):total])
  return(result)
}

#' Turn angle to target
#'
#' Calculates the turn angle between every point in the track and the last point (target).
#'
#' @param track a track data.frame containing x, y and z coordinates
#'
#' @return A numeric vector with the turn angles to target
#' @export
#'
#' @examples
#' turn2target.3d(niclas)
turn2target.3d <- function(track) {
  .is.df.xyz(track = track)
  track <- track.properties.3d(track)
  target <- Reduce(c, track[nrow(track), 1:3])
  .wrap(atan2(target[2]-track$y, target[1] - track$x) - track$a)}

#' Lift angle to target
#'
#' Calculates the lift angle between every point in the track and the last point (target).
#'
#' @param track a track data.frame containing x, y and z coordinates
#'
#' @return A numeric vector with the lift angles to target
#' @export
#'
#' @examples
#' lift2target.3d(niclas)
lift2target.3d <- function(track) {
  .is.df.xyz(track = track)
  track <- track.properties.3d(track)
  target <- Reduce(c, track[nrow(track), 1:3])
  .wrap(atan2(sqrt((target[1]-track$x) ^ 2 + (target[2]-track$y) ^ 2),
              (target[3]-track$z)) - track$g)}

#' Distance to target
#'
#' Calculates the distance between every point in the track and the last point (target).
#'
#' @param track a track data.frame containing x, y and z coordinates
#'
#' @return A numeric vector with the distances to target
#' @export
#'
#' @examples
#' dist2target.3d(niclas)
dist2target.3d <- function(track) {
  .is.df.xyz(track = track)
  target <- Reduce(c, track[nrow(track), 1:3])
  sqrt((target[1]-track$x) ^ 2 + (target[2]-track$y) ^ 2 + (target[3]-track$z) ^ 2)}

#' Distance of each track point to a given point
#'
#' @param track a list containing data.frames with x,y,z coordinates or a data.frame
#' @param point a vector with x, y or x, y, z coordinates
#' @param groundDistance logical: calculate only ground distance in x-y plane?
#'
#' @return Returns the distance of each track point to the point.
#' @export
#'
#' @examples
#' dist2point.3d(niclas, c(0,0,0))
dist2point.3d <- function(track, point, groundDistance = FALSE) {
  .is.df.xyz(track = track)
  if (groundDistance | length(point)==2) {
    sqrt((point[1]-track[,1])^2 + (point[2]-track[,2])^2)
  } else {
    sqrt((point[1]-track[,1])^2 + (point[2]-track[,2])^2 + (point[3]-track[,3])^2)
  }
}

#' Track properties of a 3-D track
#'
#' Returns the properties (distances, azimuth, polar angle,
#' turn angle & lift angle) of a track in three dimensions.
#'
#' @param track data.frame with x,y,z coordinates
#'
#' @return The data.frame with track properties
#' @export
#'
#' @examples
#' track.properties.3d(niclas)
track.properties.3d <- function(track)
{
  .is.df.xyz(track)
  # spatial coordinates
  x <- track[,1]; y <- track[,2]; z <- track[,3]
  # distance covered per axis and step
  dx <- c(NA, diff(x)); dy <- c(NA, diff(y)); dz <- c(NA, diff(z))
  d2d <- .distance.2d(dx, dy)
  # spherical coordinates
  d <- .distance.3d(dx, dy, dz)
  a <- .get.azimut(dx, dy)
  g <- .get.polar(d2d, dz)
  # guess the initial azimut/heading & polar angle/gradient (alt: sample(a,1))
  a[1] <- mean(a, na.rm=TRUE)
  g[1] <- mean(g, na.rm=TRUE)
  # Turn angle and lift angle
  t <- c(NA, diff(a))
  l <- c(NA, diff(g))
  return(data.frame(x, y, z, a, g, t, l, d))
}

#' Azimut in respect to the x-axis
#'
#' @param dx distance in x
#' @param dy distance in y
#'
#' @return The azimut in radian
#' @export
#'
#' @examples
#' .get.azimut(10, 30)
#' @noRd
.get.azimut <- function(dx, dy)
{
  .wrap(atan2(dy, dx))
}

#' Polar angle/gradient in respect to the z-axis
#'
#' @param d ground distance in xy plane
#' @param dz distance in z
#'
#' @return The polar angle in radian
#' @export
#'
#' @examples
#' .get.polar(10, 50)
#' @noRd
.get.polar <- function(d, dz)
{
  .wrap(atan2(d, dz))
}

#' Ground distance covered in XY plane
#'
#' @param dx distance in x
#' @param dy distance in y
#'
#' @return Ground distance in XY plane
#' @export
#'
#' @examples
#' .distance.2d(10, 30)
#' @noRd
.distance.2d <- function(dx, dy)
{
  sqrt(dx*dx+dy*dy)
}

#' Distance covered in 3 dimensions, radius
#'
#' @param dx distance in x
#' @param dy distamce in y
#' @param dz distance in z
#'
#' @return The radius
#' @export
#'
#' @examples
#' .distance.3d(10, 10, 10)
#' @noRd
.distance.3d <- function(dx, dy, dz)
{
  sqrt(dx*dx+dy*dy+dz*dz)
}
