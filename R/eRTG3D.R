#' Extract tldCube and autodifferences functions from a track
#'
#' Get densities creates a list consisting of the 3 dimensional
#' probability distribution cube for turning angle, lift angle and step length
#' as well as the uni-dimensional distributions of the differences
#' of the turning angles, lift angles and step lengths with a lag of 1 to maintain
#' minimal level of autocorrelation in each of the terms.
#'
#' @param track a data.frame with 3 columns containing the x,y,z coordinates
#' @param heightDistEllipsoid logical: Should a distribution of the flight height over ellipsoid be extracted and later used in the sim.cond.3d()?
#' @param DEM a raster containting a digital elevation model, covering the same extent as the track
#'
#' @return A list containing the tldCube and the autodifferences functions (and additionally the height distribution function)
#' @export
#'
#' @examples
#' get.densities.3d(track, heightDist = TRUE)
get.densities.3d <- function(track, heightDistEllipsoid = TRUE, DEM = NULL)
{
  if(any(is.na(track[,1:3]))) stop("Track 'data.frame' contains NA values.")
  # track properties
  track <- track.properties.3d(track)
  turningAngle <- track$t[2:nrow(track)]; liftAngle <- track$l[2:nrow(track)]; stepLength <- track$d3d[2:nrow(track)]
  deltaTurn <- diff(turningAngle); deltaLift <- diff(liftAngle); deltaStep <- diff(stepLength)
  # probability distribution cube for turning angle, lift angle and step length
  cubeTLD <- TurnLiftStepHist(turn = turningAngle, lift = liftAngle, step = stepLength)
  # approximate the distribution of the difference in turning angle with lag 1
  autoT <- approxfun(density.default(deltaTurn))
  # approximate the distribution of the difference in lift angle with lag 1
  autoL <- approxfun(density.default(deltaLift))
  # approximate the distribution of the difference in step length with lag 1
  autoD <- approxfun(density.default(deltaStep))
  if (heightDistEllipsoid) {hDistEllipsoid <- approxfun(density.default(track$z))} else {hDistEllipsoid <- function(x){1}}
  if (!is.null(DEM)) {
    if (!.check.extent(DEM = DEM, track = track)) stop("The track is not inside the area of the digital elevation model.")
    hDistTopo <- approxfun(density.default(track$z - raster::extract(DEM, track[,1:2])))
  } else {hDistTopo <- function(x){1}}
  return(list(tldCube = cubeTLD, autoT = autoT, autoL = autoL, autoD = autoD, hDistEllipsoid = hDistEllipsoid, hDistTopo=hDistTopo))
}

#' 3 dimensional histogram
#'
#' Derives a 3 dimensional distribution of a turn angle,
#' lift angle and step length, by using the Freedman–Diaconis rule for
#' estimating the number of bins.
#'
#' @param turn numeric vector of turn angles
#' @param lift numeric vector of lift angles
#' @param step numeric vector of step lengths
#' @param printDims logical: Should dimensions of tld-Cube be messaged?
#' @param rm.zeros logical: should combinations with zero probability be removed?
#' @param maxBin numeric scalar
#'
#' @return A 3 dimensional histogram as data.frame
#' @export
#'
#' @examples
#' TurnLiftStepHist(turn, lift, step)
TurnLiftStepHist <- function(turn, lift, step, printDims = TRUE, rm.zeros = TRUE, maxBin = 25)
{
  # define based on df rule the number of bins
  # minimally 12 bins for turn angle
  nx <- min(max(floor(2 * pi / .fd.bw(turn)), 12), maxBin)
  # minimally 12 bins for lift angle
  ny <- min(max(floor(2 * pi / .fd.bw(lift)), 12), maxBin)
  # minimally 12 bins for step lengtht
  nz <- min(max(floor(max(step) / .fd.bw(step)), 12), maxBin)
  if(printDims){message("  |TLD cube dimensions: ", nx, " x ", ny, " x ", nz)}
  # create histogram
  tCuts <- .cutMidpoints(turn, nx); lCuts <- .cutMidpoints(lift, ny); dCuts <- .cutMidpoints(step, nz)
  h <- list(turn=tCuts[[1]],
            lift=lCuts[[1]],
            step=dCuts[[1]])
  h <- do.call(data.frame, h)
  h <- as.data.frame(table(h))
  # resolutions
  tRes <- tCuts[[2]]; lRes <- lCuts[[2]]; dRes <- dCuts[[2]];
  # probabilities
  colnames(h)[4] <- "prob"
  # Remove zeros if desired
  if (rm.zeros) {h <- h[!h$prob==0, ]}
  # normalize frequency to get probabilities
  h$prob <- h$prob/sum(h$prob)
  # convert factors to numeric
  h[1:3] <- lapply(h[1:3], function(x) {as.numeric(levels(x))[x]})
  return(list(values = h, tRes = tRes, lRes = lRes, dRes = dRes))
}

#' Find corresponding midpoints of breaks
#'
#' Assigns each element of a input vector the mid point of the corresponding interval.
#' Adjusted the original cut() function from the base package.
#' The midpoints are turned to a factor which facilitates the further histogram creation.
#'
#' @param x a numeric vector
#' @param breaks integer, which defines the number of breaks
#' @param rm.empty logical: should not occuring levels in the factor be removed?
#'
#' @return Returns a factor containing the midpoint
#' @export
#'
#' @examples
#' .cutMidpoints(x, breaks)
.cutMidpoints <- function(x, breaks, rm.empty=TRUE) {
  nb <- as.integer(breaks + 1)
  dx <- diff(rx <- range(x, na.rm = TRUE))
  breaks <- seq.int(rx[1L], rx[2L], length.out = nb)
  res <- median(diff(breaks))
  breaks[c(1L, nb)] <- c(rx[1L] - dx/1000, rx[2L] +
                           dx/1000)
  code <- .bincode(x, breaks, right=TRUE, include.lowest=FALSE)
  width <- diff(breaks); minBreak <- min(breaks)
  midpoints <- sapply(1:(length(breaks)-1), function(ii) {minBreak + sum(width[1:ii-1]) + width[ii]/2})
  if(rm.empty) {list(cuts = factor(midpoints[code]), res = res)}
  else {list(cuts = factor(midpoints[code], midpoints), res = res)}
}

#' Freedman–Diaconis rule
#'
#' In statistics, this rule can be used to select the size
#' of the bins to be used in a histogram.
#'
#' @param x numeric vector
#'
#' @return The bandwith
#' @export
#'
#' @examples
#' .df.bw(x)
.fd.bw <- function(x)
{
  2 * IQR(x) / (length(x) ^ (1/3))
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
#' .check.extent(DEM, track)
.check.extent <- function(DEM, track)
{
  if(!class(DEM)=="RasterLayer") stop("'DEM' is not of type 'RasterLayer'")
  if(!class(track)=="data.frame") stop("'track' is not of type 'data.frame'")
  if(any(is.na(track[,1:2]))) stop("Track 'data.frame' contains NA values.")
  e <- raster::extent(DEM)
  return(min(track[,1]) >= e[1] && max(track[,1]) <= e[2] &&
           min(track[,2]) >= e[3] && max(track[,2]) <= e[4])
}

#' Uncontidioned Empirical Random Walk (UERW) in 3D
#'
#' This function creates unconditional walks with prescribed
#' empirical properties (turning angle, lift angle and step length
#' and the auto-differences of them. It can be used for uncon-
#' ditional walks or to seed the conditional walks with
#' comparably long simulations.
#' Simulations connecting start and end points
#' with more steps than 1/10th or more of the number of steps
#' of the empirical data should rather rely on simulated
#' unconditional walks with the same properties than on
#' the empirical data (factor 1500).
#' The conditional walk connecting a given start
#' with a certain end point by a given number of
#' steps needs an attraction term (the Q probability, see below)
#' to ensure that the target is approached and hit.
#' In order to calculate the Q probability for each step
#' the distribution of turns and lifts to target and
#' the distribution of distance to target has to be knwown.
#' They can be derived from the empirical data (ideally),
#' or estimated from an unconditional process with the same properties.
#' Creates a unconditioned empirical random walk, with a specific starting point,
#' geometrically similar to the initial trajectory.
#' For a random initial heading a0 use:
#'   sample(atan2(diff(coordinates(track)[,2]), diff(coordinates(track)[,1])),1)
#'
#' @param n.locs the number of locations for the simulated track
#' @param start vector indicating the start point c(x,y,z)
#' @param a0 initial heading in radian
#' @param g0 initial gradient/polar angle in radian
#' @param densities list object returned by get.densities.3d() function
#' @param error logical: add random noise to the turn angle, lift angle and step length to account for errors measurements?
#'
#' @return A 3 dimensional trajectory in the form of a data.frame
#' @export
#'
#' @examples
#' sim.uncond.3d(n.locs, start=c(0,0,0), a0, g0, densities)
sim.uncond.3d <- function(n.locs, start=c(0,0,0), a0, g0, densities, error = TRUE)
{
  # progress bar and time
  message(paste("  |Simulate UERW with ", n.locs, " steps", sep = ""))
  start.time <- Sys.time()
  pb <- txtProgressBar(min = 0, max = n.locs, style = 3)
  ui <- floor(n.locs/20)+1
  # get coordinates of the tldCube
  ts <- densities$tldCube$values$turn
  ls <- densities$tldCube$values$lift
  ds <- densities$tldCube$values$step
  # get probs for each turn-lift-distance combination
  tldProbs <- densities$tldCube$values$prob
  sCond <- sample(1:nrow(densities$tldCube$values), 1, prob=tldProbs)
  # "x" "y" "z" "a" "g" "t" "l" "d" "p"
  # "1" "2" "3" "4" "5" "6" "7" "8" "9"
  RTG <- matrix(0, n.locs, 9)
  RTG[1,] <- c(start[1], start[2], start[3], a0, g0, ts[sCond], ls[sCond], ds[sCond], NA)
  # Create random noise if error is TRUE (uniform distributed)
  if (error) {
    tShift <- runif(n.locs, -densities$tldCube$tRes / 2, densities$tldCube$tRes / 2)
    lShift <- runif(n.locs, -densities$tldCube$lRes / 2, densities$tldCube$lRes / 2)
    dShift <- runif(n.locs, -densities$tldCube$dRes / 2, densities$tldCube$dRes / 2)
  } else {
    tShift <- lShift <- dShift <- numeric(n.locs)
  }
  for (i in 2:n.locs)
  {
    # get influence of current autodifferences
    atProbs <- densities$autoT(RTG[i-1, 6] - ts)
    alProbs <- densities$autoL(RTG[i-1, 7] - ls)
    adProbs <- densities$autoD(RTG[i-1, 8] - ds)
    atProbs[is.na(atProbs)] <- 0
    alProbs[is.na(alProbs)] <- 0
    adProbs[is.na(adProbs)] <- 0
    atProbs <- atProbs / sum(atProbs)
    alProbs <- alProbs / sum(alProbs)
    adProbs <- adProbs / sum(adProbs)
    # multiply and take the third squareroot
    pProbs <- tldProbs * (atProbs * alProbs * adProbs)^(1/3)
    # limit gradient (0-pi)
    gAll <- (RTG[i-1, 5] + ls + lShift[i])
    pProbs <- pProbs * as.numeric((gAll > 0 & gAll < pi))
    # sample on TurnLiftStepHist = tldCube and add shifts
    rP <- sample(1:nrow(densities$tldCube$values), size = 1, prob = pProbs)
    t <- ts[rP] + tShift[i]
    l <- ls[rP] + lShift[i]
    d <- ds[rP] + dShift[i]
    p <- pProbs[rP]
    # absolute spherical orientation, wrap angles around -pi-0 & 0-pi
    a <- .wrap(RTG[i-1, 4] + t)
    g <- .wrap(RTG[i-1, 5] + l)
    # new coordinates of the next step
    x <- (d * sin(g) * cos(a)) + RTG[i-1, 1]
    y <- (d * sin(g) * sin(a)) + RTG[i-1, 2]
    z <- (d * cos(g)) + RTG[i-1, 3]
    # "x" "y" "z" "a" "g" "t" "l" "d" "p"
    RTG[i,] <- c(x, y, z, a, g, t, l, d, p)
    # update progress bar
    if ((i %% ui) == 0) {setTxtProgressBar(pb, i)}
  }
  rownames(RTG) <- c()
  colnames(RTG) <- c("x", "y", "z", "a", "g", "t", "l", "d", "p")
  # close progress bar
  setTxtProgressBar(pb, n.locs)
  close(pb)
  message(paste("  |Runtime: ", round(as.numeric(Sys.time()) - as.numeric(start.time), 2), " secs", sep = ""))
  return(as.data.frame(RTG))
}

#' Q probabilities for n steps
#'
#'  Calculates the Q probability, representing the pull to
#' the target. The number of steps on which the Q prob will be
#' quantified is number of total segments less than one
#' (the last step is defined by the target itself).
#'
#' @param sim the result of simm.uncond.3d(), or a data frame with at least
#'     x,y,z-coordinates, the arrival azimuth and the arrival gradient.
#' @param n.locs number of total segments to be modelled,
#'     the length of the desired conditioned empirical random walk
#' @param multicore logical: run computations in parallel (n-1 cores)?
#'
#' @return A list containing the Q - tldCubes for every step
#' @export
#'
#' @examples
#' qProb.3d(sim, n.locs)
qProb.3d <- function(sim, n.locs, multicore = FALSE)
{
  if (multicore) {
    if(.Platform$OS.type == "unix") {return(.qProb.3d.unix(sim, n.locs,  multicore = TRUE))}
    if(.Platform$OS.type == "windows") {return(suppressWarnings(.qProb.3d.windows(sim, n.locs)))}
  } else {
    start.time <- Sys.time()
    message(paste("  |Extracting Q probabilities for ", n.locs, " steps", sep = ""))
    # progress bar
    pb <- txtProgressBar(min = 0, max = 18, style = 3)
    # steps minus 2
    nSteps <- n.locs - 2
    # turning angles to target as a function of number of steps
    tList <- lapply(1:nSteps, function(x) .wrap(atan2(diff(sim$y, lag = x),
                                                                  diff(sim$x, lag = x)) - sim$a[1:(length(sim$a) - x)]))
    setTxtProgressBar(pb, 3)
    # lift angles to target as a function of number of steps
    lList <- lapply(1:nSteps, function(x) .wrap(atan2(sqrt(diff(sim$x, lag = x) ^ 2 + diff(sim$y, lag = x) ^ 2),
                                                                  diff(sim$z, lag = x)) - sim$g[1:(length(sim$g) - x)]))
    setTxtProgressBar(pb, 6)
    # calculate distance to target as a function of number of steps
    dList <- lapply(1:nSteps, function(x) sqrt(diff(sim$x, lag = x) ^ 2
                                                           + diff(sim$y, lag = x) ^ 2
                                                           + diff(sim$z, lag = x) ^ 2))
    setTxtProgressBar(pb, 9)
    # the Qprob is thinned to the lag that suggests breaking off of the autocorrelation
    # of the turning angle to target, the lift angle to target and the distance to target
    # for the relevant number of steps. This is mainly to reduce redundancy mainly
    # introduced by the sliding window approach adopted in estimating the relationships
    k <- suppressWarnings(cbind(unlist(lapply(lapply(lapply(lapply(lapply(tList, acf, lag.max=nSteps, plot = FALSE, mc.cores = nCores),
                                                                   '[[', 'acf'), '<', .05),
                                                     which), head, 1)) - 1,
                                unlist(lapply(lapply(lapply(lapply(lapply(lList, acf, lag.max=nSteps, plot = FALSE, mc.cores = nCores),
                                                                   '[[', 'acf'), '<', .05),
                                                     which), head, 1)) - 1,
                                unlist(lapply(lapply(lapply(lapply(lapply(dList, acf, lag.max=nSteps, plot = FALSE, mc.cores = nCores),
                                                                   '[[', 'acf'), '<', .05),
                                                     which), head, 1)) - 1))
    kk <- apply(k,1,max)
    setTxtProgressBar(pb, 14)
    tList <-mapply('[',tList,mapply(seq, 1, lapply(tList, length), by = kk))
    lList <-mapply('[',lList,mapply(seq, 1, lapply(lList, length), by = kk))
    dList <-mapply('[',dList,mapply(seq, 1, lapply(dList, length), by = kk))
    # Use multicore to speed the calculations up
    cubeList <- rev(lapply(1:nSteps, function(x) TurnLiftStepHist(turn=tList[[x]], lift=lList[[x]], step=dList[[x]], printDims = FALSE, rm.zeros = TRUE)))
    # complete progress bar and close
    setTxtProgressBar(pb, 18)
    close(pb)
    message("  |Minimum number of independent estimates: ", min(unlist(lapply(dList, length))), " for step ", which.min(unlist(lapply(dList, length))), ".")
    message(paste("  |Runtime: ", round(as.numeric(Sys.time()) - as.numeric(start.time), 2), " secs", sep = ""))
    return(cubeList)
  }
}

#' Conditioned Empirical Random Walk (CERW) in 3D
#'
#' Creates a conditioned empirical random walk, with a specific starting and ending point,
#' geometrically similar to the initial trajectory
#' (extractMethod: raster overlay method can take "simple" or "bilinear")
#'
#' @param n.locs length of the trajectory in locations
#' @param start numeric vector of length 3 with the coordinates of the start point
#' @param end numeric vector of length 3 with the coordinates of the end point
#' @param a0 initial incoming heading in radian
#' @param g0 initial incoming gradient/polar angle in radian
#' @param densities list object returned by get.densities.3d() function
#' @param qProbs list object returned by qProb.3d() function
#' @param error logical: add random noise to the turn angle, lift angle and step length to account for errors measurements?
#' @param DEM raster layer containing a digital elevation model, covering the area between start and end point
#' @param BG a background raster layer that can be used to inform the choice of steps
#'
#' @return A trajectory in the form of data.frame
#' @export
#'
#' @examples
#' sim.cond.3d(n.locs, start, end=start, a0, g0, densities, qProbs)
sim.cond.3d <- function(n.locs, start=c(0,0,0), end=start, a0, g0, densities, qProbs, error = FALSE, DEM = NULL, BG = NULL)
{
  start.time <- Sys.time()
  if(!is.null(DEM)) {
    if(!class(DEM)=="RasterLayer") stop("'DEM' is not of type 'RasterLayer'")
    if (!.check.extent(DEM = DEM, track = data.frame(rbind(start, end)))) stop("The track is not inside the area of the digital elevation model 'DEM'.")
  }
  if(!is.null(BG)) {
    if(!class(BG)=="RasterLayer") stop("'BG' is not of type 'RasterLayer'")
    if (!.check.extent(DEM = BG, track = data.frame(rbind(start, end)))) stop("The track is not inside the area of back ground layer 'BG'.")
  }
  # progress bar and time
  message(paste("  |Simulate CERW with ", n.locs, " steps", sep = ""))
  pb <- txtProgressBar(min = 0, max = n.locs-2, style = 3)
  ui <- floor(n.locs/20)+1
  # replace the probability distribution for step length 1 by the one from
  # the qProbs since that one relies on more samples derived from sim
  densities[[1]] <- tail(qProbs,1)[[1]]
  # get the coordinates of the step length and turning angle bin centres
  names(start) <- c("x", "y", "z")
  names(end) <- c("x", "y", "z")
  # get coordinates of the tldCube
  ts <- densities$tldCube$values$turn
  ls <- densities$tldCube$values$lift
  ds <- densities$tldCube$values$step
  # get probs for each combination
  tldProbs <- densities$tldCube$values$prob
  # sample one randomly to set the initial conditions
  # for the previous to first turn and previous to first step
  # for the start point, as this is needed to inform the auto-difference
  # likelihood
  sCond <- sample(1:nrow(densities$tldCube$values), 1, prob=tldProbs)
  # "x" "y" "z" "a" "g" "t" "l" "d" "p"
  # "1" "2" "3" "4" "5" "6" "7" "8" "9"
  RTG <- matrix(0, n.locs, 9)
  RTG[1, ] <- c(start[1], start[2], start[3], a0, g0, ts[sCond], ls[sCond], ds[sCond], NA)
  # Create random noise if error is TRUE
  if (error) {
    tShift <- runif(n.locs - 2, -densities$tldCube$tRes / 2, densities$tldCube$tRes / 2)
    lShift <- runif(n.locs - 2, -densities$tldCube$lRes / 2, densities$tldCube$lRes / 2)
    dShift <- runif(n.locs - 2, -densities$tldCube$dRes / 2, densities$tldCube$dRes / 2)
  } else {
    tShift <- lShift <- dShift <- numeric(n.locs - 2)
  }
  # start creating the track step for step
  for (i in 1:(n.locs - 2))
  {
    # get the auto-difference probability for turning angle
    atProbs <- densities$autoT(RTG[i, 6]- ts)
    # get the auto-difference probability for lift angle
    alProbs <- densities$autoL(RTG[i, 7] - ls)
    # get the auto-difference probability for step length
    adProbs <- densities$autoD(RTG[i, 8] - ds)
    # set NAs to zero probability
    atProbs[is.na(atProbs)] <- 0
    alProbs[is.na(alProbs)] <- 0
    adProbs[is.na(adProbs)] <- 0
    # standardize the probabilities to sum to one
    atProbs <- atProbs / sum(atProbs)
    alProbs <- alProbs / sum(alProbs)
    adProbs <- adProbs / sum(adProbs)
    # calculate the probability to make a step forward. The auto-difference probabilities are
    # calculated as one jointly contributing probability and therefore square rooted befor
    # multiplication with the two dimensional probability distribution
    P <- (tldProbs) * (atProbs * alProbs * adProbs)^(1/3)
    # calculate the azimuth
    a <- .wrap(RTG[i, 4] + ts + tShift[i])
    # calculate the gradient
    g <- .wrap(RTG[i, 5] + ls + lShift[i])
    # convert the coordinates from step length turning angle dimension
    x1 <- ((ds + dShift[i]) * sin(g) * cos(a)) + RTG[i, 1]
    y1 <- ((ds + dShift[i]) * sin(g) * sin(a)) + RTG[i, 2]
    z1 <- ((ds + dShift[i]) * cos(g)) + RTG[i, 3]
    # calculate the distances of the cell centers in the spatial domain
    # to the target (last location of the empirical track)
    endD <- as.numeric(sqrt((end[1] - x1) ^ 2 + (end[2] - y1) ^ 2 + (end[3] - z1) ^ 2))
    # calculate the azimuth of the cell centres to the target and substract from it the direction of arrival
    # resulting in turning angle towards target
    endT <- as.numeric(.wrap(atan2(as.numeric(end[2] - y1), as.numeric(end[1] - x1)) - a))
    # calculate the gradient of the possibilite steps to the target and substract from it the angle of arrival
    # resulting in turning angle towards target
    endL <- as.numeric(.wrap(atan2(as.numeric(sqrt((end[1] - x1) ^ 2 + (end[2] - y1) ^ 2)), as.numeric(end[3] - z1))) - g)
    # get the probabilities of making it distance and turning angle wise
    # which is derived from the two dimensional probability distribution for the
    # appropriate step being modelled
    # get possible coordinates
    qCube <- qProbs[[i]]
    tVal <- unique(qCube$values$turn)
    lVal <- unique(qCube$values$lift)
    dVal <- unique(qCube$values$step)
    # find closest coordinates
    tCoords <- unlist(lapply(endT, function(x) tVal[which.min(abs(tVal-x))]))
    lCoords <- unlist(lapply(endL, function(x) lVal[which.min(abs(lVal-x))]))
    dCoords <- unlist(lapply(endD, function(x) dVal[which.min(abs(dVal-x))]))
    # extract Q
    Q <- unlist(lapply(1:length(tCoords), function(x){
      test <- (qCube$values$turn == tCoords[x] & qCube$values$lift == lCoords[x] & qCube$values$step == dCoords[x]);
      if(any(test==TRUE)) {return(qCube$values$prob[test])} else {return(0)}
    }))
    # the overall probability is the product of the probability
    # of making a step forward and the probability of making it to the
    # target. The weight of the target probability needs to be adjusted
    # by division by end distance, because the number of cells to choose from
    # are increasing with distance to target, which needs to be accounted
    # for prior to sampling based on overall probability
    Probs <- P * Q / endD
    # limit gradient to 0 and pi
    Probs <- Probs * as.numeric((g > 0) & (g < pi))
    # Account for probable flight height, if a DEM is provided the relative flight height is taken
    # Otherwise only the absoulte ellipsoid height.
    if(!is.null(DEM))
    {
      surface <- raster::extract(DEM, cbind(x1, y1))
      demP <-  densities$hDistTopo(z1 - surface) * as.numeric(z1 >= surface)
      demP[is.na(demP)] <- 0
      demP <- demP / sum(demP)
      hProb <- densities$hDistEllipsoid(z1)
      hProb[is.na(hProb)] <- 0
      hProb <- hProb / sum(hProb)
      Probs <- Probs * sqrt(demP * hProb)
    } else {
      hProb <- densities$hDistEllipsoid(z1)
      hProb[is.na(hProb)] <- 0
      hProb <- hProb / sum(hProb)
      Probs <- Probs * hProb
    }
    # use the coordinates of the spatial grid for which the probabilities are calculated
    # and use it to overlay a background raster for example to avoid walks in certain areas
    if(!is.null(BG))
    {
      bgP <- raster::extract(BG, cbind(x1, y1))
      Probs <- Probs * bgP
    }
    # make sure we have no missing nor negative probabilities
    Probs[is.na(Probs)] <- 0
    Probs[Probs <= 0] <- 0
    # check whether the run might have ended up in a dead-end,
    # which will set the zero probability status to TRUE
    if(all(Probs==0)){
      RTG <- NULL
      close(pb)
      message(paste("  |Runtime: ", round(as.numeric(Sys.time()) - as.numeric(start.time), 2), " secs", sep = ""))
      warning("Dead end encountered.")
      return(RTG)
    }else{
      # draw a point randomly based on the probability
      rP <- sample.int(nrow(densities$tldCube$values), size = 1, prob = Probs)
      # "x" "y" "z" "a" "g" "t" "l" "d" "p"
      # "1" "2" "3" "4" "5" "6" "7" "8" "9"
      RTG[i + 1, ] <- c(x1[rP], y1[rP], z1[rP], a[rP], g[rP], ts[rP], ls[rP], ds[rP], Probs[rP])
      # update progress bar
      if ((i %% ui) == 0) {setTxtProgressBar(pb, i)}
    }
  }
  # the track is forced to target location and the appropriate distance is added
  RTG[1, 8] <- NA
  RTG[n.locs,] <- c(end[1], end[2], end[3], NA, NA, NA, NA, NA, NA)
  RTG[n.locs, 8] <- sqrt((RTG[n.locs, 1] - RTG[n.locs-1, 1])^2 + 
                         (RTG[n.locs, 2] - RTG[n.locs-1, 2])^2 + 
                         (RTG[n.locs, 3] - RTG[n.locs-1, 3])^2)
  rownames(RTG) <- c()
  colnames(RTG) <- c("x", "y", "z", "a", "g", "t", "l", "d", "p")
  # close progress bar
  setTxtProgressBar(pb, i)
  close(pb)
  message(paste("  |Runtime: ", round(as.numeric(Sys.time()) - as.numeric(start.time), 2), " secs", sep = ""))
  return(as.data.frame(RTG))
}


#' Conditioned Empirical Random Walks (CERW) in 3D
#'
#' Creates n conditioned empirical random walks, with a specific starting and ending point,
#' geometrically similar to the initial trajectory
#'
#' @param n.sim number of CERWs to simulate
#' @param n.locs length of the trajectory in locations
#' @param start numeric vector of length 3 with the coordinates of the start point
#' @param end numeric vector of length 3 with the coordinates of the end point
#' @param a0 initial incoming heading in radian
#' @param g0 initial incoming gradient/polar angle in radian
#' @param densities list object returned by get.densities.3d() function
#' @param qProbs list object returned by qProb.3d() function
#' @param error logical: add random noise to the turn angle, lift angle and step length to account for errors measurements?
#' @param DEM raster layer containing a digital elevation model, covering the area between start and end point
#' @param BG a background raster layer that can be used to inform the choice of steps
#' @param multicore logical: run computations in parallel (n-1 cores)?
#'
#' @return
#' A list containing the CERWs or NULLs if dead ends have been encountered.
#' @export
#'
#' @examples
#' n.sim.cond.3d(n.sim, n.locs, start = c(0,0,0), end=start, a0, g0, densities, qProbs)
n.sim.cond.3d <- function(n.sim, n.locs, start = c(0,0,0), end=start, a0, g0, densities, qProbs, error = FALSE, multicore = FALSE, DEM = NULL, BG = NULL)
{
  start.time <- Sys.time()
  n.sim <- round(n.sim)
  if (n.sim <= 1) {return(sim.cond.3d(n.locs=n.locs, start=start, end=end, a0=a0, g0=g0, densities=densities, qProbs=qProbs, error=error, DEM=DEM, BG=BG))}
  message(paste("  |Simulate ", n.sim ," CERWs with ", n.locs, " steps", sep = ""))
  if (multicore) {
    # nCores <- parallel::detectCores()-1
    if(.Platform$OS.type == "unix") {
      .n.sim.cond.3d.unix(n.locs=n.locs, start=start, end=end, a0=a0, g0=g0, densities=densities, qProbs=qProbs, error=error, DEM=DEM, BG=BG)
      # message(paste("  |Running on nCores = ", nCores, sep=""))
      # message("  |...")
      # cerwList <- parallel::mclapply(X = 1:n.sim, FUN = function(x){sim.cond.3d(n.locs, start, end, a0, g0, densities, qProbs, error, DEM, BG)}, mc.cores = nCores)
    }
    if(.Platform$OS.type == "windows") {
      .n.sim.cond.3d.windows(n.locs=n.locs, start=start, end=end, a0=a0, g0=g0, densities=densities, qProbs=qProbs, error=error, DEM=DEM, BG=BG)
      #stop("Parallel version not yet supported on Windows. Please set 'multicore' to 'FALSE' or change to a unix system.")
      }
  } else {
    cerwList <- suppressMessages(lapply(X = 1:n.sim, FUN = function(x){sim.cond.3d(n.locs=n.locs, start=start, end=end, a0=a0, g0=g0, densities=densities, qProbs=qProbs, error=error, DEM=DEM, BG=BG)}))
  }
  message(paste("  |Runtime: ", round(as.numeric(Sys.time()) - as.numeric(start.time), 2), " secs", sep = ""))
  return(cerwList)
}

#' Wraps angles to [0,2pi] degrees
#'
#' @param x the angle in radians
#'
#' @return The wrapped angle
#' @export
#'
#' @examples
#' .wrap(x)
.wrap <- function(x)
{
  (x + pi) %% (2 * pi) - pi
}

#' Simulation of a three dimensional Correlated Random Walk
#'
#' @param nStep the number of steps of the simulated trajectory
#' @param rTurn the correlation on the turn angle
#' @param rLift the correlation of the lift angle
#' @param meanStep the mean step length
#' @param start a vector of length 3 containing the coordinates of the startpoint of the trajectory
#'
#' @return A trajectory in the form of data.frame
#' @export
#'
#' @examples
#' sim.crw.3d(nStep, rTurn, rLift, meanStep, start = c(0,0,0))
sim.crw.3d <- function(nStep, rTurn, rLift, meanStep, start = c(0,0,0))
{
  # correlated angles and distance
  t <- CircStats::rwrpnorm(n = nStep - 2, mu = 0, rho = rTurn)
  a <- .wrap(cumsum(c(runif(1, 0, 2 * pi), t)))
  l <- CircStats::rwrpnorm(n = nStep - 2, mu = 0, rho = rLift)
  g <- .wrap(cumsum(c(runif(1, 0, pi), l)))
  f <- abs(scale(CircStats::rwrpnorm(n = nStep - 1, mu = 0, rho = (rTurn+rLift)/2))[,1])
  d <- rep(meanStep, nStep-1) * f
  # deltas in all 3 directions
  dx <- (d * sin(g) * cos(a))
  dy <- (d * sin(g) * sin(a))
  dz <- (d * cos(g))
  # generate track
  t <- data.frame(
    x = cumsum(c(start[1], dx)),
    y = cumsum(c(start[2], dy)),
    z = cumsum(c(start[3], dz))
  )
  return(t)
}

#' Track properties of a 3D track
#'
#' Returns the properties (distances, azimut, polar angle,
#' turn angle & lift angle) of a track in three dimensions.
#'
#' @param track data.frame with x,y,z coordinates
#'
#' @return The data.frame with track properties
#' @export
#'
#' @examples
#' track.properties.3d(track)
track.properties.3d <- function(track)
{
  if(!is.data.frame(track)) stop("Input is not of type 'data.frame'.")
  if(any(is.na(track[,1:3]))) stop("Track 'data.frame' contains NA values.")
  # spatial coordinates
  x <- track[,1]; y <- track[,2]; z <- track[,3]
  # distance covered per axis and step
  dx <- c(NA, diff(x)); dy <- c(NA, diff(y)); dz <- c(NA, diff(z))
  # spherical coordinates
  d2d <- .distance.2d(dx, dy)
  d3d <- .distance.3d(dx, dy, dz)
  a <- .get.azimut(dx, dy)
  g <- .get.polar(d2d, dz)
  # guess the initial azimut/heading & polar angle/gradient (alt: sample(a,1))
  a[1] <- mean(a, na.rm=TRUE)
  g[1] <- mean(g, na.rm=TRUE)
  # Turn angle and lift angle
  t <- c(NA, diff(a))
  l <- c(NA, diff(g))
  return(data.frame(x, y, z, a, g, dx, dy, dz, d2d, d3d, t, l))
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
#' .get.azimut(dx, dy)
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
#' .get.polar(d, dz)
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
#' .distance.2d(dx, dy)
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
#' .distance.3d(dx, dy, dz)
.distance.3d <- function(dx, dy, dz)
{
  sqrt(dx*dx+dy*dy+dz*dz)
}

#' Internally verification of the simulated track
#'
#' Uses two-sample Kolmogorov-Smirnov test to compare the geometric characteristics of the orginal track
#' with the characteristics of the simulated track.
#'
#' @param track1 data.frame with x,y,z coordinates of the original track
#' @param track2 data.frame or list of data.frames with x,y,z coordinates of the simulated track
#' @param alpha scalar: significance level, default alpha = 0.05
#' @param plotDensities logical: plot the densites of turn angle, lift angle and step length of the two tracks?
#'
#' @return Test objects of the 6 two-sample Kolmogorov-Smirnov test conducted.
#' @export
#'
#' @examples
#' test.verification.3d(track1, track2)
test.verification.3d <- function(track1, track2, alpha = 0.05, plotDensities = FALSE)
{
  message("  |*** Two-sample Kolmogorov-Smirnov test ***")
  track1 <- track.properties.3d(track1)[2:nrow(track1), ]
  t1 <- track1$t; l1 <- track1$l; d1 <- track1$d3d;
  diffT1 <- diff(track1$t); diffL1 <- diff(track1$l); diffD1 <- diff(track1$d3d);
  if(is.data.frame(track2)){
    track2 <- track.properties.3d(track2)[2:nrow(track2), ]
    t2 <- track2$t; l2 <- track2$l; d2 <- track2$d3d;
    diffT2 <- diff(track2$t); diffL2 <- diff(track2$l); diffD2 <- diff(track2$d3d);
  } else {
    track2 <- filter.dead.ends(track2)
    track2 <- lapply(track2, function(x){track.properties.3d(x)[2:nrow(x), ]})
    diffTrack2 <- lapply(track2, function(x){data.frame(diffT = diff(x$t), diffL = diff(x$l), diffD = diff(x$d3d))})
    track2 <- do.call("rbind", track2)
    diffTrack2 <- do.call("rbind", diffTrack2)
    t2 <- track2$t; l2 <- track2$l; d2 <- track2$d3d;
    diffT2 <- diffTrack2$diffT; diffL2 <- diffTrack2$diffL; diffD2 <- diffTrack2$diffD;
  }
  # turn
  turnT <- suppressWarnings(ks.test(t1, t2, alternative = "two.sided"))
  diffTurnT <- suppressWarnings(ks.test(diffT1, diffT2, alternative = "two.sided"))
  # lift
  liftT <- suppressWarnings(ks.test(l1, l2, alternative = "two.sided"))
  diffLiftT <- suppressWarnings(ks.test(diffL1, diffL2, alternative = "two.sided"))
  # step
  stepT <- suppressWarnings(ks.test(d1, d2, alternative = "two.sided"))
  diffStepT <- suppressWarnings(ks.test(diffD1, diffD2, alternative = "two.sided"))
  message("  |H0: Probability distributions do not differ significantly")
  message("  |H1: Probability distributions differ significantly")
  message(paste("  |Turn angle  – ", .test2text(turnT, alpha), ", autodifferences – ", .test2text(diffTurnT, alpha), sep=""))
  message(paste("  |Lift angle  – ", .test2text(liftT, alpha), ", autodifferences – ", .test2text(diffLiftT, alpha), sep=""))
  message(paste("  |Step length – ", .test2text(stepT, alpha), ", autodifferences – ", .test2text(diffStepT, alpha), sep=""))
  if (plotDensities) {
    suppressWarnings(plot2d.multiplot(
      .plot2d.density(t1, t2, titleText = "Turn angle"),
      .plot2d.density(l1, l2, titleText = "Lift angle"),
      .plot2d.density(d1, d2, titleText = "Step length"),
      cols = 1
    ))
  }
  return(list(turnT, liftT, stepT, diffTurnT, diffLiftT, diffStepT))
}

#' Extract test results as string
#'
#' @param test object of type 'htest'
#' @param alpha scalar: significance level, default alpha = 0.05
#'
#' @return A character describing the results.
#' @export
#'
#' @examples
#' .test2text(test, alpha)
.test2text <- function(test, alpha)
{
  p <- test$p.value
  paste("p-value: ", round(p,3) , if(p<alpha){
    paste(" < ", alpha, ", *H1*", sep = "")
  } else {
    paste(" > ", alpha, ", *H0*", sep = "")
  }, sep = "")
}

#' Function to filter out tracks that have found a dead end (=NULL)
#'
#' @param cerwList list of data.frames and NULL entries
#'
#' @return A list that is only containing valid tracks.
#' @export
#'
#' @examples
#' filter.dead.ends(cerwList)
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

#' Crops the DEM to the extent of the track with a buffer
#'
#' @param DEM a raster containing a digital elevation model, covering the extent as the track
#' @param track data.frame with x,y,z coordinates of the original track
#' @param buffer bufferwith, by default set to 100
#'
#' @return A the cropped digital elevation model as a raster layer.
#' @export
#'
#' @examples
#' dem.track.extent(DEM, track)
dem.track.extent <- function(DEM, track, buffer=100)
{
  if (!.check.extent(DEM = DEM, track = track)) stop("The track is not inside the area of the digital elevation model.")
  return(raster::crop(dem, extent(min(track$x)-buffer, max(track$x)+buffer, min(track$y)-buffer, max(track$y)+buffer)))
}

#' Test the functionality of the eRTG3D
#'
#' The test simulates a CRW with given parameters and reconstructs it by using the eRTG3D
#'
#' @param multicore logical: test with multicore?
#' @param returnResult logical: return tracks generated?
#' @param plot2d logical: plot tracks on 2d plane?
#' @param plot3d logical: plot tracks in 3D?
#'
#' @return A list containing the original CRW and the simulated track (CERW).
#' @export
#'
#' @examples
#' test.eRTG3D.3d()
test.eRTG.3d <- function(multicore = FALSE, returnResult = FALSE, plot2d = FALSE, plot3d = FALSE)
{
  message("  |*** Testing eRTG3D ***")
  set.seed(123)
  nStep <- 25
  crw <- track.properties.3d(
    sim.crw.3d(nStep = nStep, rTurn = 0.99, rLift = 0.99, meanStep = 1, start = c(0, 0, 10)))
  D <- get.densities.3d(crw, heightDistEllipsoid = FALSE)
  uerw <- sim.uncond.3d(nStep*1500, start = c(crw$x[1],crw$y[1],crw$z[1]),
                        a0 = crw$a[1], g0 = crw$g[1], densities = D)
  tests.uerw <- test.verification.3d(crw, uerw, alpha = 0.05, plotDensities = FALSE)
  Q <- qProb.3d(uerw, nStep, multicore = multicore)
  cerw <- sim.cond.3d(nStep, start=c(crw$x[1],crw$y[1],crw$z[1]), end=c(crw$x[nStep],crw$y[nStep],crw$z[nStep]),
                      a0 = crw$a[1], g0 = crw$g[1], densities=D, qProbs=Q)
  tests.cerw <- test.verification.3d(crw, cerw, alpha = 0.05, plotDensities = FALSE)
  message("  |*** Test passed successfully ***")
  if(plot2d){plot2d(crw, cerw)}
  if(plot3d){plot3d(crw, cerw)}
  if(returnResult){return(list(crw = crw, cerw = cerw))}
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
#' @param multicore logical: run calculations on multiple cores?
#' @param filterDeadEnds: logical: remove tracks (='NULL') that ended in a dead end?
#' @param plot2d logical: plot tracks on 2d plane?
#' @param plot3d logical: plot tracks in 3D?
#'
#' @return A list or data.frame containing the simulated track(s) (CERW).
#' @export
#'
#' @examples
#' reproduce.track.3d(track)
reproduce.track.3d <- function(track, n.sim = 1, multicore = FALSE, error = TRUE, DEM = NULL, BG = NULL,  plot2d = FALSE, plot3d = FALSE, filterDeadEnds = TRUE)
{
  n.locs <- nrow(track)
  track <- track.properties.3d(track)
  if (n.locs>1500) stop("Track is to long (>1500 steps).")
  D <- get.densities.3d(track, heightDistEllipsoid = TRUE, DEM = DEM)
  uerw <- sim.uncond.3d(n.locs*1500, start = c(track$x[1],track$y[1],track$z[1]),
                        a0 = track$a[1], g0 = track$g[1], densities = D, error = error)
  Q <- qProb.3d(uerw, n.locs, multicore = multicore)
  cerwList <- suppressWarnings(n.sim.cond.3d(n.sim = n.sim, n.locs <- n.locs, start=c(track$x[1],track$y[1],track$z[1]), end=c(track$x[n.locs],track$y[n.locs],track$z[n.locs]),
                                             a0 = track$a[1], g0 = track$g[1], densities=D, qProbs=Q, error = error, multicore = multicore, DEM = DEM, BG = BG))
  if(filterDeadEnds){cerwList <- filter.dead.ends(cerwList)}
  if(plot2d){plot2d(origTrack = track, cerwList = cerwList, DEM = DEM)}
  if(plot3d){plot3d(origTrack = track, cerwList = cerwList, surface = TRUE, DEM = DEM)}
  return(cerwList)
}
