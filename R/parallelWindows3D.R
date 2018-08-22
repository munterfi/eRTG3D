#' Parallel Q probabilities for n steps on windows
#'
#' Calculates the Q probability, representing the pull to
#' the target. The number of steps on which the Q prob will be
#' quantified is number of total segments less than one
#' (the last step is defined by the target itself).
#'
#' @param sim the result of simm.uncond.3d(), or a data frame with at least
#'     x,y,z-coordinates, the arrival azimuth and the arrival gradient.
#' @param n.locs number of total segments to be modelled,
#'     the length of the desired conditioned empirical random walk
#' @param maxBin numeric scalar, maximum number of bins per dimension of the tld-cube (\link[eRTG3D]{turnLiftStepHist})
#'
#' @return A list containing the Q - tldCubes for every step
#' @export
#'
#' @examples
#' qProb.3d.windows(sim, n.locs)
#' @noRd
.qProb.3d.windows <- function(sim, n.locs, maxBin = 25)
{
  start.time <- Sys.time()
  # set up cluster (parallel, doParallel and later plyr)
  nCores <- parallel::detectCores()-1
  message(paste("  |Extracting Q probabilities for ", n.locs, " steps (Parallel on nCores = ", nCores, ")", sep = ""))
  # progressbar not working with plyr::alply, .progress = "text"
  message("  |...")
  # set up cluster
  cl <- parallel::makeCluster(nCores)
  doParallel::registerDoParallel(nCores)
  opts <- list(preschedule=FALSE)
  # define variables and functions needed later to pass them to the cluster
  sim <- sim; n.locs <- n.locs; maxBin <- maxBin
  wrap <- function(x) {(x + pi) %% (2 * pi) - pi}
  turnLiftStepHist <- function(turn, lift, step, printDims = TRUE, rm.zeros = TRUE, maxBin = 25)
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
    tCuts <- cutMidpoints(turn, nx); lCuts <- cutMidpoints(lift, ny); dCuts <- cutMidpoints(step, nz)
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
  
  cutMidpoints <- function(x, breaks, rm.empty=TRUE) {
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
  fd.bw <- function(x){2 * IQR(x) / (length(x) ^ (1/3))}
  # steps minus 2
  nSteps <- n.locs - 2
  sim <- track.properties.3d(sim)
  # lift angles to target as a function of number of steps
  cubeList <- plyr::alply(1:nSteps, function(x) {
    # turn angle, lift angles and distance to target as a function of number of steps
    t <- wrap(atan2(diff(sim$y, lag = x), diff(sim$x, lag = x)) - sim$a[1:(length(sim$a) - x)])
    l <- wrap(atan2(sqrt(diff(sim$x, lag = x) ^ 2 + diff(sim$y, lag = x) ^ 2),
                     diff(sim$z, lag = x)) - sim$g[1:(length(sim$g) - x)])
    d <- sqrt(diff(sim$x, lag = x) ^ 2 + diff(sim$y, lag = x) ^ 2 + diff(sim$z, lag = x) ^ 2)
    # the Qprob is thinned to the lag that suggests breaking off of the autocorrelation
    # of the turning angle to target, the lift angle to target and the distance to target
    # for the relevant number of steps. This is mainly to reduce redundancy mainly
    # introduced by the sliding window approach adopted in estimating the relationships
    k <- max(head(which(acf(t, lag.max = nSteps, plot = FALSE)$acf < 0.05),1)-1,
             head(which(acf(l, lag.max = nSteps, plot = FALSE)$acf < 0.05),1)-1,
             head(which(acf(d, lag.max = nSteps, plot = FALSE)$acf < 0.05),1)-1)
    t <- t[seq(1, length(t), by = k)]
    l <- l[seq(1, length(l), by = k)]
    d <- d[seq(1, length(d), by = k)]
    # get stepTurnLiftHistograms
    return(turnLiftStepHist(turn=t, lift=l, step=d, printDims = FALSE, rm.zeros = TRUE, maxBin = maxBin))
  }, .margins = 1, .parallel = TRUE, .paropts = list(.options.snow=opts))
  # stop cluster
  parallel::stopCluster(cl)
  message(paste("  |Runtime: ", round(as.numeric(Sys.time()) - as.numeric(start.time), 2), " secs", sep = ""))
  return(rev(cubeList))
}

#' Parallel computation of n Conditioned Empirical Random Walks (CERW) in 3-D on Windows
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
#'
#' @return
#' A list containing the CERWs or NULLs if dead ends have been encountered.
#' @export
#'
#' @examples
#' .n.sim.cond.3d.windows(n.sim, n.locs, start = c(0,0,0), end=start, a0, g0, densities, qProbs)
#' @noRd
.n.sim.cond.3d.windows <- function(n.sim, n.locs, start = c(0,0,0), end=start, a0, g0, densities, qProbs, error = FALSE, DEM = NULL, BG = NULL)
{
  start.time <- Sys.time()
  # set up cluster (parallel, doParallel and later plyr)
  nCores <- parallel::detectCores()-1
  message(paste("  |Running on nCores = ", nCores, sep=""))
  # progressbar not working with plyr::alply, .progress = "text"
  message("  |...")
  # set up cluster
  cl <- parallel::makeCluster(nCores)
  doParallel::registerDoParallel(nCores)
  opts <- list(preschedule=FALSE)
  # set up all variables and functions again for Windows
  wrap <- function(x) {(x + pi) %% (2 * pi) - pi}
  
  check.extent <- function(DEM, track)
  {
    if(!class(DEM)=="RasterLayer") stop("'DEM' is not of type 'RasterLayer'")
    e <- raster::extent(DEM)
    if(!(min(track[,1]) >= e[1] && max(track[,1]) <= e[2] &&
         min(track[,2]) >= e[3] && max(track[,2]) <= e[4])) stop("The track is not inside the area of the digital elevation model.")
  }
  
  sim.cond.3d <- function(n.locs, start=c(0,0,0), end=start, a0, g0, densities, qProbs, error = FALSE, DEM = NULL, BG = NULL)
  {
    start.time <- Sys.time()
    if(!is.null(DEM)) {
      check.extent(DEM = DEM, track = data.frame(rbind(start, end)))
    }
    if(!is.null(BG)) {
      check.extent(DEM = BG, track = data.frame(rbind(start, end)))
    }
    # progress bar and time
    message(paste("  |Simulate CERW with ", n.locs, " steps", sep = ""))
    ui <- floor(n.locs/20)+1
    # replace the probability distribution for step length 1 by the one from
    # the qProbs since that one relies on more samples derived from sim
    densities[[1]] <- tail(qProbs,1)[[1]]
    # get the coordinates of the step length and turning angle bin centres
    start <- Reduce(c, start); names(start) <- c("x", "y", "z")
    end <- Reduce(c, end); names(end) <- c("x", "y", "z")
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
      atProbs <- densities$autoT(RTG[i, 6] - ts + tShift[i])
      # get the auto-difference probability for lift angle
      alProbs <- densities$autoL(RTG[i, 7] - ls + lShift[i])
      # get the auto-difference probability for step length
      adProbs <- densities$autoD(RTG[i, 8] - ds + dShift[i])
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
      a <- wrap(RTG[i, 4] + ts + tShift[i])
      # calculate the gradient
      g <- wrap(RTG[i, 5] + ls + lShift[i])
      # convert the coordinates from step length turning angle dimension
      x1 <- ((ds + dShift[i]) * sin(g) * cos(a)) + RTG[i, 1]
      y1 <- ((ds + dShift[i]) * sin(g) * sin(a)) + RTG[i, 2]
      z1 <- ((ds + dShift[i]) * cos(g)) + RTG[i, 3]
      # calculate the distances of the cell centers in the spatial domain
      # to the target (last location of the empirical track)
      endD <- as.numeric(sqrt((end[1] - x1) ^ 2 + (end[2] - y1) ^ 2 + (end[3] - z1) ^ 2))
      # calculate the azimuth of the cell centres to the target and substract from it the direction of arrival
      # resulting in turning angle towards target
      endT <- as.numeric(wrap(atan2(as.numeric(end[2] - y1), as.numeric(end[1] - x1)) - a))
      # calculate the gradient of the possibilite steps to the target and substract from it the angle of arrival
      # resulting in turning angle towards target
      endL <- as.numeric(wrap(atan2(as.numeric(sqrt((end[1] - x1) ^ 2 + (end[2] - y1) ^ 2)), as.numeric(end[3] - z1))) - g)
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
      # Apply gradient distribution or if gradientDensity is set to FALSE in get.densities.3d(), limit gradient to 0-pi.
      gProbs <- densities$gDens(g)
      gProbs[is.na(gProbs)] <- 0
      gProbs <- gProbs / sum(gProbs)
      Probs <- Probs * gProbs
      # Account for probable flight height, if a DEM is provided the relative flight height is taken
      # Otherwise only the absolute ellipsoid height.
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
        warning("Dead end encountered.")
        return(RTG)
      }else{
        # draw a point randomly based on the probability
        rP <- sample.int(nrow(densities$tldCube$values), size = 1, prob = Probs)
        # "x" "y" "z" "a" "g" "t" "l" "d" "p"
        # "1" "2" "3" "4" "5" "6" "7" "8" "9"
        RTG[i + 1, ] <- c(x1[rP], y1[rP], z1[rP], a[rP], g[rP], ts[rP], ls[rP], ds[rP], Probs[rP])
      }
    }
    # the track is forced to target location and the appropriate distance is added
    RTG[1, 8] <- NA
    RTG[n.locs,] <- c(end[1], end[2], end[3], NA, NA, NA, NA, NA, NA)
    RTG[n.locs, 8] <- sqrt((RTG[n.locs, 1] - RTG[n.locs-1, 1])^2 + 
                             (RTG[n.locs, 2] - RTG[n.locs-1, 2])^2 + 
                             (RTG[n.locs, 3] - RTG[n.locs-1, 3])^2)
    # Stop if the step length of the last step is larger than the largest possible step
    if(RTG[n.locs, 8] > max(densities$tldCube$values$step, na.rm = TRUE)*sqrt(2)) {
      RTG <- NULL
      warning("Dead end encountered in last step.")
      return(RTG)
    }
    rownames(RTG) <- c()
    colnames(RTG) <- c("x", "y", "z", "a", "g", "t", "l", "d", "p")
    return(as.data.frame(RTG))
  }
  
  # submit simulations to cluster
  simTracks <- plyr::alply(1:n.sim, function(x) {
    return(sim.cond.3d(n.locs, start, end, a0, g0, densities, qProbs, error, DEM, BG))
  }, .margins = 1, .parallel = TRUE, .paropts = list(.options.snow=opts))
  # stop cluster
  parallel::stopCluster(cl)
  return(simTracks)
}

#' Simulates multiple 'gliding & soaring' tracks with a given number of gliding steps
#'
#' Creates conditional empirical random walks in gliding mode, between a start and end point.
#' The walk is performed on a MODE layer and, if provided, additionally on a background and digital elevation layer.
#' The gliding is simulated with \link[eRTG3D]{sim.cond.3d} and soaring with \link[eRTG3D]{sim.uncond.3d},
#' therefore soaring is not restricted towards the target and can happen completly free as long as there are good thermal conditions.
#' It is important to extract for every mode in the MODE raster layer a corresponding densities object with \link[eRTG3D]{get.densities.3d} 
#' and pass them to the function.
#'
#' @param MODE raster layer containing the number/index of the mode, which should be used at each location
#' @param dGliding density object returned by the \link[eRTG3D]{get.densities.3d} function for gliding mode
#' @param dSoaring density object returned by the \link[eRTG3D]{get.densities.3d} function for soaring mode
#' @param qGliding the Q probabilites for the steps in gliding mode (\link[eRTG3D]{qProb.3d})
#' @param start numeric vector of length 3 with the coordinates of the start point
#' @param end numeric vector of length 3 with the coordinates of the end point
#' @param a0 initial incoming heading in radian
#' @param g0 initial incoming gradient/polar angle in radian
#' @param error logical: add random noise to the turn angle, lift angle and step length to account for errors measurements?
#' @param glideRatio ratio between vertical and horizontal movement, by default set to 15 meters forward movement per meter vertical movement
#' @param DEM raster layer containing a digital elevation model, covering the area between start and end point
#' @param BG a background raster layer that can be used to inform the choice of steps
#' @param smoothTransition logical: should the transitions between soaring and the following gliding sections be smoothed? Recommended to avoid dead ends 
#' @param n.sim number of simulations to produce
#' @param multicore logical: should simulations be spread to the available number of cores?
#'
#' @return A list containing 'soaring-gliding' trajectories or \code{NULL}s if dead ends have been encountered.
#' @export
#'
#' @note The MODE raster layer must be in the following structure: Gliding pixels have the value 1 and soaring pixel the values 2. \code{NA}'s are not allowed in the raster.
#'
#' @examples
#' n.sim.glidingSoaring.3d(locsVec, start = c(0,0,0), end=start, a0, g0, dList, qList, MODE)
#' @noRd
.n.sim.glidingSoaring.3d.windows <- function(MODE, dGliding, dSoaring, qGliding, start, end=end, a0, g0,
                                          error = TRUE, smoothTransition = TRUE, glideRatio, DEM = NULL, BG = NULL)
{
  start.time <- Sys.time()
  # set up cluster (parallel, doParallel and later plyr)
  nCores <- parallel::detectCores()-1
  message(paste("  |Running on nCores = ", nCores, sep=""))
  # progressbar not working with plyr::alply, .progress = "text"
  message("  |...")
  # set up cluster
  cl <- parallel::makeCluster(nCores)
  doParallel::registerDoParallel(nCores)
  opts <- list(preschedule=FALSE)
  # set up all variables and functions again for Windows
  wrap <- function(x) {(x + pi) %% (2 * pi) - pi}
  
  check.extent <- function(DEM, track)
  {
    if(!class(DEM)=="RasterLayer") stop("'DEM' is not of type 'RasterLayer'")
    e <- raster::extent(DEM)
    if(!(min(track[,1]) >= e[1] && max(track[,1]) <= e[2] &&
         min(track[,2]) >= e[3] && max(track[,2]) <= e[4])) stop("The track is not inside the area of the digital elevation model.")
  }
  
  sim.glidingSoaring.3d <- function(MODE, dGliding, dSoaring, qGliding, start=c(0,0,0), end=start, a0, g0,
                                    error = TRUE, smoothTransition = TRUE, glideRatio = 15, DEM = NULL, BG = NULL)
  {
    start.time <- Sys.time()
    # Check extent
    check.extent(DEM = MODE, track = data.frame(rbind(start, end)))
    if(!is.null(DEM)) {
      check.extent(DEM = DEM, track = data.frame(rbind(start, end)))
    }
    if(!is.null(BG)) {
      check.extent(DEM = BG, track = data.frame(rbind(start, end)))
    }
    # create dList
    dList <- list(dGliding, dSoaring)
    n.locs <- length(qGliding) + 1
    # progress bar and time
    message(paste("  |Simulate 'gliding & soaring' with ", n.locs , " gliding steps", sep = ""))
    pb <- txtProgressBar(min = 0, max = n.locs, style = 3)
    ui <- floor(n.locs/20)+1
    # replace the probability distribution for step length 1 by the one from
    # the qProbs since that one relies on more samples derived from sim
    dList[[1]][[1]] <- tail(qGliding,1)[[1]]
    # MODE: 1 = gliding, 2 = soaring, extract it.
    modeInd <- matrix(0, n.locs * 5, 2)
    m <- raster::extract(MODE, cbind(start[1], start[2]))
    # get coordinates of the tldCube
    ts <- dList[[m]]$tldCube$values$turn
    ls <- dList[[m]]$tldCube$values$lift
    ds <- dList[[m]]$tldCube$values$step
    # get probs for each combination
    tldProbs <- dList[[m]]$tldCube$values$prob
    # sample one randomly to set the initial conditions
    # for the previous to first turn and previous to first step
    # for the start point, as this is needed to inform the auto-difference
    # likelihood
    sCond <- sample(1:nrow(dList[[m]]$tldCube$values), 1, prob=tldProbs)
    # "x" "y" "z" "a" "g" "t" "l" "d" "p"
    # "1" "2" "3" "4" "5" "6" "7" "8" "9"
    RTG <- matrix(0, n.locs * 25, 9)
    RTG[1, ] <- c(start[1], start[2], start[3], a0, g0, ts[sCond], ls[sCond], ds[sCond], NA)
    # Create random noise if error is TRUE
    tShift <- lShift <- dShift <- matrix(0, n.locs * 25, 2)
    if (error) {
      for (j in 1:2){
        tShift[, j] <- runif(n.locs * 25, -dList[[j]]$tldCube$tRes / 2, dList[[j]]$tldCube$tRes / 2)
        lShift[, j] <- runif(n.locs * 25, -dList[[j]]$tldCube$lRes / 2, dList[[j]]$tldCube$lRes / 2)
        dShift[, j] <- runif(n.locs * 25, -dList[[j]]$tldCube$dRes / 2, dList[[j]]$tldCube$dRes / 2)
      }
    }
    # start creating the track step for step
    i <- 1
    deadEnd <- FALSE
    while(!deadEnd)
    {
      # If needed steps are reached, the track is stopped
      if (sum(modeInd[ ,1]) == (n.locs-1))
      {
        # update total number of steps (inluding the soaring steps!)
        n.locs <- i + 1
        RTG <- RTG[1:n.locs, ]
        # the track is forced to target location and the appropriate distance is added
        RTG[1, 8] <- NA
        RTG[n.locs,] <- c(end[1], end[2], end[3], NA, NA, NA, NA, NA, NA)
        RTG[n.locs, 8] <- sqrt((RTG[n.locs, 1] - RTG[n.locs-1, 1])^2 + 
                                 (RTG[n.locs, 2] - RTG[n.locs-1, 2])^2 + 
                                 (RTG[n.locs, 3] - RTG[n.locs-1, 3])^2)
        # Stop if the step length of the last step is larger than the largest possible step
        if(RTG[n.locs, 8] > max(dList[[m]]$tldCube$values$step, na.rm = TRUE)*3) {
          deadEnd <- TRUE
        }
        rownames(RTG) <- c()
        colnames(RTG) <- c("x", "y", "z", "a", "g", "t", "l", "d", "p")
        # close progress bar
        setTxtProgressBar(pb, sum(modeInd[ ,1])+1)
        close(pb)
        message(paste("  |Runtime: ", round(as.numeric(Sys.time()) - as.numeric(start.time), 2), " secs", sep = ""))
        return(as.data.frame(RTG))
      }
      # Exract mode
      m <- raster::extract(MODE, cbind(RTG[i, 1], RTG[i, 2]))
      if (m == 2) # Soaring
      {
        # Check glide ratio mode
        dz <- (RTG[i, 3] - end[3])
        dxy <- sqrt((RTG[i, 1] - end[1])^2 + (RTG[i, 2] - end[2])^2)
        gr <- dxy/dz
        if(gr < glideRatio && gr > 0)
        {
          # Glide ratio correction
          m <- 1
          # ... leave soaring
        } else {
          # UERW in soaring
          modeInd[i, m] <- 1
          # write on console
          cat("\r  |Mode:", m, "\r")
          flush.console()
          # get coordinates of the tldCube
          ts <- dList[[m]]$tldCube$values$turn
          ls <- dList[[m]]$tldCube$values$lift
          ds <- dList[[m]]$tldCube$values$step
          # get probs for each combination
          tldProbs <- dList[[m]]$tldCube$values$prob
          # get the auto-difference probability for turning angle
          atProbs <- dList[[m]]$autoT(RTG[i, 6] - ts + tShift[i, m])
          # get the auto-difference probability for lift angle
          alProbs <- dList[[m]]$autoL(RTG[i, 7] - ls + lShift[i, m])
          # get the auto-difference probability for step length
          adProbs <- dList[[m]]$autoD(RTG[i, 8] - ds + dShift[i, m])
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
          Probs <- (tldProbs) * (atProbs * alProbs * adProbs)^(1/3)
          # calculate the azimuth
          a <- wrap(RTG[i, 4] + ts + tShift[i, m])
          # calculate the gradient
          g <- wrap(RTG[i, 5] + ls + lShift[i, m])
          # convert the coordinates from step length turning angle dimension
          x1 <- ((ds + dShift[i, m]) * sin(g) * cos(a)) + RTG[i, 1]
          y1 <- ((ds + dShift[i, m]) * sin(g) * sin(a)) + RTG[i, 2]
          z1 <- ((ds + dShift[i, m]) * cos(g)) + RTG[i, 3]
          # Account for probable flight height, if a DEM is provided the relative flight height is taken
          # Otherwise only the absolute ellipsoid height.
          if(!is.null(DEM))
          {
            surface <- raster::extract(DEM, cbind(x1, y1))
            demP <-  dList[[m]]$hDistTopo(z1 - surface) * as.numeric(z1 >= surface)
            demP[is.na(demP)] <- 0
            demP <- demP / sum(demP)
            hProb <- dList[[m]]$hDistEllipsoid(z1)
            hProb[is.na(hProb)] <- 0
            hProb <- hProb / sum(hProb)
            Probs <- Probs * sqrt(demP * hProb)
          } else {
            hProb <- dList[[m]]$hDistEllipsoid(z1)
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
          # Apply gradient distribution or if gradientDensity is set to FALSE in get.densities.3d(), limit gradient to 0-pi.
          gProbs <- dList[[m]]$gDens(g)
          gProbs[is.na(gProbs)] <- 0
          gProbs <- gProbs / sum(gProbs)
          Probs <- Probs * gProbs
          # make sure we have no missing nor negative probabilities
          Probs[is.na(Probs)] <- 0
          Probs[Probs <= 0] <- 0
          # check whether the run might have ended up in a dead-end,
          # which will set the zero probability status to TRUE
          if(all(Probs==0)){
            deadEnd <- TRUE
          }else{
            # draw a point randomly based on the probability
            rP <- sample.int(nrow(dList[[m]]$tldCube$values), size = 1, prob = Probs)
            # "x" "y" "z" "a" "g" "t" "l" "d" "p"
            # "1" "2" "3" "4" "5" "6" "7" "8" "9"
            RTG[i + 1, ] <- c(x1[rP], y1[rP], z1[rP], a[rP], g[rP], ts[rP], ls[rP], ds[rP], Probs[rP])
            # increase counter
            i <- i + 1
            # update progress bar
            #if ((sum(modeInd[ ,1]) %% ui) == 0) {setTxtProgressBar(pb, sum(modeInd[ ,1]))}
          }
        }
      }
      if (m == 1) # Gliding
      {
        # CERW in gliding mode
        modeInd[i, m] <- 1
        # write on console
        cat("\r  |Mode:", m, "\r")
        flush.console()
        # get coordinates of the tldCube
        ts <- dList[[m]]$tldCube$values$turn
        ls <- dList[[m]]$tldCube$values$lift
        ds <- dList[[m]]$tldCube$values$step
        # get probs for each combination
        tldProbs <- dList[[m]]$tldCube$values$prob
        # get the auto-difference probability for turning angle
        atProbs <- dList[[m]]$autoT(RTG[i, 6] - ts + tShift[i, m])
        # get the auto-difference probability for lift angle
        alProbs <- dList[[m]]$autoL(RTG[i, 7] - ls + lShift[i, m])
        # get the auto-difference probability for step length
        adProbs <- dList[[m]]$autoD(RTG[i, 8] - ds + dShift[i, m])
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
        # smooths the transition between a soaring and the following gliding section
        # to avoid dead ends, but it can distort the statistical distribution of the track
        if(smoothTransition){
          # test if mode change happened (from soaring to gliding)
          if(all(colSums(rbind(c(0,0),modeInd[(i-1):i,])) == c(1,1)))
          {
            # calculate the azimuth towards target (end)
            aTarget <- wrap(atan2(end[2]-RTG[i, 2], end[1] - RTG[i, 1]))
            # calculate the azimuth towards target
            gTarget <- wrap(atan2(sqrt((end[1]-RTG[i, 1]) ^ 2 + (end[2]-RTG[i, 2]) ^ 2),
                                   (end[3]-RTG[i, 3])))
            # get difference and check for influence
            if(abs(diffAC <- aTarget-RTG[i, 4]) > pi/12) {
              # allowing from 5 to 9 steps
              timesC <- max(min(round(abs(diffAC/RTG[i, 6])), 9), 5)
              diffGC <- gTarget-RTG[i, 5]
              liftsC <- rep(diffGC/timesC, timesC) + lShift[i:(i+timesC-1), m]
              distsC <- rep(RTG[i, 8], timesC) + dShift[i:(i+timesC-1), m]
              turnsC <- rep(diffAC/timesC, timesC) + tShift[i:(i+timesC-1), m]
              aC <- cumsum(c(RTG[i, 4], turnsC))
              gC <- cumsum(c(RTG[i, 5], liftsC))
              xC <- cumsum(c(RTG[i, 1], ((distsC) * sin(gC[2:(timesC+1)]) * cos(aC[2:(timesC+1)]))))
              yC <- cumsum(c(RTG[i, 2], ((distsC) * sin(gC[2:(timesC+1)]) * sin(aC[2:(timesC+1)]))))
              zC <- cumsum(c(RTG[i, 3], ((distsC) * cos(gC[2:(timesC+1)]))))
              # "x" "y" "z" "a" "g" "t" "l" "d" "p"
              # "1" "2" "3" "4" "5" "6" "7" "8" "9"
              RTG[((i+1):(i+timesC+1)), 1:5] <- cbind(xC,yC,zC,aC,gC)
              i <- i + timesC
            }
          }
        }
        # calculate the azimuth
        a <- wrap(RTG[i, 4] + ts + tShift[i, m])
        # calculate the gradient
        g <- wrap(RTG[i, 5] + ls + lShift[i, m])
        # convert the coordinates from step length turning angle dimension
        x1 <- ((ds + dShift[i, m]) * sin(g) * cos(a)) + RTG[i, 1]
        y1 <- ((ds + dShift[i, m]) * sin(g) * sin(a)) + RTG[i, 2]
        z1 <- ((ds + dShift[i, m]) * cos(g)) + RTG[i, 3]
        # calculate the distances of the cell centers in the spatial domain
        # to the target (last location of the empirical track)
        endD <- as.numeric(sqrt((end[1] - x1) ^ 2 + (end[2] - y1) ^ 2 + (end[3] - z1) ^ 2))
        # calculate the azimuth of the cell centres to the target and substract from it the direction of arrival
        # resulting in turning angle towards target
        endT <- as.numeric(wrap(atan2(as.numeric(end[2] - y1), as.numeric(end[1] - x1)) - a))
        # calculate the gradient of the possibilite steps to the target and substract from it the angle of arrival
        # resulting in lift angle towards target
        endL <- as.numeric(wrap(atan2(as.numeric(sqrt((end[1] - x1) ^ 2 + (end[2] - y1) ^ 2)), as.numeric(end[3] - z1))) - g)
        # get the probabilities of making it distance and turning angle wise
        # which is derived from the two dimensional probability distribution for the
        # appropriate step being modelled
        # get possible coordinates
        qCube <- qGliding[[sum(modeInd[, m])]]
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
        # Apply gradient distribution or if gradientDensity is set to FALSE in get.densities.3d(), limit gradient to 0-pi.
        gProbs <- dList[[m]]$gDens(g)
        gProbs[is.na(gProbs)] <- 0
        gProbs <- gProbs / sum(gProbs)
        Probs <- Probs * gProbs
        # Account for probable flight height, if a DEM is provided the relative flight height is taken
        # Otherwise only the absolute ellipsoid height.
        if(!is.null(DEM))
        {
          surface <- raster::extract(DEM, cbind(x1, y1))
          demP <-  dList[[m]]$hDistTopo(z1 - surface) * as.numeric(z1 >= surface)
          demP[is.na(demP)] <- 0
          demP <- demP / sum(demP)
          hProb <- dList[[m]]$hDistEllipsoid(z1)
          hProb[is.na(hProb)] <- 0
          hProb <- hProb / sum(hProb)
          Probs <- Probs * sqrt(demP * hProb)
        } else {
          hProb <- dList[[m]]$hDistEllipsoid(z1)
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
          deadEnd <- TRUE
        }else{
          # draw a point randomly based on the probability
          rP <- sample.int(nrow(dList[[m]]$tldCube$values), size = 1, prob = Probs)
          # "x" "y" "z" "a" "g" "t" "l" "d" "p"
          # "1" "2" "3" "4" "5" "6" "7" "8" "9"
          RTG[i + 1, ] <- c(x1[rP], y1[rP], z1[rP], a[rP], g[rP], ts[rP], ls[rP], ds[rP], Probs[rP])
          # increase counter
          i <- i + 1
          # update progress bar
          if ((sum(modeInd[ ,1]) %% ui) == 0) {setTxtProgressBar(pb, sum(modeInd[ ,1]))}
        }
      }
    }
    warning("Dead end encountered.")
    return(RTG <- NULL)
  }
  
  # submit simulations to cluster
  simTracks <- plyr::alply(1:n.sim, function(x) {
    return(sim.glidingSoaring.3d(MODE = MODE, dGliding = dGliding, dSoaring = dSoaring, qGliding = qGliding, start=start, end=end, a0=a0, g0=g0,
                            error = error, smoothTransition = smoothTransition, glideRatio = glideRatio, DEM = DEM, BG = BG))
  }, .margins = 1, .parallel = TRUE, .paropts = list(.options.snow=opts))
  # stop cluster
  parallel::stopCluster(cl)
  return(simTracks)
}
