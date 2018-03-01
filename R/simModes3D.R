#' Multiple Conditioned Empirical Random Walks (CERW) with modes in 3D
#'
#' Creates n conditioned empirical random walks using different modes, with a specific starting and ending point,
#' geometrically similar to the initial trajectory by applying \link[eRTG3D]{sim.cond.modes.3d} multiple times.
#'
#' @param n.sim number of CERWs to simulate
#' @param locsVec length of the trajectory in locations
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
#' n.sim.cond.modes.3d(n.sim, locsVec, start = c(0,0,0), end=start, a0, g0, dList, qList, MODE)
n.sim.cond.modes.3d <- function(n.sim, locsVec, start = c(0,0,0), end = start, a0, g0, dList, qList, glideRatio = 20,
                                error = FALSE, multicore = FALSE, MODE, DEM = NULL, BG = NULL)
{
  start.time <- Sys.time()
  n.sim <- round(n.sim)
  if (n.sim <= 1) {return(sim.cond.modes.3d(locsVec, start = start, end = end, MODE = MODE, a0 = a0, g0 = g0,
                                           dList = dList, qList = qList, glideRatio = glideRatio, error = error, DEM = DEM, BG = BG))}
  message(paste("  |Simulate ", n.sim ," CERWs with ca. ", round(mean(locsVec)), " steps", sep = ""))
  if (multicore) {
    if(.Platform$OS.type == "unix") {
      nCores <- parallel::detectCores()-1
      message(paste("  |Running on nCores = ", nCores, sep=""))
      message("  |...")
      cerwList <- parallel::mclapply(X = 1:n.sim, FUN = function(X) {
        sim.cond.modes.3d(locsVec, start = start, end = end, MODE = MODE, a0 = a0, g0 = g0,
                          dList = dList, qList = qList, glideRatio = glideRatio, error = error, DEM = DEM, BG = BG)},
        mc.cores = nCores)
    }
    if(.Platform$OS.type == "windows") {
      stop("Multicore not yet implemented on Windows system, please use use a unix based system.")
    }
  } else {
    cerwList <- suppressMessages(lapply(X = 1:n.sim, FUN = function(X) {
      sim.cond.modes.3d(locsVec, start=start, end=end, MODE = MODE, a0 = a0, g0 = g0,
                        dList=dList, qList=qList, error = error, DEM = DEM, BG = BG)}))
  }
  message(paste("  |Runtime: ", round(as.numeric(Sys.time()) - as.numeric(start.time), 2), " secs", sep = ""))
  return(cerwList)
}

#' Conditioned Empirical Random Walk (CERW) in 3D with different movement modes
#'
#' Creates a conditioned empirical random walk, with a specific starting and ending point,
#' geometrically similar to the initial trajectory, by using different modes of movement.
#' It is important to extract for every mode in the MODE raster layer a corresponding D object \link[eRTG3D]{get.densities.3d}
#' and the Q probabilities \link[eRTG3D]{qProb.3d}. The number of steps has to be the same for every mode that is used.
#'
#' @param locsVec vector of maximum lengths of modes
#' @param start numeric vector of length 3 with the coordinates of the start point
#' @param end numeric vector of length 3 with the coordinates of the end point
#' @param a0 initial incoming heading in radian
#' @param g0 initial incoming gradient/polar angle in radian
#' @param dList list of list objects returned by the \link[eRTG3D]{get.densities.3d} function, one for each MODE
#' @param qList list of lists with objects returned by the \link[eRTG3D]{qProb.3d} function, one for each MODE
#' @param error logical: add random noise to the turn angle, lift angle and step length to account for errors measurements?
#' @param MODE raster layer containing the number/index of the mode, which should be used at each location
#' @param DEM raster layer containing a digital elevation model, covering the area between start and end point
#' @param BG a background raster layer that can be used to inform the choice of steps
#'
#' @return A trajectory in the form of data.frame
#' @export
#'
#' @examples
#' sim.cond.modes.3d(locsVec, start = c(0,0,0), end=start, a0, g0, dList, qList, MODE)
sim.cond.modes.3d <- function(locsVec, start=c(0,0,0), end=start, a0, g0, dList, qList, error = FALSE, glideRatio = 20, MODE, DEM = NULL, BG = NULL)
{
  start.time <- Sys.time()
  n.locs <- max(locsVec)
  .check.extent(DEM = MODE, track = data.frame(rbind(start, end)))
  if(!is.null(DEM)) {
    .check.extent(DEM = DEM, track = data.frame(rbind(start, end)))
  }
  if(!is.null(BG)) {
    .check.extent(DEM = BG, track = data.frame(rbind(start, end)))
  }
  # progress bar and time
  message(paste("  |Simulate CERW with ca. ", round(mean(locsVec)) , " steps", sep = ""))
  pb <- txtProgressBar(min = 0, max = n.locs-2, style = 3)
  ui <- floor(n.locs/20)+1
  # replace the probability distribution for step length 1 by the one from
  # the qProbs since that one relies on more samples derived from sim
  dList[[1]][[1]] <- tail(qList[[1]],1)[[1]]
  dList[[2]][[1]] <- tail(qList[[2]],1)[[1]]
  # get the coordinates of the step length and turning angle bin centres
  names(start) <- c("x", "y", "z")
  names(end) <- c("x", "y", "z")
  # Exract mode
  modes <- unique(raster::values(MODE))
  m <- raster::extract(MODE, cbind(start[1], start[2]))
  densities <- dList[[m]]
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
  tShift <- lShift <- dShift <- matrix(0, n.locs-2, length(modes))
  if (error) {
    for (j in 1:length(modes)){
      tShift[, m] <- runif(n.locs - 2, -dList[[m]]$tldCube$tRes / 2, dList[[m]]$tldCube$tRes / 2)
      lShift[, m] <- runif(n.locs - 2, -dList[[m]]$tldCube$lRes / 2, dList[[m]]$tldCube$lRes / 2)
      dShift[, m] <- runif(n.locs - 2, -dList[[m]]$tldCube$dRes / 2, dList[[m]]$tldCube$dRes / 2)
    }
  }
  # start creating the track step for step
  i <- 1
  modeInd <- matrix(0, n.locs, 2)
  for (i in 1:(n.locs - 2))
  {
    # Exract mode
    m <- raster::extract(MODE, cbind(RTG[i, 1], RTG[i, 2]))
    # Check glide ratio
    dz <- (RTG[i, 3] - end[3])
    dxy <- sqrt((RTG[i, 1] - end[1])^2 + (RTG[i, 2] - end[2])^2)
    gr <- dxy/dz
    if(gr < glideRatio && gr > 0)
    {
      # Glide ratio correction
      m <- 1
    }
    modeInd[i, m] <- 1
    cat("\r  |Mode:", m, "\r")
    flush.console()
    densities <- dList[[m]]
    qProbs <- qList[[m]]
    ###
    # get coordinates of the tldCube
    ts <- densities$tldCube$values$turn
    ls <- densities$tldCube$values$lift
    ds <- densities$tldCube$values$step
    # get probs for each combination
    tldProbs <- densities$tldCube$values$prob
    ###
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
    a <- .wrap(RTG[i, 4] + ts + tShift[i, m])
    # calculate the gradient
    g <- .wrap(RTG[i, 5] + ls + lShift[i, m])
    # convert the coordinates from step length turning angle dimension
    x1 <- ((ds + dShift[i, m]) * sin(g) * cos(a)) + RTG[i, 1]
    y1 <- ((ds + dShift[i, m]) * sin(g) * sin(a)) + RTG[i, 2]
    z1 <- ((ds + dShift[i, m]) * cos(g)) + RTG[i, 3]
    # calculate the distances of the cell centers in the spatial domain
    # to the target (last location of the empirical track)
    endD <- as.numeric(sqrt((end[1] - x1) ^ 2 + (end[2] - y1) ^ 2 + (end[3] - z1) ^ 2))
    # calculate the azimuth of the cell centres to the target and substract from it the direction of arrival
    # resulting in turning angle towards target
    endT <- as.numeric(.wrap(atan2(as.numeric(end[2] - y1), as.numeric(end[1] - x1)) - a))
    # calculate the gradient of the possibilite steps to the target and substract from it the angle of arrival
    # resulting in lift angle towards target
    endL <- as.numeric(.wrap(atan2(as.numeric(sqrt((end[1] - x1) ^ 2 + (end[2] - y1) ^ 2)), as.numeric(end[3] - z1))) - g)
    # get the probabilities of making it distance and turning angle wise
    # which is derived from the two dimensional probability distribution for the
    # appropriate step being modelled
    # get possible coordinates
    qCube <- qProbs[[sum(modeInd[, m])]]
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
    
    # # limit gradient to p/2 and pi for gliding (m=1) 0 and pi/2 for soaring (m=2)
    # if(m == 1){
    #   #Probs <- Probs * as.numeric((g > 0) & (g < (pi))) # gliding
    #   Probs <- Probs * as.numeric((g > (pi/2)) & (g < pi)) # gliding
    # } else {
    #   Probs <- Probs * as.numeric((g > 0) & (g < (pi/2))) # soaring
    # }
    
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
      if(i == 1) {
        RTG <- NULL
        warning("Dead end encountered in first step.")
        return(RTG)
      }
      RTG <- RTG[1:i, ]
      rownames(RTG) <- c()
      colnames(RTG) <- c("x", "y", "z", "a", "g", "t", "l", "d", "p")
      RTG <- as.data.frame(RTG)
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
    # If needed steps are reached, the track is stopped
    if (any(colSums(modeInd, dims = 1) == (locsVec-2)))
    {
      n.locs <- i+2
      RTG <- RTG[1:n.locs, ]
      # the track is forced to target location and the appropriate distance is added
      RTG[1, 8] <- NA
      RTG[n.locs,] <- c(end[1], end[2], end[3], NA, NA, NA, NA, NA, NA)
      RTG[n.locs, 8] <- sqrt((RTG[n.locs, 1] - RTG[n.locs-1, 1])^2 + 
                               (RTG[n.locs, 2] - RTG[n.locs-1, 2])^2 + 
                               (RTG[n.locs, 3] - RTG[n.locs-1, 3])^2)
      # Stop if the step length of the last step is larger than the largest possible step
      if(RTG[n.locs, 8] > max(densities$tldCube$values$step, na.rm = TRUE)*3) {
        RTG <- RTG[1:n.locs-1, ]
        rownames(RTG) <- c()
        colnames(RTG) <- c("x", "y", "z", "a", "g", "t", "l", "d", "p")
        RTG <- as.data.frame(RTG)
        close(pb)
        message(paste("  |Runtime: ", round(as.numeric(Sys.time()) - as.numeric(start.time), 2), " secs", sep = ""))
        warning("Dead end encountered in last step.")
        return(RTG)
      }
      rownames(RTG) <- c()
      colnames(RTG) <- c("x", "y", "z", "a", "g", "t", "l", "d", "p")
      # close progress bar
      setTxtProgressBar(pb, max(locsVec)-2)
      close(pb)
      message(paste("  |Runtime: ", round(as.numeric(Sys.time()) - as.numeric(start.time), 2), " secs", sep = ""))
      return(as.data.frame(RTG))
    }
  }
}