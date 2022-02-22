#' Calculate glide ratio
#'
#' Calculates the ratio between horizontal movement and vertical movement.
#' The value expresses the distance covered forward movement per distance movement in sinking.
#'
#' @param track a track data.frame containing x, y and z coordinates of a gliding section
#'
#' @return The ratio between horizontal and vertical movement.
#' @export
#'
#' @examples
#' get.glideRatio.3d(niclas)
get.glideRatio.3d <- function(track) {
  start <- track[2, 1:3]
  end <- track[nrow(track), 1:3]
  return(-1 * (sqrt((end[1] - start[1])^2 + (end[2] - start[2])^2) / as.numeric(end[3] - start[3])))
}

#' Simulates 'gliding & soaring' track with a given number of gliding steps
#'
#' Creates a conditional empirical random walk in gliding mode, between a start and end point.
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
#' @param verbose logical: print current mode used?
#'
#' @return A 'soaring-gliding' trajectory in the form of data.frame
#' @export
#'
#' @note The MODE raster layer must be in the following structure: Gliding pixels have the value 1 and soaring pixel the values 2. \code{NA}'s are not allowed in the raster.
#'
#' @examples
#' print("tbd.")
sim.glidingSoaring.3d <- function(MODE, dGliding, dSoaring, qGliding, start = c(0, 0, 0), end = start, a0, g0,
                                  error = TRUE, smoothTransition = TRUE, glideRatio = 15, DEM = NULL, BG = NULL,
                                  verbose = FALSE) {
  start.time <- Sys.time()
  # Check extent
  .check.extent(DEM = MODE, track = data.frame(rbind(start, end)))
  if (!is.null(DEM)) {
    .check.extent(DEM = DEM, track = data.frame(rbind(start, end)))
  }
  if (!is.null(BG)) {
    .check.extent(DEM = BG, track = data.frame(rbind(start, end)))
  }
  # create dList
  dList <- list(dGliding, dSoaring)
  n.locs <- length(qGliding) + 1
  # progress bar and time
  message(paste("  |Simulate 'gliding & soaring' with ", n.locs, " gliding steps", sep = ""))
  pb <- utils::txtProgressBar(min = 0, max = n.locs, style = 3)
  ui <- floor(n.locs / 20) + 1
  # replace the probability distribution for step length 1 by the one from
  # the qProbs since that one relies on more samples derived from sim
  dList[[1]][[1]] <- utils::tail(qGliding, 1)[[1]]
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
  sCond <- sample(seq_len(nrow(dList[[m]]$tldCube$values)), 1, prob = tldProbs)
  # "x" "y" "z" "a" "g" "t" "l" "d" "p"
  # "1" "2" "3" "4" "5" "6" "7" "8" "9"
  RTG <- matrix(0, n.locs * 25, 9)
  RTG[1, ] <- c(start[1], start[2], start[3], a0, g0, ts[sCond], ls[sCond], ds[sCond], NA)
  # Create random noise if error is TRUE
  tShift <- lShift <- dShift <- matrix(0, n.locs * 25, 2)
  if (error) {
    for (j in 1:2) {
      tShift[, j] <- stats::runif(n.locs * 25, -dList[[j]]$tldCube$tRes / 2, dList[[j]]$tldCube$tRes / 2)
      lShift[, j] <- stats::runif(n.locs * 25, -dList[[j]]$tldCube$lRes / 2, dList[[j]]$tldCube$lRes / 2)
      dShift[, j] <- stats::runif(n.locs * 25, -dList[[j]]$tldCube$dRes / 2, dList[[j]]$tldCube$dRes / 2)
    }
  }
  # start creating the track step for step
  i <- 1
  deadEnd <- FALSE
  while (!deadEnd) {
    # If needed steps are reached, the track is stopped
    if (sum(modeInd[, 1]) == (n.locs - 1)) {
      # update total number of steps (inluding the soaring steps!)
      n.locs <- i + 1
      RTG <- RTG[1:n.locs, ]
      # the track is forced to target location and the appropriate distance is added
      RTG[1, 8] <- NA
      RTG[n.locs, ] <- c(end[1], end[2], end[3], NA, NA, NA, NA, NA, NA)
      RTG[n.locs, 8] <- sqrt((RTG[n.locs, 1] - RTG[n.locs - 1, 1])^2 +
        (RTG[n.locs, 2] - RTG[n.locs - 1, 2])^2 +
        (RTG[n.locs, 3] - RTG[n.locs - 1, 3])^2)
      # Stop if the step length of the last step is larger than the largest possible step
      if (RTG[n.locs, 8] > max(dList[[m]]$tldCube$values$step, na.rm = TRUE) * 3) {
        deadEnd <- TRUE
      }
      rownames(RTG) <- c()
      colnames(RTG) <- c("x", "y", "z", "a", "g", "t", "l", "d", "p")
      # close progress bar
      utils::setTxtProgressBar(pb, sum(modeInd[, 1]) + 1)
      close(pb)
      message(paste("  |Runtime: ", round(as.numeric(Sys.time()) - as.numeric(start.time), 2), " secs", sep = ""))
      return(as.data.frame(RTG))
    }
    # Exract mode
    m <- raster::extract(MODE, cbind(RTG[i, 1], RTG[i, 2]))
    if (m == 2) { # Soaring
      # Check glide ratio mode
      dz <- (RTG[i, 3] - end[3])
      dxy <- sqrt((RTG[i, 1] - end[1])^2 + (RTG[i, 2] - end[2])^2)
      gr <- dxy / dz
      if (gr < glideRatio && gr > 0) {
        # Glide ratio correction
        m <- 1
        # ... leave soaring
      } else {
        # UERW in soaring
        modeInd[i, m] <- 1
        # write on console
        if (verbose) {
          cat("\r  |Mode:", m, "\r")
          utils::flush.console()
        }
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
        Probs <- (tldProbs) * (atProbs * alProbs * adProbs)^(1 / 3)
        # calculate the azimuth
        a <- .wrap(RTG[i, 4] + ts + tShift[i, m])
        # calculate the gradient
        g <- .wrap(RTG[i, 5] + ls + lShift[i, m])
        # convert the coordinates from step length turning angle dimension
        x1 <- ((ds + dShift[i, m]) * sin(g) * cos(a)) + RTG[i, 1]
        y1 <- ((ds + dShift[i, m]) * sin(g) * sin(a)) + RTG[i, 2]
        z1 <- ((ds + dShift[i, m]) * cos(g)) + RTG[i, 3]
        # Account for probable flight height, if a DEM is provided the relative flight height is taken
        # Otherwise only the absolute ellipsoid height.
        if (!is.null(DEM)) {
          surface <- raster::extract(DEM, cbind(x1, y1))
          demP <- dList[[m]]$hDistTopo(z1 - surface) * as.numeric(z1 >= surface)
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
        if (!is.null(BG)) {
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
        if (all(Probs == 0)) {
          deadEnd <- TRUE
        } else {
          # draw a point randomly based on the probability
          rP <- sample.int(nrow(dList[[m]]$tldCube$values), size = 1, prob = Probs)
          # "x" "y" "z" "a" "g" "t" "l" "d" "p"
          # "1" "2" "3" "4" "5" "6" "7" "8" "9"
          RTG[i + 1, ] <- c(x1[rP], y1[rP], z1[rP], a[rP], g[rP], ts[rP], ls[rP], ds[rP], Probs[rP])
          # increase counter
          i <- i + 1
          # update progress bar
          # if ((sum(modeInd[ ,1]) %% ui) == 0) {utils::setTxtProgressBar(pb, sum(modeInd[ ,1]))}
        }
      }
    }
    if (m == 1) # Gliding
      {
        # CERW in gliding mode
        modeInd[i, m] <- 1
        # write on console
        if (verbose) {
          cat("\r  |Mode:", m, "\r")
          utils::flush.console()
        }
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
        P <- (tldProbs) * (atProbs * alProbs * adProbs)^(1 / 3)
        # smooths the transition between a soaring and the following gliding section
        # to avoid dead ends, but it can distort the statistical distribution of the track
        if (smoothTransition) {
          # test if mode change happened (from soaring to gliding)
          if (all(colSums(rbind(c(0, 0), modeInd[(i - 1):i, ])) == c(1, 1))) {
            # calculate the azimuth towards target (end)
            aTarget <- .wrap(atan2(end[2] - RTG[i, 2], end[1] - RTG[i, 1]))
            # calculate the azimuth towards target
            gTarget <- .wrap(atan2(
              sqrt((end[1] - RTG[i, 1])^2 + (end[2] - RTG[i, 2])^2),
              (end[3] - RTG[i, 3])
            ))
            # get difference and check for influence
            if (abs(diffAC <- aTarget - RTG[i, 4]) > pi / 12) {
              # allowing from 5 to 9 steps
              timesC <- max(min(round(abs(diffAC / RTG[i, 6])), 9), 5)
              diffGC <- gTarget - RTG[i, 5]
              liftsC <- rep(diffGC / timesC, timesC) + lShift[i:(i + timesC - 1), m]
              distsC <- rep(RTG[i, 8], timesC) + dShift[i:(i + timesC - 1), m]
              turnsC <- rep(diffAC / timesC, timesC) + tShift[i:(i + timesC - 1), m]
              aC <- cumsum(c(RTG[i, 4], turnsC))
              gC <- cumsum(c(RTG[i, 5], liftsC))
              xC <- cumsum(c(RTG[i, 1], ((distsC) * sin(gC[2:(timesC + 1)]) * cos(aC[2:(timesC + 1)]))))
              yC <- cumsum(c(RTG[i, 2], ((distsC) * sin(gC[2:(timesC + 1)]) * sin(aC[2:(timesC + 1)]))))
              zC <- cumsum(c(RTG[i, 3], ((distsC) * cos(gC[2:(timesC + 1)]))))
              # "x" "y" "z" "a" "g" "t" "l" "d" "p"
              # "1" "2" "3" "4" "5" "6" "7" "8" "9"
              RTG[((i + 1):(i + timesC + 1)), 1:5] <- cbind(xC, yC, zC, aC, gC)
              i <- i + timesC
            }
          }
        }
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
        endD <- as.numeric(sqrt((end[1] - x1)^2 + (end[2] - y1)^2 + (end[3] - z1)^2))
        # calculate the azimuth of the cell centres to the target and substract from it the direction of arrival
        # resulting in turning angle towards target
        endT <- as.numeric(.wrap(atan2(as.numeric(end[2] - y1), as.numeric(end[1] - x1)) - a))
        # calculate the gradient of the possibilite steps to the target and substract from it the angle of arrival
        # resulting in lift angle towards target
        endL <- as.numeric(.wrap(atan2(as.numeric(sqrt((end[1] - x1)^2 + (end[2] - y1)^2)), as.numeric(end[3] - z1))) - g)
        # get the probabilities of making it distance and turning angle wise
        # which is derived from the two dimensional probability distribution for the
        # appropriate step being modelled
        # get possible coordinates
        qCube <- qGliding[[sum(modeInd[, m])]]
        tVal <- unique(qCube$values$turn)
        lVal <- unique(qCube$values$lift)
        dVal <- unique(qCube$values$step)
        # find closest coordinates
        tCoords <- unlist(lapply(endT, function(x) tVal[which.min(abs(tVal - x))]))
        lCoords <- unlist(lapply(endL, function(x) lVal[which.min(abs(lVal - x))]))
        dCoords <- unlist(lapply(endD, function(x) dVal[which.min(abs(dVal - x))]))
        # extract Q
        Q <- unlist(lapply(seq_len(length(tCoords)), function(x) {
          test <- (qCube$values$turn == tCoords[x] & qCube$values$lift == lCoords[x] & qCube$values$step == dCoords[x])
          if (any(test == TRUE)) {
            return(qCube$values$prob[test])
          } else {
            return(0)
          }
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
        if (!is.null(DEM)) {
          surface <- raster::extract(DEM, cbind(x1, y1))
          demP <- dList[[m]]$hDistTopo(z1 - surface) * as.numeric(z1 >= surface)
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
        if (!is.null(BG)) {
          bgP <- raster::extract(BG, cbind(x1, y1))
          Probs <- Probs * bgP
        }
        # make sure we have no missing nor negative probabilities
        Probs[is.na(Probs)] <- 0
        Probs[Probs <= 0] <- 0
        # check whether the run might have ended up in a dead-end,
        # which will set the zero probability status to TRUE
        if (all(Probs == 0)) {
          deadEnd <- TRUE
        } else {
          # draw a point randomly based on the probability
          rP <- sample.int(nrow(dList[[m]]$tldCube$values), size = 1, prob = Probs)
          # "x" "y" "z" "a" "g" "t" "l" "d" "p"
          # "1" "2" "3" "4" "5" "6" "7" "8" "9"
          RTG[i + 1, ] <- c(x1[rP], y1[rP], z1[rP], a[rP], g[rP], ts[rP], ls[rP], ds[rP], Probs[rP])
          # increase counter
          i <- i + 1
          # update progress bar
          if ((sum(modeInd[, 1]) %% ui) == 0) {
            utils::setTxtProgressBar(pb, sum(modeInd[, 1]))
          }
        }
      }
  }
  warning("Dead end encountered.")
  # RTG <- as.data.frame(RTG[1:i-1,])
  # colnames(RTG) <- c("x", "y", "z", "a", "g", "t", "l", "d", "p")
  # return(as.data.frame(RTG[1:i-1,]))
  return(RTG <- NULL)
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
#' @param verbose logical: print current mode used?
#' @param n.sim number of simulations to produce
#' @param parallel logical: run computations in parallel (n-1 cores)? Or numeric: the number of nodes (maximum: n - 1 cores)
#'
#' @return A list containing 'soaring-gliding' trajectories or \code{NULL}s if dead ends have been encountered.
#' @export
#'
#' @note The MODE raster layer must be in the following structure: Gliding pixels have the value 1 and soaring pixel the values 2. \code{NA}'s are not allowed in the raster.
#'
#' @examples
#' print("tbd.")
n.sim.glidingSoaring.3d <- function(n.sim = 1, parallel = FALSE, MODE, dGliding, dSoaring, qGliding, start = c(0, 0, 0), end = start, a0, g0,
                                    error = TRUE, smoothTransition = TRUE, glideRatio = 20, DEM = NULL, BG = NULL, verbose = FALSE) {
  n.sim <- round(n.sim)
  if (n.sim <= 1) {
    return(sim.glidingSoaring.3d(
      MODE = MODE, dGliding = dGliding, dSoaring = dSoaring, qGliding = qGliding, start = start, end = end, a0 = a0, g0 = g0,
      error = error, smoothTransition = smoothTransition, glideRatio = glideRatio, DEM = DEM, BG = BG, verbose = verbose
    ))
  }
  nNodes <- .nNodes(parallel)
  message(paste("  |Simulate ", n.sim, " 'gliding & soaring' with ", (length(qGliding) + 1), " gliding steps", sep = ""))
  cerwList <- parpblapply(
    X = 1:n.sim, FUN = function(x) {
      n.sim.glidingSoaring.3d(
        MODE = MODE, dGliding = dGliding, dSoaring = dSoaring,
        qGliding = qGliding, start = start, end = end, a0 = a0, g0 = g0,
        error = error, smoothTransition = smoothTransition,
        glideRatio = glideRatio, DEM = DEM, BG = BG, verbose = verbose
      )
    },
    export = c(
      "MODE", "dGliding", "dSoaring", "qGliding",
      "start", "end", "a0", "g0", "error",
      "smoothTransition", "glideRatio", "DEM", "BG", "verbose"
    ),
    packages = c("eRTG3D"),
    nNodes = nNodes, envir = environment()
  )
  return(cerwList)
}
