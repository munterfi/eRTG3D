#' Extract tldCube and autodifference approximation functions
#'
#' Creates a list consisting of the three dimensional
#' probability distribution cube for turning angle, lift angle and step length (\link[eRTG3D]{turnLiftStepHist})
#' as well as the uni-dimensional distributions of the differences
#' of the turn angles, lift angles and step lengths with a lag of 1 to maintain
#' minimal level of autocorrelation in each of the terms.
#' Additionally also the distribution of the flight height over the ellipsoid (absolute)
#' and the distribution of flight height over the topography (relative) can be included.
#'
#' @param turnAngle turn angles of the track (t)
#' @param liftAngle lift angles of the track (l)
#' @param stepLength stepLength of the track (d)
#' @param deltaLift auto differences of the turn angles (diff(t))
#' @param deltaTurn auto differences of the lift angles (diff(l))
#' @param deltaStep auto differences of the step length (diff(d))
#' @param gradientAngle \code{NULL} or the gardient angles of the track
#' @param heightEllipsoid flight height over the ellipsoid (absolute) or \code{NULL} to exclude this distribution
#' @param heightTopo flight height over the topography (relative) or \code{NULL} to exclude this distribution
#' @param maxBin numeric scalar, maximum number of bins per dimension of the tld-cube (\link[eRTG3D]{turnLiftStepHist})
#'
#' @return A list containing the tldCube and the autodifferences functions (and additionally the flight height distribution functions)
#' @export
#'
#' @examples
#' niclas <- track.properties.3d(niclas)[2:nrow(niclas), ]
#' P <- get.densities.3d(
#'   turnAngle = niclas$t, liftAngle = niclas$l, stepLength = niclas$d,
#'   deltaLift = diff(niclas$t), deltaTurn = diff(niclas$l), deltaStep = diff(niclas$d),
#'   gradientAngle = NULL, heightEllipsoid = NULL, heightTopo = NULL, maxBin = 25
#' )
get.densities.3d <- function(turnAngle, liftAngle, stepLength, deltaLift, deltaTurn, deltaStep, gradientAngle = NULL, heightEllipsoid = NULL, heightTopo = NULL, maxBin = 25) {
  # probability distribution cube for turning angle, lift angle and step length
  cubeTLD <- turnLiftStepHist(turn = turnAngle, lift = liftAngle, step = stepLength, maxBin = maxBin)
  # approximate the distribution of the difference in turning angle with lag 1
  autoT <- stats::approxfun(stats::density.default(deltaTurn))
  # approximate the distribution of the difference in lift angle with lag 1
  autoL <- stats::approxfun(stats::density.default(deltaLift))
  # approximate the distribution of the difference in step length with lag 1
  autoD <- stats::approxfun(stats::density.default(deltaStep))
  if (!is.null(gradientAngle)) {
    gDens <- stats::approxfun(stats::density.default(gradientAngle[gradientAngle > 0 & gradientAngle < pi]))
  } else {
    gDens <- function(x) {
      return(as.numeric(x > 0 & x < pi))
    }
  }
  if (!is.null(heightEllipsoid)) {
    hDistEllipsoid <- stats::approxfun(stats::density.default(heightEllipsoid))
  } else {
    hDistEllipsoid <- function(x) {
      1
    }
  }
  if (!is.null(heightTopo)) {
    hDistTopo <- stats::approxfun(stats::density.default(heightTopo))
  } else {
    hDistTopo <- function(x) {
      1
    }
  }
  return(list(tldCube = cubeTLD, autoT = autoT, autoL = autoL, autoD = autoD, gDens = gDens, hDistEllipsoid = hDistEllipsoid, hDistTopo = hDistTopo))
}

#' Three dimensional histogram
#'
#' Derives a three dimensional distribution of a turn angle,
#' lift angle and step length, using the Freedman–Diaconis rule for
#' estimating the number of bins.
#'
#' @param turn numeric vector of turn angles
#' @param lift numeric vector of lift angles
#' @param step numeric vector of step lengths
#' @param printDims logical: should dimensions of tld-Cube be messaged?
#' @param rm.zeros logical: should combinations with zero probability be removed?
#' @param maxBin numeric scalar, maximum number of bins per dimension of the tld-cube.
#'
#' @return A three dimensional histogram as data.frame
#' @export
#'
#' @examples
#' niclas <- track.properties.3d(niclas)[2:nrow(niclas), ]
#' turnLiftStepHist(niclas$t, niclas$l, niclas$d)
turnLiftStepHist <- function(turn, lift, step, printDims = TRUE, rm.zeros = TRUE, maxBin = 25) {
  # define number of bins based on Freedman-Diaconis
  nx <- min(grDevices::nclass.FD(turn), maxBin)
  ny <- min(grDevices::nclass.FD(lift), maxBin)
  nz <- min(grDevices::nclass.FD(step), maxBin)
  if (printDims) {
    message("  |TLD cube dimensions: ", nx, " x ", ny, " x ", nz)
  }
  # create histogram
  tCuts <- .cutMidpoints(turn, nx)
  lCuts <- .cutMidpoints(lift, ny)
  dCuts <- .cutMidpoints(step, nz)
  h <- list(
    turn = tCuts[[1]],
    lift = lCuts[[1]],
    step = dCuts[[1]]
  )
  h <- do.call(data.frame, h)
  h <- as.data.frame(table(h))
  # resolutions
  tRes <- tCuts[[2]]
  lRes <- lCuts[[2]]
  dRes <- dCuts[[2]]
  # probabilities
  colnames(h)[4] <- "prob"
  # Remove zeros if desired
  if (rm.zeros) {
    h <- h[!h$prob == 0, ]
  }
  # normalize frequency to get probabilities
  h$prob <- h$prob / sum(h$prob)
  # convert factors to numeric
  h[1:3] <- lapply(h[1:3], function(x) {
    as.numeric(levels(x))[x]
  })
  return(list(values = h, tRes = tRes, lRes = lRes, dRes = dRes))
}

#' Find corresponding midpoints of breaks
#'
#' Assigns each element of a input vector the mid point of the corresponding interval.
#' Adjusted the original \code{cut()} function from the base package.
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
#' @noRd
.cutMidpoints <- function(x, breaks, rm.empty = TRUE) {
  if (breaks <= 1) {
    return(list(cuts = factor(rep(mean(x), length(x))), res = 0))
  }
  nb <- as.integer(breaks + 1)
  dx <- diff(rx <- range(x, na.rm = TRUE))
  breaks <- seq.int(rx[1L], rx[2L], length.out = nb)
  res <- stats::median(diff(breaks))
  breaks[c(1L, nb)] <- c(rx[1L] - dx / 1000, rx[2L] +
    dx / 1000)
  code <- .bincode(x, breaks, right = TRUE, include.lowest = FALSE)
  width <- diff(breaks)
  minBreak <- min(breaks)
  midpoints <- sapply(1:(length(breaks) - 1), function(ii) {
    minBreak + sum(width[1:ii - 1]) + width[ii] / 2
  })
  if (rm.empty) {
    list(cuts = factor(midpoints[code]), res = res)
  } else {
    list(cuts = factor(midpoints[code], midpoints), res = res)
  }
}

#' Unconditional Empirical Random Walk (UERW) in 3-D
#'
#' This function creates unconditional walks with prescribed
#' empirical properties (turning angle, lift angle and step length
#' and the auto-differences of them. It can be used for uncon-
#' ditional walks or to seed the conditional walks with
#' comparably long simulations.
#' The conditional walk connecting a given start
#' with a certain end point by a given number of
#' steps needs an attraction term (the Q probability, see \link[eRTG3D]{qProb.3d})
#' to ensure that the target is approached and hit.
#' In order to calculate the Q probability for each step
#' the distribution of turns and lifts to target and
#' the distribution of distance to target has to be known.
#' They can be derived from the empirical data (ideally),
#' or estimated from an unconditional process with the same properties.
#' Creates a unconditional empirical random walk, with a specific starting point,
#' geometrically similar to the initial trajectory.
#'
#' @section Note:
#' Simulations connecting start and end points
#' with more steps than 1/10th or more of the number of steps
#' of the empirical data should rather rely on simulated
#' unconditional walks with the same properties than on
#' the empirical data (\code{factor = 1500}).
#'
#' @section Random initial heading:
#' For a random initial heading a0 use:
#'   \code{sample(atan2(diff(coordinates(track)[,2]), diff(coordinates(track)[,1])),1)}
#'
#' @param n.locs the number of locations for the simulated track
#' @param start vector indicating the start point \code{c(x,y,z)}
#' @param a0 initial heading in radian
#' @param g0 initial gradient/polar angle in radian
#' @param densities list object returned by the \link[eRTG3D]{get.densities.3d} function
#' @param error logical: add random noise to the turn angle, lift angle and step length to account for errors measurements?
#'
#' @return A 3 dimensional trajectory in the form of a data.frame
#' @export
#'
#' @examples
#' sim.uncond.3d(
#'   10, start = c(0, 0, 0), a0 = pi / 2, g0 = pi / 2,
#'   densities = get.track.densities.3d(niclas)
#' )
sim.uncond.3d <- function(n.locs, start = c(0, 0, 0), a0, g0, densities, error = TRUE) {
  # progress bar and time
  message(paste("  |Simulate UERW with ", n.locs, " steps", sep = ""))
  start.time <- Sys.time()
  pb <- utils::txtProgressBar(min = 0, max = n.locs, style = 3)
  ui <- floor(n.locs / 20) + 1
  # get coordinates of the tldCube
  ts <- densities$tldCube$values$turn
  ls <- densities$tldCube$values$lift
  ds <- densities$tldCube$values$step
  # get probs for each turn-lift-distance combination
  tldProbs <- densities$tldCube$values$prob
  sCond <- sample(seq_len(nrow(densities$tldCube$values)), 1, prob = tldProbs)
  # "x" "y" "z" "a" "g" "t" "l" "d" "p"
  # "1" "2" "3" "4" "5" "6" "7" "8" "9"
  RTG <- matrix(0, n.locs, 9)
  RTG[1, ] <- c(start[1], start[2], start[3], a0, g0, ts[sCond], ls[sCond], ds[sCond], NA)
  # Create random noise if error is TRUE (uniform distributed)
  if (error) {
    tShift <- stats::runif(n.locs, -densities$tldCube$tRes / 2, densities$tldCube$tRes / 2)
    lShift <- stats::runif(n.locs, -densities$tldCube$lRes / 2, densities$tldCube$lRes / 2)
    dShift <- stats::runif(n.locs, -densities$tldCube$dRes / 2, densities$tldCube$dRes / 2)
  } else {
    tShift <- lShift <- dShift <- numeric(n.locs)
  }
  for (i in 2:n.locs) {
    # get influence of current autodifferences
    atProbs <- densities$autoT(RTG[i - 1, 6] - ts + tShift[i])
    alProbs <- densities$autoL(RTG[i - 1, 7] - ls + lShift[i])
    adProbs <- densities$autoD(RTG[i - 1, 8] - ds + dShift[i])
    atProbs[is.na(atProbs)] <- 0
    alProbs[is.na(alProbs)] <- 0
    adProbs[is.na(adProbs)] <- 0
    atProbs <- atProbs / sum(atProbs)
    alProbs <- alProbs / sum(alProbs)
    adProbs <- adProbs / sum(adProbs)
    # multiply and take the third squareroot
    pProbs <- tldProbs * (atProbs * alProbs * adProbs)^(1 / 3)
    # Apply gradient distribution or if gradientDensity is set to FALSE in get.densities.3d(), limit gradient to 0-pi.
    gAll <- (RTG[i - 1, 5] + ls + lShift[i])
    gProbs <- densities$gDens(gAll)
    gProbs[is.na(gProbs)] <- 0
    gProbs <- gProbs / sum(gProbs)
    pProbs <- pProbs * gProbs
    # sample on turnLiftStepHist = tldCube and add shifts
    rP <- sample(seq_len(nrow(densities$tldCube$values)), size = 1, prob = pProbs)
    t <- ts[rP] + tShift[i]
    l <- ls[rP] + lShift[i]
    d <- ds[rP] + dShift[i]
    p <- pProbs[rP]
    # absolute spherical orientation, wrap angles around -pi-0 & 0-pi
    a <- .wrap(RTG[i - 1, 4] + t)
    g <- .wrap(RTG[i - 1, 5] + l)
    # new coordinates of the next step
    x <- round((d * sin(g) * cos(a)) + RTG[i - 1, 1], 12)
    y <- round((d * sin(g) * sin(a)) + RTG[i - 1, 2], 12)
    z <- round((d * cos(g)) + RTG[i - 1, 3], 12)
    # "x" "y" "z" "a" "g" "t" "l" "d" "p"
    RTG[i, ] <- c(x, y, z, a, g, t, l, d, p)
    # update progress bar
    if ((i %% ui) == 0) {
      utils::setTxtProgressBar(pb, i)
    }
  }
  rownames(RTG) <- c()
  colnames(RTG) <- c("x", "y", "z", "a", "g", "t", "l", "d", "p")
  # close progress bar
  utils::setTxtProgressBar(pb, n.locs)
  close(pb)
  message(sprintf("  |Elapsed time: %ss", round(as.numeric(Sys.time()) - as.numeric(start.time), 1)))
  return(as.data.frame(RTG))
}

#' Q probabilities for n steps
#'
#' Calculates the Q probability, representing the pull to
#' the target. The number of steps on which the Q prob will be
#' quantified is number of total segments less than one
#' (the last step is defined by the target itself).
#'
#' @param sim the result of \link[eRTG3D]{sim.uncond.3d}, or a data frame with at least
#'     x,y,z-coordinates, the arrival azimuth and the arrival gradient.
#' @param n.locs number of total segments to be modeled,
#'     the length of the desired conditional empirical random walk
#' @param parallel logical: run computations in parallel (n-1 cores)? Or numeric: the number of nodes (maximum: n - 1 cores)
#' @param maxBin numeric scalar, maximum number of bins per dimension of the tld-cube (\link[eRTG3D]{turnLiftStepHist})
#'
#' @return A list containing the Q - tldCubes for every step
#' @export
#'
#' @examples
#' qProb.3d(niclas, 3)
qProb.3d <- function(sim, n.locs, parallel = FALSE, maxBin = 25) {
  nNodes <- .nNodes(parallel)
  message(paste("  |Extracting Q probabilities for ", n.locs, " steps", sep = ""))
  sim <- track.properties.3d(sim)
  nSteps <- n.locs - 2
  # lift angles to target as a function of number of steps
  cubeList <- parpblapply(1:nSteps, function(x) {
    # turn angle, lift angles and distance to target as a function of number of steps
    t <- .wrap(atan2(diff(sim$y, lag = x), diff(sim$x, lag = x)) - sim$a[1:(length(sim$a) - x)])
    l <- .wrap(atan2(
      sqrt(diff(sim$x, lag = x)^2 + diff(sim$y, lag = x)^2),
      diff(sim$z, lag = x)
    ) - sim$g[1:(length(sim$g) - x)])
    d <- sqrt(diff(sim$x, lag = x)^2 + diff(sim$y, lag = x)^2 + diff(sim$z, lag = x)^2)
    # the Qprob is thinned to the lag that suggests breaking off of the autocorrelation
    # of the turning angle to target, the lift angle to target and the distance to target
    # for the relevant number of steps. This is mainly to reduce redundancy mainly
    # introduced by the sliding window approach adopted in estimating the relationships
    k <- max(
      utils::head(which(stats::acf(t, lag.max = nSteps, plot = FALSE)$acf < 0.05), 1) - 1,
      utils::head(which(stats::acf(l, lag.max = nSteps, plot = FALSE)$acf < 0.05), 1) - 1,
      utils::head(which(stats::acf(d, lag.max = nSteps, plot = FALSE)$acf < 0.05), 1) - 1
    )
    t <- t[seq(1, length(t), by = k)]
    l <- l[seq(1, length(l), by = k)]
    d <- d[seq(1, length(d), by = k)]
    # get stepTurnLiftHistograms
    return(turnLiftStepHist(turn = t, lift = l, step = d, printDims = FALSE, rm.zeros = TRUE, maxBin = maxBin))
  }, export = c("sim", "nSteps", "maxBin"), packages = c("eRTG3D"), envir = environment(), nNodes = nNodes)
  return(rev(cubeList))
}

#' Conditional Empirical Random Walk (CERW) in 3-D
#'
#' Creates a conditional empirical random walk, with a specific starting and ending point,
#' geometrically similar to the initial trajectory
#' (extractMethod: raster overlay method can take "simple" or "bilinear")
#'
#' @param n.locs length of the trajectory in locations
#' @param start numeric vector of length 3 with the coordinates of the start point
#' @param end numeric vector of length 3 with the coordinates of the end point
#' @param a0 initial incoming heading in radian
#' @param g0 initial incoming gradient/polar angle in radian
#' @param densities list object returned by the \link[eRTG3D]{get.densities.3d} function
#' @param qProbs list object returned by the \link[eRTG3D]{qProb.3d} function
#' @param error logical: add random noise to the turn angle, lift angle and step length to account for errors measurements?
#' @param DEM raster layer containing a digital elevation model, covering the area between start and end point
#' @param BG a background raster layer that can be used to inform the choice of steps
#'
#' @return A trajectory in the form of data.frame
#' @export
#'
#' @examples
#' niclas <- track.properties.3d(niclas)
#' n.locs <- 3
#' P <- get.track.densities.3d(niclas)
#' f <- 1500
#' start <- Reduce(c, niclas[1, 1:3])
#' end <- Reduce(c, niclas[n.locs, 1:3])
#' a0 <- niclas$a[1]
#' g0 <- niclas$g[1]
#' uerw <- sim.uncond.3d(
#'   n.locs * f, start = start, a0 = a0, g0 = g0, densities = P
#' )
#' Q <- qProb.3d(uerw, n.locs)
#' sim.cond.3d(
#'   n.locs = n.locs, start = start, end = end,
#'   a0 = a0, g0 = g0, densities = P, qProbs = Q
#' )
sim.cond.3d <- function(n.locs, start = c(0, 0, 0), end = start, a0, g0, densities, qProbs, error = FALSE, DEM = NULL, BG = NULL) {
  start.time <- Sys.time()
  if (!is.null(DEM)) {
    .check.extent(DEM = DEM, track = data.frame(rbind(start, end)))
  }
  if (!is.null(BG)) {
    .check.extent(DEM = BG, track = data.frame(rbind(start, end)))
  }
  # progress bar and time
  message(paste("  |Simulate CERW with ", n.locs, " steps", sep = ""))
  pb <- utils::txtProgressBar(min = 0, max = n.locs - 2, style = 3)
  ui <- floor(n.locs / 20) + 1
  # replace the probability distribution for step length 1 by the one from
  # the qProbs since that one relies on more samples derived from sim
  densities[[1]] <- utils::tail(qProbs, 1)[[1]]
  # get the coordinates of the step length and turning angle bin centres
  start <- Reduce(c, start)
  names(start) <- c("x", "y", "z")
  end <- Reduce(c, end)
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
  sCond <- sample(seq_len(nrow(densities$tldCube$values)), 1, prob = tldProbs)
  # "x" "y" "z" "a" "g" "t" "l" "d" "p"
  # "1" "2" "3" "4" "5" "6" "7" "8" "9"
  RTG <- matrix(0, n.locs, 9)
  RTG[1, ] <- c(start[1], start[2], start[3], a0, g0, ts[sCond], ls[sCond], ds[sCond], NA)
  # Create random noise if error is TRUE
  if (error) {
    tShift <- stats::runif(n.locs - 2, -densities$tldCube$tRes / 2, densities$tldCube$tRes / 2)
    lShift <- stats::runif(n.locs - 2, -densities$tldCube$lRes / 2, densities$tldCube$lRes / 2)
    dShift <- stats::runif(n.locs - 2, -densities$tldCube$dRes / 2, densities$tldCube$dRes / 2)
  } else {
    tShift <- lShift <- dShift <- numeric(n.locs - 2)
  }
  # start creating the track step for step
  for (i in 1:(n.locs - 2)) {
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
    P <- (tldProbs) * (atProbs * alProbs * adProbs)^(1 / 3)
    # calculate the azimuth
    a <- .wrap(RTG[i, 4] + ts + tShift[i])
    # calculate the gradient
    g <- .wrap(RTG[i, 5] + ls + lShift[i])
    # convert the coordinates from step length turning angle dimension
    x1 <- round(((ds + dShift[i]) * sin(g) * cos(a)) + RTG[i, 1], 12)
    y1 <- round(((ds + dShift[i]) * sin(g) * sin(a)) + RTG[i, 2], 12)
    z1 <- round(((ds + dShift[i]) * cos(g)) + RTG[i, 3], 12)
    # calculate the distances of the cell centers in the spatial domain
    # to the target (last location of the empirical track)
    endD <- as.numeric(sqrt((end[1] - x1)^2 + (end[2] - y1)^2 + (end[3] - z1)^2))
    # calculate the azimuth of the cell centres to the target and substract from it the direction of arrival
    # resulting in turning angle towards target
    endT <- as.numeric(.wrap(atan2(as.numeric(end[2] - y1), as.numeric(end[1] - x1)) - a))
    # calculate the gradient of the possibilite steps to the target and substract from it the angle of arrival
    # resulting in turning angle towards target
    endL <- as.numeric(.wrap(atan2(as.numeric(sqrt((end[1] - x1)^2 + (end[2] - y1)^2)), as.numeric(end[3] - z1))) - g)
    # get the probabilities of making it distance and turning angle wise
    # which is derived from the two dimensional probability distribution for the
    # appropriate step being modelled
    # get possible coordinates
    qCube <- qProbs[[i]]
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
    gProbs <- densities$gDens(g)
    gProbs[is.na(gProbs)] <- 0
    gProbs <- gProbs / sum(gProbs)
    Probs <- Probs * gProbs
    # Account for probable flight height, if a DEM is provided the relative flight height is taken
    # Otherwise only the absolute ellipsoid height.
    if (!is.null(DEM)) {
      surface <- raster::extract(DEM, cbind(x1, y1))
      demP <- densities$hDistTopo(z1 - surface) * as.numeric(z1 >= surface)
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
      RTG <- NULL
      close(pb)
      message(sprintf("  |Elapsed time: %ss", round(as.numeric(Sys.time()) - as.numeric(start.time), 1)))
      warning("Dead end encountered.")
      return(RTG)
    } else {
      # draw a point randomly based on the probability
      rP <- sample.int(nrow(densities$tldCube$values), size = 1, prob = Probs)
      # "x" "y" "z" "a" "g" "t" "l" "d" "p"
      # "1" "2" "3" "4" "5" "6" "7" "8" "9"
      RTG[i + 1, ] <- c(x1[rP], y1[rP], z1[rP], a[rP], g[rP], ts[rP], ls[rP], ds[rP], Probs[rP])
      # update progress bar
      if ((i %% ui) == 0) {
        utils::setTxtProgressBar(pb, i)
      }
    }
  }
  # the track is forced to target location and the appropriate distance is added
  RTG[1, 8] <- NA
  RTG[n.locs, ] <- c(end[1], end[2], end[3], NA, NA, NA, NA, NA, NA)
  RTG[n.locs, 8] <- sqrt((RTG[n.locs, 1] - RTG[n.locs - 1, 1])^2 +
    (RTG[n.locs, 2] - RTG[n.locs - 1, 2])^2 +
    (RTG[n.locs, 3] - RTG[n.locs - 1, 3])^2)
  # Stop if the step length of the last step is larger than the largest possible step
  if (RTG[n.locs, 8] > max(densities$tldCube$values$step, na.rm = TRUE) * sqrt(2)) {
    RTG <- NULL
    close(pb)
    message(sprintf("  |Elapsed time: %ss", round(as.numeric(Sys.time()) - as.numeric(start.time), 1)))
    warning("Dead end encountered in last step.")
    return(RTG)
  }
  rownames(RTG) <- c()
  colnames(RTG) <- c("x", "y", "z", "a", "g", "t", "l", "d", "p")
  # close progress bar
  utils::setTxtProgressBar(pb, i)
  close(pb)
  message(sprintf("  |Elapsed time: %ss", round(as.numeric(Sys.time()) - as.numeric(start.time), 1)))
  return(as.data.frame(RTG))
}

#' Conditional Empirical Random Walks (CERW) in 3-D
#'
#' Creates multiple conditional empirical random walks, with a specific starting and ending point,
#' geometrically similar to the initial trajectory by applying \link[eRTG3D]{sim.cond.3d} multiple times.
#'
#' @param n.sim number of CERWs to simulate
#' @param n.locs length of the trajectory in locations
#' @param start numeric vector of length 3 with the coordinates of the start point
#' @param end numeric vector of length 3 with the coordinates of the end point
#' @param a0 initial incoming heading in radian
#' @param g0 initial incoming gradient/polar angle in radian
#' @param densities list object returned by the \link[eRTG3D]{get.densities.3d} function
#' @param qProbs list object returned by the \link[eRTG3D]{qProb.3d} function
#' @param error logical: add random noise to the turn angle, lift angle and step length to account for errors measurements?
#' @param DEM raster layer containing a digital elevation model, covering the area between start and end point
#' @param BG a background raster layer that can be used to inform the choice of steps
#' @param parallel logical: run computations in parallel (n-1 cores)? Or numeric: the number of nodes (maximum: n - 1 cores)
#'
#' @return
#' A list containing the CERWs or \code{NULL}s if dead ends have been encountered.
#' @export
#'
#' @examples
#' niclas <- track.properties.3d(niclas)
#' n.locs <- 3
#' P <- get.track.densities.3d(niclas)
#' f <- 1500
#' start <- Reduce(c, niclas[1, 1:3])
#' end <- Reduce(c, niclas[n.locs, 1:3])
#' a0 <- niclas$a[1]
#' g0 <- niclas$g[1]
#' uerw <- sim.uncond.3d(
#'   n.locs * f, start = start, a0 = a0, g0 = g0, densities = P
#' )
#' Q <- qProb.3d(uerw, n.locs)
#' n.sim.cond.3d(
#'   n.sim = 2, n.locs = n.locs,
#'   start = start, end = end,
#'   a0 = a0, g0 = g0,
#'   densities = P, qProbs = Q
#' )
n.sim.cond.3d <- function(n.sim, n.locs, start = c(0, 0, 0), end = start, a0, g0, densities, qProbs, error = FALSE, parallel = FALSE, DEM = NULL, BG = NULL) {
  n.sim <- round(n.sim)
  if (n.sim <= 1) {
    return(sim.cond.3d(n.locs = n.locs, start = start, end = end, a0 = a0, g0 = g0, densities = densities, qProbs = qProbs, error = error, DEM = DEM, BG = BG))
  }
  nNodes <- .nNodes(parallel)
  message(paste("  |Simulate ", n.sim, " CERWs with ", n.locs, " steps", sep = ""))
  cerwList <- parpblapply(
    X = 1:n.sim, FUN = function(x) {
      sim.cond.3d(n.locs = n.locs, start = start, end = end, a0 = a0, g0 = g0, densities = densities, qProbs = qProbs, error = error, DEM = DEM, BG = BG)
    },
    export = c("n.locs", "start", "end", "a0", "g0", "densities", "qProbs", "error", "DEM", "BG"), packages = c("eRTG3D"), nNodes = nNodes, envir = environment()
  )
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
#' @noRd
.wrap <- function(x) {
  (x + pi) %% (2 * pi) - pi
}
