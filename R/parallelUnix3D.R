#' Parallel or single core Q probabilities for n steps on unix (and windows single core)
#'
#'  Calculates the Q probability, representing the pull to
#' the target. The number of steps on which the Q prob will be
#' quantified is number of total segments less than one
#' (the last step is defined by the target itself).
#'
#' @param sim the result of simm.uncond.3d(), or a data frame with at least
#'     x,y,z-coordinates, the arrival azimuth and the arrival gradient.
#' @param n.locs number of total segments to be modelled,
#'     the length of the desired conditional empirical random walk
#' @param maxBin numeric scalar, maximum number of bins per dimension of the tld-cube (\link[eRTG3D]{turnLiftStepHist})
#'
#' @return A list containing the Q - tldCubes for every step
#' @export
#'
#' @examples
#' .qProb.3d.unix(sim, n.locs)
.qProb.3d.unix <- function(sim, n.locs, maxBin = 25)
{
  start.time <- Sys.time()
  nCores <- parallel::detectCores()-1
  message(paste("  |Extracting Q probabilities for ", n.locs, " steps (Parallel on nCores = ", nCores, ")", sep = ""))
  # steps minus 2
  nSteps <- n.locs - 2
  sim <- track.properties.3d(sim)
  # lift angles to target as a function of number of steps
  cubeList <- pbmcapply::pbmclapply(1:nSteps, function(x) {
    # turn angle, lift angles and distance to target as a function of number of steps
    t <- .wrap(atan2(diff(sim$y, lag = x), diff(sim$x, lag = x)) - sim$a[1:(length(sim$a) - x)])
    l <- .wrap(atan2(sqrt(diff(sim$x, lag = x) ^ 2 + diff(sim$y, lag = x) ^ 2),
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
  }, mc.cores = nCores, mc.style = "txt")
  message(paste("  |Runtime: ", round(as.numeric(Sys.time()) - as.numeric(start.time), 2), " secs", sep = ""))
  return(rev(cubeList))
}

#' Parallel computation of n Conditional Empirical Random Walks (CERW) in 3-D on Unix
#'
#' Creates n conditional empirical random walks, with a specific starting and ending point,
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
#' .n.sim.cond.3d.unix(n.sim, n.locs, start = c(0,0,0), end=start, a0, g0, densities, qProbs)
.n.sim.cond.3d.unix <- function(n.sim, n.locs, start = c(0,0,0), end=start, a0, g0, densities, qProbs, error = FALSE, DEM = NULL, BG = NULL)
{
  nCores <- parallel::detectCores()-1
  message(paste("  |Running on nCores = ", nCores, sep=""))
  return(pbmcapply::pbmclapply(X = 1:n.sim, FUN = function(x){sim.cond.3d(n.locs, start, end, a0, g0, densities, qProbs, error, DEM, BG)}, mc.cores = nCores, mc.style = "txt"))
}