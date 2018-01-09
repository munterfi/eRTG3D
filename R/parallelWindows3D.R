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
.qProb.3d.windows <- function(sim, n.locs, maxBin = 25)
{
  start.time <- Sys.time()
  # set up cluster (parallel, doParallel and later plyr)
  nCores <- parallel::detectCores()-1
  message(paste("  |Extracting Q probabilities for ", n.locs, " steps (Parallel on nCores = ", nCores, ")", sep = ""))
  cl <- parallel::makeCluster(nCores)
  doParallel::registerDoParallel(nCores)
  opts <- list(preschedule=FALSE)
  # define variables and functions needed later to pass them to the cluster
  sim <- sim
  wrap <- function(x) {(x + pi) %% (2 * pi) - pi}
  turnLiftStepHist <- function(turn, lift, step, printDims = TRUE, rm.zeros = TRUE, maxBin = 25)
  {
    nx <- min(max(floor(2 * pi / fd.bw(turn)), 12), maxBin)
    ny <- min(max(floor(2 * pi / fd.bw(lift)), 12), maxBin)
    nz <- min(max(floor(max(step) / fd.bw(step)), 12), maxBin)
    if(printDims){message("  |Dims = x: ", nx, ", y: ", ny, ", z: ", nz)}
    tCuts <- cutMidpoints(turn, nx); lCuts <- cutMidpoints(lift, ny); dCuts <- cutMidpoints(step, nz)
    h <- list(turn=tCuts[[1]],
              lift=lCuts[[1]],
              step=dCuts[[1]])
    h <- do.call(data.frame, h)
    h <- as.data.frame(table(h))
    tRes <- tCuts[[2]]; lRes <- lCuts[[2]]; dRes <- dCuts[[2]];
    colnames(h)[4] <- "prob"
    if (rm.zeros) {h <- h[!h$prob==0, ]}
    h$prob <- h$prob/sum(h$prob)
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
  # progress bar
  pb <- txtProgressBar(min = 0, max = 18, style = 3)
  # steps minus 2
  nSteps <- n.locs - 2
  # turning angles to target as a function of number of steps
  tList <- plyr::alply(1:nSteps, function(x) wrap(atan2(diff(sim$y, lag = x),
                                                        diff(sim$x, lag = x)) - sim$a[1:(length(sim$a) - x)]), .margins = 1, .parallel = TRUE,
                       .paropts = list(.options.snow=opts))
  setTxtProgressBar(pb, 3)
  # lift angles to target as a function of number of steps
  lList <- plyr::alply(1:nSteps, function(x) wrap(atan2(sqrt(diff(sim$x, lag = x) ^ 2 + diff(sim$y, lag = x) ^ 2),
                                                        diff(sim$z, lag = x)) - sim$g[1:(length(sim$g) - x)]), .margins = 1, .parallel = TRUE,
                       .paropts = list(.options.snow=opts))
  setTxtProgressBar(pb, 6)
  # calculate distance to target as a function of number of steps
  dList <- plyr::alply(1:nSteps, function(x) sqrt(diff(sim$x, lag = x) ^ 2
                                                  + diff(sim$y, lag = x) ^ 2
                                                  + diff(sim$z, lag = x) ^ 2), .margins = 1, .parallel = TRUE,
                       .paropts = list(.options.snow=opts))
  setTxtProgressBar(pb, 9)
  # the Qprob is thinned to the lag that suggests breaking off of the autocorrelation
  # of the turning angle to target, the lift angle to target and the distance to target
  # for the relevant number of steps. This is mainly to reduce redundancy mainly
  # introduced by the sliding window approach adopted in estimating the relationships
  k <- cbind(unlist(lapply(lapply(lapply(lapply(plyr::llply(tList, acf, lag.max=nSteps, plot = FALSE, .parallel = TRUE,
                                                            .paropts = list(.options.snow=opts)), '[[','acf'),'<',.05),
                                  which), head, 1)) - 1,
             unlist(lapply(lapply(lapply(lapply(plyr::llply(lList, acf, lag.max=nSteps, plot = FALSE, .parallel = TRUE,
                                                            .paropts = list(.options.snow=opts)), '[[','acf'),'<',.05),
                                  which), head, 1)) - 1,
             unlist(lapply(lapply(lapply(lapply(plyr::llply(dList, acf, lag.max=nSteps, plot = FALSE, .parallel = TRUE,
                                                            .paropts = list(.options.snow=opts)), '[[','acf'),'<',.05),
                                  which), head, 1)) - 1)
  kk <- apply(k,1,max)
  setTxtProgressBar(pb, 14)
  tList <-mapply('[',tList,mapply(seq, 1, lapply(tList, length), by = kk))
  lList <-mapply('[',lList,mapply(seq, 1, lapply(lList, length), by = kk))
  dList <-mapply('[',dList,mapply(seq, 1, lapply(dList, length), by = kk))
  # Use multicore to speed the calculations up
  cubeList <- rev(plyr::llply(1:nSteps, function(x) turnLiftStepHist(turn=tList[[x]], lift=lList[[x]], step=dList[[x]], rm.zeros = TRUE, maxBin = maxBin),
                              .parallel = TRUE,  .paropts = list(.options.snow=opts)))
  # stop cluster
  parallel::stopCluster(cl)
  # complete progress bar and close
  setTxtProgressBar(pb, 18)
  close(pb)
  message("  |Minimum number of independent estimates: ", min(unlist(lapply(dList, length))), " for step ", which.min(unlist(lapply(dList, length))), ".")
  message(paste("  |Runtime: ", round(as.numeric(Sys.time()) - as.numeric(start.time), 2), " secs", sep = ""))
  return(cubeList)
}

#' Parallel computation of n Conditioned Empirical Random Walks (CERW) in 3D on Windows
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
.n.sim.cond.3d.windows <- function(n.sim, n.locs, start = c(0,0,0), end=start, a0, g0, densities, qProbs, error = FALSE, DEM = NULL, BG = NULL)
{
  warning("Parallel version not yet supported on Windows. Please set 'multicore' to 'FALSE' or change to a unix system.")
  return(NULL)
}
