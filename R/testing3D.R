#' Statistical Verification of the simulated track
#'
#' Uses two-sample Kolmogorov-Smirnov test or the one-sample t-test to compare the geometric characteristics of the original track
#' with the characteristics of the simulated track.
#'
#' @param track1 data.frame or list of data.frames with x,y,z coordinates of the original track
#' @param track2 data.frame or list of data.frames with x,y,z coordinates of the simulated track
#' @param alpha scalar: significance level, default \code{alpha = 0.05}
#' @param plot logical: plot the densities or differences of turn angle, lift angle and step length of the two tracks?
#' @param test character: either \code{"ks"} or \code{"ttest"} to choose the kind of test procedure.
#'
#' @return Test objects of the 6 two-sample Kolmogorov-Smirnov test conducted.
#' @export
#' 
#' @note By choosing \code{test = "ttest"} a random sample, without replacement is taken from the longer track,
#' to shorten it to the length of the longer track. The order of the shorter track is also sampled randomly.
#' Then the two randomly ordered vectors of turn angles, lift angles and step lengths are substracted from each other.
#' If the both tracks stem from the same distributions the the mean deviatio should tend to towards zero, therefore the 
#' difference is tested two-sided against \code{mu = 0} with a one-sample t-test.
#' 
#' By setting \code{test = "ks"} a two-sample Kolmogorov-Smirnov test is carried out on the distributions of turn angles,
#' lift angles and step lengths of the two tracks.
#'
#' @examples
#' test.verification.3d(niclas, niclas)
test.verification.3d <- function(track1, track2, alpha = 0.05, plot = FALSE, test = "ks")
{
  if (!any(test == c("ks", "ttest"))) stop("The variable 'test' must either be 'ks' or 'ttest'.")
  if (!is.list(track1) || !is.list(track2)) stop("Track input has to be of type list or data.frame.")
  if (is.list(track1) && is.data.frame(track1)) {track1 <- list(track1)}
  if (is.list(track2) && is.data.frame(track2)) {track2 <-list(track2)}
  track1 <- filter.dead.ends(track1); track2 <- filter.dead.ends(track2)
  # track(s) 1
  track1 <- lapply(track1, function(x){track.properties.3d(x)[2:nrow(x), ]})
  difftrack1 <- do.call("rbind", lapply(track1, function(x){data.frame(diffT = diff(x$t), diffL = diff(x$l), diffD = diff(x$d))}))
  track1 <- do.call("rbind", track1)
  t1 <- track1$t; l1 <- track1$l; d1 <- track1$d;
  diffT1 <- difftrack1$diffT; diffL1 <- difftrack1$diffL; diffD1 <- difftrack1$diffD;
  # track(s) 2
  track2 <- lapply(track2, function(x){track.properties.3d(x)[2:nrow(x), ]})
  diffTrack2 <- do.call("rbind", lapply(track2, function(x){data.frame(diffT = diff(x$t), diffL = diff(x$l), diffD = diff(x$d))}))
  track2 <- do.call("rbind", track2)
  t2 <- track2$t; l2 <- track2$l; d2 <- track2$d;
  diffT2 <- diffTrack2$diffT; diffL2 <- diffTrack2$diffL; diffD2 <- diffTrack2$diffD;
  if (test == "ks"){
    message("  |*** Two-sample Kolmogorov-Smirnov test ***")
    message("  |H0: Probability distributions do not differ significantly")
    message("  |H1: Probability distributions differ significantly")
    # turn
    turnT <- suppressWarnings(stats::ks.test(t1, t2, alternative = "two.sided"))
    diffTurnT <- suppressWarnings(stats::ks.test(diffT1, diffT2, alternative = "two.sided"))
    # lift
    liftT <- suppressWarnings(stats::ks.test(l1, l2, alternative = "two.sided"))
    diffLiftT <- suppressWarnings(stats::ks.test(diffL1, diffL2, alternative = "two.sided"))
    # step
    stepT <- suppressWarnings(stats::ks.test(d1, d2, alternative = "two.sided"))
    diffStepT <- suppressWarnings(stats::ks.test(diffD1, diffD2, alternative = "two.sided"))
    # print results
    message(paste("  |Turn angle  - ", .test2text(turnT, alpha), ", autodifferences - ", .test2text(diffTurnT, alpha), sep=""))
    message(paste("  |Lift angle  - ", .test2text(liftT, alpha), ", autodifferences - ", .test2text(diffLiftT, alpha), sep=""))
    message(paste("  |Step length - ", .test2text(stepT, alpha), ", autodifferences - ", .test2text(diffStepT, alpha), sep=""))
    if (plot) {
      suppressWarnings(plot3d.multiplot(
        .plot3d.density(t1, t2, titleText = "Turn angle"),
        .plot3d.density(l1, l2, titleText = "Lift angle"),
        .plot3d.density(d1, d2, titleText = "Step length"),
        cols = 1
      ))
    }
    return(list(turnT, liftT, stepT, diffTurnT, diffLiftT, diffStepT))
  }
  if (test == "ttest"){
    message("  |*** One Sample t-test ***")
    message("  |H0: Difference between tracks does not differ significantly from 0")
    message("  |H1: Difference between tracks differs significantly from 0")
    nSample <- min(nrow(track1), nrow(track2))
    # turn
    turnT <- suppressWarnings(stats::t.test(diffT <- (sample(t1, nSample)-sample(t2, nSample)), mu = 0, alternative = "two.sided"))
    # lift
    liftT <- suppressWarnings(stats::t.test(diffL <- (sample(l1, nSample)-sample(l2, nSample)), mu = 0, alternative = "two.sided"))
    # step
    stepT <- suppressWarnings(stats::t.test(diffD <- (sample(d1, nSample)-sample(d2, nSample)), mu = 0, alternative = "two.sided"))
    # print results
    message(paste("  |Turn angle  - ", .test2text(turnT, alpha), sep=""))
    message(paste("  |Lift angle  - ", .test2text(liftT, alpha), sep=""))
    message(paste("  |Step length - ", .test2text(stepT, alpha), sep=""))
    if (plot) {
      suppressWarnings(plot3d.multiplot(
        .plot3d.density(diffT, titleText = "Mean difference turn angle"),
        .plot3d.density(diffL, titleText = "Mean difference Lift angle"),
        .plot3d.density(diffD, titleText = "Mean difference Step length"),
        cols = 1
      ))
    }
    return(list(turnT, liftT, stepT))
  }
}

#' Extract test results as string
#'
#' @param test object of type \code{htest}
#' @param alpha scalar: significance level, default \code{alpha = 0.05}
#'
#' @return A character describing the results.
#' @export
#'
#' @examples
#' .test2text(test, alpha)
#' @noRd
.test2text <- function(test, alpha)
{
  p <- test$p.value
  paste("p-value: ", round(p,3) , if(p<alpha){
    paste(" < ", alpha, ", *H1*", sep = "")
  } else {
    paste(" > ", alpha, ", *H0*", sep = "")
  }, sep = "")
}

#' Test the functionality of the eRTG3D
#'
#' The test simulates a CRW with given parameters and reconstructs it by using the eRTG3D
#'
#' @param parallel logical: test running parallel?
#' @param returnResult logical: return tracks generated?
#' @param plot2d logical: plot tracks on 2-D plane?
#' @param plot3d logical: plot tracks in 3-D?
#' @param plotDensities logical: plot densities of turning angle, lift angle and step length?
#'
#' @return A list containing the original CRW and the simulated track (CERW).
#' @export
#'
#' @examples
#' \donttest{
#' test.eRTG.3d()
#' }
test.eRTG.3d <- function(parallel = FALSE, returnResult = FALSE, plot2d = FALSE, plot3d = TRUE, plotDensities = TRUE)
{
  message("  |*** Testing eRTG3D ***")
  set.seed(1)
  nStep <- 25
  crw <- track.properties.3d(
    sim.crw.3d(nStep = nStep, rTurn = 0.99, rLift = 0.99, meanStep = 1, start = c(0, 0, 10)))
  turnAngle <- crw$t[2:nrow(crw)]; liftAngle <- crw$l[2:nrow(crw)]; stepLength <- crw$d[2:nrow(crw)]
  deltaTurn <- diff(turnAngle); deltaLift <- diff(liftAngle); deltaStep <- diff(stepLength)
  heightEllipsoid <- crw$z
  D <- get.densities.3d(liftAngle = liftAngle, turnAngle = turnAngle, stepLength = stepLength,
                        deltaLift = deltaLift, deltaTurn = deltaTurn, deltaStep = deltaStep,
                        heightEllipsoid = heightEllipsoid, heightTopo = NULL)
  uerw <- sim.uncond.3d(nStep*1500, start = c(crw$x[1],crw$y[1],crw$z[1]),
                        a0 = crw$a[1], g0 = crw$g[1], densities = D)
  tests.uerw <- test.verification.3d(crw, uerw, alpha = 0.05)
  if(parallel) {
    nNodes <- .nNodes(parallel)
    Q <- qProb.3d(uerw, nStep, parallel = nNodes)
    cerw <- n.sim.cond.3d(n.sim = 100, n.locs = nStep, start=c(crw$x[1],crw$y[1],crw$z[1]), end=c(crw$x[nStep],crw$y[nStep],crw$z[nStep]),
                        a0 = crw$a[1], g0 = crw$g[1], densities=D, qProbs=Q, parallel = nNodes)
    cerw <- filter.dead.ends(cerw)
  } else {
    Q <- qProb.3d(uerw, nStep, parallel = FALSE)
    cerw <- sim.cond.3d(nStep, start=c(crw$x[1],crw$y[1],crw$z[1]), end=c(crw$x[nStep],crw$y[nStep],crw$z[nStep]),
                        a0 = crw$a[1], g0 = crw$g[1], densities=D, qProbs=Q)
  }
  tests.cerw <- test.verification.3d(crw, cerw, alpha = 0.05)
  if(plot2d){print(plot2d(crw, cerw))}
  if(plotDensities){plot3d.densities(crw, cerw)}
  if(plot3d){plot3d(crw, cerw)}
  if(returnResult){return(list(crw = crw, cerw = cerw))}
  message("  |*** Test passed successfully ***")
}

#' Simulation of a three dimensional Correlated Random Walk
#'
#' @param nStep the number of steps of the simulated trajectory
#' @param rTurn the correlation on the turn angle
#' @param rLift the correlation of the lift angle
#' @param meanStep the mean step length
#' @param start a vector of length 3 containing the coordinates of the start point of the trajectory
#'
#' @return A trajectory in the form of data.frame
#' @export
#'
#' @examples
#' sim.crw.3d(nStep=10, rTurn=0.9, rLift=0.9, meanStep=1, start = c(0,0,0))
sim.crw.3d <- function(nStep, rTurn, rLift, meanStep, start = c(0,0,0))
{
  # correlated angles and distance
  t <- CircStats::rwrpnorm(n = nStep - 2, mu = 0, rho = rTurn)
  a <- .wrap(cumsum(c(stats::runif(1, 0, 2 * pi), t)))
  l <- CircStats::rwrpnorm(n = nStep - 2, mu = 0, rho = rLift)
  g <- abs(.wrap(cumsum(c(stats::runif(1, 0, pi), l))))
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
