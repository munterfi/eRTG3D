## ---- echo = FALSE, eval=TRUE, include=FALSE-----------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
options(knitr.table.format = "html")

library(eRTG3D)
set.seed(123)
cerwList <- reproduce.track.3d(n.sim = 100, multicore=TRUE, niclas, DEM = dem, BG = (dem<650), filterDeadEnds = TRUE)
cerw <- filter.dead.ends(cerwList)[[1]]

## ---- eval = FALSE, fig.show='hold'--------------------------------------
#  library(eRTG3D)
#  test.eRTG.3d()

## ---- eval = FALSE, fig.show='hold'--------------------------------------
#  test.eRTG.3d(multicore = TRUE)

## ---- eval = FALSE, fig.show='hold'--------------------------------------
#  results <- test.eRTG.3d(returnResult = TRUE, plot2d = TRUE, plot3d = TRUE)

## ----eval=FALSE----------------------------------------------------------
#  niclas <- track.properties.3d(niclas)

## ----echo=FALSE, results = "asis", fig.height=5, fig.width=10------------
niclas <- track.properties.3d(niclas)
niclas[c("dx", "dy", "dz", "d2d")] <- NULL
pander::pandoc.table(head(round(niclas, 2),5))

## ----eval=FALSE----------------------------------------------------------
#  D <- get.densities.3d(niclas, heightDistEllipsoid = TRUE, DEM = dem)

## ----eval=FALSE----------------------------------------------------------
#  sim.locs <- nrow(niclas)
#  f <- 1500
#  uerw <- sim.uncond.3d(sim.locs*f, start = c(niclas$x[1], niclas$y[1], niclas$z[1]),
#                        a0 = niclas$a[1], g0 = niclas$g[1], densities = D)
#  Q <- qProb.3d(uerw, sim.locs)

## ----eval=FALSE----------------------------------------------------------
#  start=c(niclas$x[1],niclas$y[1],niclas$z[1])
#  end=c(niclas$x[sim.locs],niclas$y[sim.locs],niclas$z[sim.locs])
#  a0 = niclas$a[1]
#  g0 = niclas$g[1]

## ----eval=FALSE----------------------------------------------------------
#  cerw <- sim.cond.3d(sim.locs, start=start, end=end, a0 = a0, g0 = g0, densities=D, qProbs=Q, DEM = dem)

## ----eval=FALSE----------------------------------------------------------
#  cerw <- n.sim.cond.3d(n.sim = 100, sim.locs, start=start, end=end, a0 = a0, g0 = g0, densities=D, qProbs=Q, DEM = dem)

## ---- fig.show='hold', eval=FALSE----------------------------------------
#  cerwList <- reproduce.track.3d(n.sim = 100, niclas, DEM = dem)

## ----eval=TRUE, fig.height=5, fig.width=7--------------------------------
tests <- test.verification.3d(niclas, cerw, plotDensities = FALSE)

## ----eval=TRUE, fig.height=5, fig.width=7--------------------------------
tests <- test.verification.3d(niclas, cerwList, plotDensities = FALSE)

## ----eval=TRUE, fig.height=5, fig.width=7--------------------------------
plot2d.densities(niclas, cerwList)

## ----eval=TRUE, fig.height=5, fig.width=7--------------------------------
plot2d(niclas, cerwList, titleText=paste("Steps: ", nrow(niclas), ", Niclas", sep=""), DEM=dem)

## ----eval=FALSE, fig.height=5, fig.width=7-------------------------------
#  plot3d(niclas, cerw, titleText=paste("Steps: ", nrow(niclas), ", Niclas", sep=""), surface=TRUE, DEM=dem)

## ----eval=FALSE, fig.height=5, fig.width=7-------------------------------
#  Q <- qProb.3d(uerw, sim.locs, multicore=TRUE)
#  cerwList <- reproduce.track.3d(n.sim = 100, niclas, DEM = dem, multicore=TRUE)

## ----eval=TRUE, echo=FALSE, fig.height=5, fig.width=7, include=FALSE-----
D <- get.densities.3d(niclas, heightDistEllipsoid = TRUE, DEM = dem)
sim.locs <- nrow(niclas)
f <- 1500
uerw <- sim.uncond.3d(sim.locs*f, start = c(niclas$x[1], niclas$y[1], niclas$z[1]), 
                      a0 = niclas$a[1], g0 = niclas$g[1], densities = D)
Q <- qProb.3d(uerw, sim.locs)
start=c(niclas$x[1],niclas$y[1],niclas$z[1])
end=c(niclas$x[sim.locs],niclas$y[sim.locs],niclas$z[sim.locs])
a0 = niclas$a[1]
g0 = niclas$g[1]

## ----eval=TRUE, fig.height=5, fig.width=7--------------------------------
cerw <- n.sim.cond.3d(n.sim = 100, sim.locs, start=start, end=end, a0 = a0, g0 = g0, densities=D, qProbs=Q, DEM = dem, multicore=TRUE)

