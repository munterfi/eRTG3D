## ---- echo = FALSE, eval=TRUE, include=FALSE-----------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
options(knitr.table.format = "html")

library(eRTG3D)
set.seed(123)

## ----eval=FALSE----------------------------------------------------------
#  niclas <- track.properties.3d(niclas)

## ----echo=FALSE, results = "asis", fig.height=5, fig.width=10------------
niclas <- track.properties.3d(niclas)
pander::pandoc.table(head(round(niclas, 2),5))

## ----eval=FALSE----------------------------------------------------------
#  P <- get.track.densities.3d(niclas, heightDistEllipsoid = TRUE, DEM = dem)

## ----eval=FALSE----------------------------------------------------------
#  sim.locs <- nrow(niclas)
#  f <- 1500
#  uerw <- sim.uncond.3d(sim.locs*f, start = c(niclas$x[1], niclas$y[1], niclas$z[1]),
#                        a0 = niclas$a[1], g0 = niclas$g[1], densities = P)
#  Q <- qProb.3d(uerw, sim.locs)

## ----eval=FALSE----------------------------------------------------------
#  start <- c(niclas$x[1], niclas$y[1], niclas$z[1])
#  end <- c(niclas$x[nrow(niclas)], niclas$y[nrow(niclas)], niclas$z[nrow(niclas)])
#  a0 <- niclas$a[1]
#  g0 <- niclas$g[1]

## ----eval=FALSE----------------------------------------------------------
#  cerw <- sim.cond.3d(sim.locs, start = start, end = end, a0 = a0, g0 = g0,
#                      densities = P, qProbs = Q, DEM = dem)

## ----eval=FALSE----------------------------------------------------------
#  cerwList <- n.sim.cond.3d(n.sim = 100, sim.locs,
#                            start = start, end = end, a0 = a0, g0 = g0,
#                            densities = P, qProbs = Q, DEM = dem)

## ----eval=FALSE----------------------------------------------------------
#  trackSections <- track.split.3d(track, timeLag)
#  P <- get.section.densities.3d(trackSections, DEM = dem)

## ---- fig.show='hold', eval=FALSE----------------------------------------
#  cerwList <- reproduce.track.3d(n.sim = 100, niclas, DEM = dem)

