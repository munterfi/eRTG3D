## ---- echo = FALSE, eval=TRUE, include=FALSE-----------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
options(knitr.table.format = "html")

library(eRTG3D)
set.seed(123)

## ----eval=FALSE, echo=FALSE, fig.height=5, fig.width=7, include=FALSE----
#  niclas <- track.properties.3d(niclas)
#  sim.locs <- nrow(niclas)
#  f <- 1500
#  a0 = niclas$a[1]
#  g0 = niclas$g[1]

## ----eval=FALSE----------------------------------------------------------
#  P <- get.track.densities.3d(niclas, heightDistEllipsoid = TRUE, DEM = dem)
#  uerw <- sim.uncond.3d(sim.locs*f, start = c(niclas$x[1], niclas$y[1], niclas$z[1]),
#                        a0 = a0, g0 = g0, densities = P)

## ----eval=FALSE, fig.height=5, fig.width=7-------------------------------
#  Q <- qProb.3d(uerw, sim.locs, parallel = TRUE)
#  cerwList <- reproduce.track.3d(n.sim = 100, niclas, DEM = dem, parallel = TRUE)

## ----eval=FALSE, fig.height=5, fig.width=7, quiet = TRUE-----------------
#  cerwList <- n.sim.cond.3d(n.sim = 100, sim.locs, start=start, end=end,a0 = a0, g0 = g0,
#                        densities=P, qProbs=Q, DEM = dem, parallel = TRUE)

## ----eval=FALSE, fig.height=5, fig.width=7, quiet = TRUE-----------------
#  cerwList <- n.sim.cond.3d(n.sim = 100, sim.locs, start=start, end=end,a0 = a0, g0 = g0,
#                        densities=P, qProbs=Q, DEM = dem, parallel = 4)

