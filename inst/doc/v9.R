## ---- echo = FALSE, eval=TRUE, include=FALSE-----------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
options(knitr.table.format = "html")

library(eRTG3D)
library(raster)
library(ggplot2)
set.seed(123)

## ----eval=TRUE, echo=TRUE, fig.align="center"----------------------------
crws <- lapply(X=seq(1:100), FUN = function(X) {
  sim.crw.3d(nStep = 100, rTurn = 0.99, rLift = 0.99, meanStep = 0.1)
  })
plot2d(crws)

## ----eval=TRUE, echo=TRUE, fig.align="center"----------------------------
points <- do.call("rbind", crws)
extent <- extent(c(-10, 10, -10, 10))
ud <- voxelCount(points, extent, xyRes=5, zMin=-10, zMax=10)
plotRaster(ud)

## ----eval=TRUE, echo=TRUE, fig.align="center"----------------------------
chi <- chiMaps(ud)
plotRaster(chi, centerColorBar=TRUE)

## ----eval = FALSE--------------------------------------------------------
#  saveImageSlices(ud, filename = "utilization-distribution", dir="folder/path")
#  saveImageSlices(chi, filename = "chi-map-cube", dir="folder/path")

