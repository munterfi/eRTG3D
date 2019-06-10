## ---- echo = FALSE, eval=TRUE, include=FALSE-----------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
options(knitr.table.format = "html")

library(eRTG3D)
library(ggplot2)
set.seed(123)

## ------------------------------------------------------------------------
trajectory.3D <- sim.crw.3d(nStep = 100, rTurn = 0.99, rLift = 0.99, meanStep = 0.1)

## ------------------------------------------------------------------------
trajectory.2D <- trajectory.3D
trajectory.2D$z <- 0
head(trajectory.2D)

## ---- echo = FALSE, eval=TRUE, include=FALSE-----------------------------
simulation.2D <- reproduce.track.3d(trajectory.2D)

## ----eval=FALSE----------------------------------------------------------
#  simulation.2D <- reproduce.track.3d(trajectory.2D)

## ----eval=TRUE, fig.align="center"---------------------------------------
plot2d(trajectory.2D, simulation.2D, titleText = "2-D trajectory simulation")

