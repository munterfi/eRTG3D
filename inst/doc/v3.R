## ---- echo = FALSE, eval=TRUE, include=FALSE-----------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
options(knitr.table.format = "html")

library(eRTG3D)
library(ggplot2)
set.seed(1234)
cerwList <- reproduce.track.3d(n.sim = 5, niclas, DEM = dem)
# <img src="figs/plot3d.png" alt="Drawing" style="width: 96.5%;"/>

## ----eval=FALSE----------------------------------------------------------
#  cerwList <- reproduce.track.3d(n.sim = 5, niclas, DEM = dem)

## ----eval=TRUE, fig.align="center"---------------------------------------
plot2d(niclas, cerwList, DEM=dem,
       titleText=paste("Steps: ", nrow(niclas), ", Niclas", sep=""))

## ----eval=TRUE, echo=FALSE, fig.align="center"---------------------------
#p <- plot3d(crws)
#htmltools::tags$img(p, timesHeight=1)

## ----eval=FALSE----------------------------------------------------------
#  plot3d(niclas, cerwList, DEM=dem,
#         titleText=paste("Steps: ", nrow(niclas), ", Niclas", sep=""))

## ----eval=TRUE, fig.align="center"---------------------------------------
plot3d.densities(niclas, cerwList)

