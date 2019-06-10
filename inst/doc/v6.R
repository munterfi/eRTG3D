## ---- echo = FALSE, eval=TRUE, include=FALSE-----------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
options(knitr.table.format = "html")

library(eRTG3D)
set.seed(123)
WGS84 <- "+init=epsg:4326"
LV95 <- "+init=epsg:2056"
track.wgs84 <- transformCRS.3d(niclas, fromCRS=LV95, toCRS=WGS84)[ ,1:3]

## ----eval=TRUE-----------------------------------------------------------
head(track.wgs84)

## ----eval=TRUE-----------------------------------------------------------
WGS84 <- "+init=epsg:4326"
LV95 <- "+init=epsg:2056"
track <- transformCRS.3d(track.wgs84, fromCRS=WGS84, toCRS=LV95)
head(track)

