## ---- echo = FALSE, eval=TRUE, include=FALSE-----------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
options(knitr.table.format = "html")

library(eRTG3D)
set.seed(123)

## ----eval=TRUE, fig.height=5, fig.width=7--------------------------------
LV95 <- "+init=epsg:2056"
niclas <- track2sf.3d(niclas, CRS = LV95)
is.sf.3d(niclas)
head(niclas, 3)

## ----eval=TRUE, fig.height=5, fig.width=7--------------------------------
niclas <- sf2df.3d(niclas)
is.sf.3d(niclas)

