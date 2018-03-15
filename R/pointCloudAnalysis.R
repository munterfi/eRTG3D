voxelcount <- function(points, extent, xyRes, zRes = xyRes, zMin, zMax, standartize = FALSE){
  rTem <- raster::raster(extent, res=xyRes)
  rTem[] <- 0
  rStack <- raster::stack()
  for (i in 1:round((zMax-zMin)/zRes)) {
    cat('\r', paste("  |Counting points in Voxels for height: ", zMin+(i-1)*zRes, "m - ", (zMin+i*zRes), "m ...", sep = ""))
    flush.console()
    p <- points[points[,3] > (zMin+(i-1)*zRes) & points[,3] < (zMin+i*zRes), ]
    if (!nrow(p) == 0) {
      p <- sp::SpatialPoints(coords = cbind(p[,1], p[,2]))
      r <- raster::rasterize(p, rTem, fun='count')
      r[is.na(r[])] <- 0
      if(standartize){
        r <- ((r - cellStats(r, "min")) / (cellStats(r, "max") - cellStats(r, "min")))
      }
      rStack <- raster::stack(rStack, r)
    } else {
      rStack <- raster::stack(rStack, rTem)
    }
    names(rStack)[i] <- c(paste("m", zMin+(i-1)*zRes, "-", (zMin+i*zRes), sep = ""))
  }
  cat('\r', "  |Done.                                                             ")
  flush.console()
  return(rStack)
}

chimaps <- function(stack1, stack2) {
  if(length(stack1@layers) != length(stack2@layers)) stop("Stack need to have the same number of rasters")
  rStack <- raster::stack()
  for (i in 1:length(stack1@layers)){
    cat('\r', paste("  |Calcuate chi-map for raster:", i, "..."))
    flush.console()
    r1 <- stack1@layers[[i]]; r1[r1 == 0] <- NA
    r2 <- stack2@layers[[i]]; r2[r2 == 0] <- NA
    rChi <- (r1 - r2) / sqrt(r2)
    rChi[is.na(rChi)] <- 0
    rStack <- raster::stack(rStack, rChi)
    names(rStack)[i] <- c(paste("chimap",names(stack1)[i], sep = ))
  }
  cat('\r', "  |Done.                                                             ")
  flush.console()
  return(rStack)
}

plotRaster <- function(r, title = character(0)){
  if (class(r)[1] == "RasterLayer"){
    print(rasterVis::levelplot(r, par.settings=rasterVis::BuRdTheme(), interpolate = TRUE,
                               margin=FALSE, at=seq(-1.01*maxValue(abs(r)), 1.01*maxValue(abs(r)), len=100),
                               main=title, xlab="Easting", ylab="Northing"))
  }
  if (class(r)[1] == "RasterStack"){
    plotList <- list()
    for (i in 1:length(r@layers)){
      if(cellStats(r@layers[[i]], "sum") != 0) {
        p <- rasterVis::levelplot(r@layers[[i]], par.settings=rasterVis::BuRdTheme(), interpolate = TRUE,
                                  margin=FALSE, at=seq(-1.01*maxValue(abs(r@layers[[i]])), 1.01*maxValue(abs(r@layers[[i]])), len=100),
                                  main = names(r)[i], xlab="Easting", ylab="Northing")
        plotList <- append(plotList, list(p))
      }
    }
    do.call(eval(parse(text="gridExtra::grid.arrange")), c(plotList, ncol = round(sqrt(length(r@layers)))))
  }
}