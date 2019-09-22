#' Apply voxel counting on a point cloud
#' 
#' A \code{rasterStack} object is created, representing the 3–D voxel cube.
#' The z axis is sliced into regular sections between the maximum and minimum value.
#' For every height slice a raster with points per cell counts is created. Additionally
#' the voxels can be standartized between 0 and 1.
#'
#' @param points a x, y, z data.frame
#' @param extent a raster extent object of the extent to create the rasters
#' @param xyRes resolution in the ground plane of the created rasters
#' @param zRes resolution in the z axis (by default \code{zRes = xyRes})
#' @param zMin minimum z value
#' @param zMax maximum height value
#' @param standartize logical: standartize the values?
#' @param verbose logical: print currently processed height band in raster stack?
#'
#' @return A \code{rasterStack} object, representing the 3–D voxel cube.
#' @export
#'
#' @examples
#' voxelCount(niclas, raster::extent(dem), 100, 100, 1000, 1400, standartize = TRUE)
voxelCount <- function(points, extent, xyRes, zRes = xyRes, zMin, zMax, standartize = FALSE, verbose = FALSE){
  rTem <- raster::raster(extent, res=xyRes)
  rTem[] <- 0
  rStack <- raster::stack()
  for (i in 1:round((zMax-zMin)/zRes)) {
    if (verbose) {
      cat('\r', paste("  |Counting points in Voxels for height: ", zMin+(i-1)*zRes, "m - ", (zMin+i*zRes), "m ...", sep = ""))
      utils::flush.console()
    }
    p <- points[points[,3] > (zMin+(i-1)*zRes) & points[,3] < (zMin+i*zRes), ]
    if (!nrow(p) == 0) {
      r <- raster::rasterize(cbind(p[,1], p[,2]), rTem, fun='count')
      r[is.na(r[])] <- 0
      rStack <- raster::stack(rStack, r)
    } else {
      rStack <- raster::stack(rStack, rTem)
    }
    names(rStack)[i] <- c(paste("m", zMin+(i-1)*zRes, "-", (zMin+i*zRes), sep = ""))
  }
  if(standartize){
    maxR <- max(raster::maxValue(rStack))
    minR <- min(raster::minValue(rStack))
    for(i in 1:length(rStack@layers)) {
      rStack@layers[[i]] <- (rStack@layers[[i]] - minR) / (maxR - minR)
    }
  }
  if (verbose) {
    cat('\r', "  |Done.                                                             \n")
    utils::flush.console()
  }
  return(rStack)
}

#' Chi maps of two variables
#' 
#' Calculates the chi maps for one \code{rasterStack} or all raster all the raster pairs stored in two \code{rasterStack}s.
#' As observed values, the first stack is used. The expected value is either set to the mean of the first stack, or if given 
#' to be the values of the second stack.
#' 
#' @param stack1 \code{rasterStack}
#' @param stack2 \code{rasterStack}  \code{NULL} or containing the same number of \code{rasterLayer}s and has euqal extent and resolution.
#' @param verbose logical: print currently processed height band in raster stack?
#'
#' @return A \code{rasterStack} containing the chi maps.
#' @export
#'
#' @examples
#' chiMaps(raster::stack(dem))
chiMaps <- function(stack1, stack2 = NULL, verbose = FALSE) {
  if (is.null(stack2)){
    rStack <- raster::stack()
    expMean <- mean(raster::values(stack1))
    for (i in 1:length(stack1@layers)){
      if (verbose) {
        cat('\r', paste("  |Calcuate chi map for raster:", i, "..."))
        utils::flush.console()
      }
      r1 <- stack1@layers[[i]]; r1[r1 == 0] <- NA
      rChi <- (r1 - expMean) / sqrt(expMean)
      rChi[is.na(rChi)] <- 0
      rStack <- raster::stack(rStack, rChi)
      names(rStack)[i] <- c(paste("chiMap",names(stack1)[i], sep = ))
    }
    if (verbose) {
      cat('\r', "  |Done.                                                             \n")
      utils::flush.console()
    }
    return(rStack)
  } else {
    if(length(stack1@layers) != length(stack2@layers)) stop("Stack need to have the same number of rasters")
    rStack <- raster::stack()
    for (i in 1:length(stack1@layers)){
      if (verbose) {
        cat('\r', paste("  |Calcuate chi map for raster:", i, "..."))
        utils::flush.console()
      }
      r1 <- stack1@layers[[i]]; r1[r1 == 0] <- NA
      r2 <- stack2@layers[[i]]; r2[r2 == 0] <- NA
      rChi <- (r1 - r2) / sqrt(r2)
      rChi[is.na(rChi)] <- 0
      rStack <- raster::stack(rStack, rChi)
      names(rStack)[i] <- c(paste("chiMap",names(stack1)[i], sep = ))
    }
    if (verbose) {
      cat('\r', "  |Done.                                                             \n")
      utils::flush.console()
    }
    return(rStack)
  }
}

#' Plots a rasterLayer or rasterStack
#'
#' @param r \code{rasterLayer} or \code{rasterStack}
#' @param title title text of plot(s)
#' @param centerColorBar logical: center colobar around 0 and use \code{RdBuTheme()}?
#' @param ncol number of columns to plot a stack, by default estimated by the square root
#'
#' @return Plots the rasters
#' @export
#'
#' @examples
#' plotRaster(dem)
plotRaster <- function(r, title = character(0), centerColorBar = FALSE, ncol = NULL){
  if (centerColorBar) {
    colTheme <- rasterVis::BuRdTheme()
    maxVal <- max(raster::values(r)[is.finite(raster::values(r))], na.rm = TRUE)
    colSeq <- seq(-1*maxVal, maxVal, len=100)
  } else {
    colTheme <- rasterVis::YlOrRdTheme()
    minVal <- min(raster::values(r)[is.finite(raster::values(r))], na.rm = TRUE)
    maxVal <- max(raster::values(r)[is.finite(raster::values(r))], na.rm = TRUE)
    colSeq <- seq(minVal, maxVal, len=100)
  }
  if (class(r)[1] == "RasterLayer"){
    print(rasterVis::levelplot(r, par.settings=colTheme, interpolate = TRUE,
                               margin=FALSE, at=colSeq,
                               main=title, xlab="Easting", ylab="Northing"))
  }
  if (class(r)[1] == "RasterBrick"){r <- raster::stack(r)}
  if (class(r)[1] == "RasterStack"){
    plotList <- list()
    for (i in 1:length(r@layers)){
      if(raster::cellStats(r@layers[[i]], "sum") != 0) {
        p <- rasterVis::levelplot(r@layers[[i]], par.settings=colTheme, interpolate = TRUE,
                                  margin=FALSE, at=colSeq,
                                  main = names(r)[i], xlab="Easting", ylab="Northing")
        plotList <- append(plotList, list(p))
      }
    }
    if(is.null(ncol)){ncol <- round(sqrt(length(r@layers)))}
    do.call(eval(parse(text="gridExtra::grid.arrange")), c(plotList, ncol = ncol))
  }
}

#' Converts a rasterStack to logarithmic scale
#' 
#' Avoids the problem of -Inf occuring for log(0).
#'
#' @param rStack rasterStack to convert to logarithmic scale
#' @param standartize logical: standartize cube between 0 and 1
#' @param InfVal the value that \code{Inf} and \code{-Inf} should be rpeplaced with
#'
#' @return A rasterStack in logarithmic scale
#' @export
#'
#' @examples
#' logRasterStack(raster::stack(dem))
logRasterStack <- function(rStack, standartize = FALSE, InfVal = NA)
{
  rStack <- log(rStack)
  rStack[is.infinite(rStack)] <- InfVal
  if(standartize) {
    naInd <- is.na(rStack)
    rStack[naInd] <- 0
    maxR <- max(raster::maxValue(rStack))
    minR <- min(raster::minValue(rStack))
    rStack <- raster::stack(rStack)
    for(i in 1:length(rStack@layers)) {
      rStack@layers[[i]] <- (rStack@layers[[i]] - minR) / (maxR - minR)
    }
    rStack[naInd] <- NA
  }
  return(raster::stack(rStack))
}

#' Export a dataCube as image slice sequence
#' 
#' Exports a dataCube of type \code{rasterStack} as Tiff image sequence.
#' Image sequences are a common structure to represent voxel data and
#' most of the specific software to visualize voxel data is able to read it (e.g. blender)
#'
#' @param rStack rasterStack to be saved to Tiff image slices
#' @param filename name of the image slices
#' @param dir directory, where the slices should be stored
#' @param NaVal numeric value that should represent NA values in the Tiff image, default is \code{NaVal = 0}
#'
#' @return Saves the Tiff image files.
#' @export
#'
#' @examples
#' crws <- lapply(X=seq(1:100), FUN = function(X) {
#'   sim.crw.3d(nStep = 100, rTurn = 0.99, rLift = 0.99, meanStep = 0.1)
#' })
#' points <- do.call("rbind", crws)
#' extent <- raster::extent(c(-10, 10, -10, 10))
#' ud <- voxelCount(points, extent, xyRes=5,
#'                  zMin=-10, zMax=10, standartize = TRUE)
#' saveImageSlices(ud, filename = "saveImageSlices_test", dir = tempdir())
saveImageSlices <- function(rStack, filename, dir, NaVal = 0)
{
  rStack[is.na(rStack)] <- NaVal
  rStack <- raster::stack(rStack)
  for (i in 1:length(rStack@layers)) {
    tiff::writeTIFF(raster::as.matrix(rStack[[i]]),
                    file.path(dir, paste(sprintf("%03d", i), ".", filename, ".tif", sep = "")))
  }
}
