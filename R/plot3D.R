#' Plot 3D track(s) with a surface
#'
#' @param origTrack a data.frame with x,y,z coordinates
#' @param cerwList a list containing a data.frame with x,y,z coordinates or a data.frame
#' @param titleText string with title of the plot
#' @param surface logical: should the surface layer be plotted? If no surface raster is provided, a zero plane is created.
#' @param DEM an object of type 'RasterLayer', needs overlapping extent with the lines
#' @param maxHeight Maximum plot height, default 8000m
#'
#' @return
#' Plots a 2D ggplot2 object
#' @export
#'
#' @examples
#' plot3d(track)
plot3d <- function(origTrack, cerwList=NULL, titleText = character(1), surface=FALSE, DEM=NULL, maxHeight=8000) {
  multipleTrack <- TRUE
  singleTrack <- FALSE
  if (is.data.frame(cerwList)){
    multipleTrack <- FALSE
    cerwList <- list(cerwList)
  }
  if (is.null(cerwList)){
    cerwList <- list(origTrack)
    multipleTrack <- FALSE
    singleTrack <- TRUE
  }
  cerwList <- cerwList[!unlist(lapply(cerwList, is.null))]
  if (surface) {
    minX <- min(floor(min(origTrack$x)), min(unlist(lapply(X = cerwList, FUN = function(track){floor(min(track$x))}))))
    maxX <- max(floor(max(origTrack$x)), max(unlist(lapply(X = cerwList, FUN = function(track){floor(max(track$x))}))))+1
    minY <- min(floor(min(origTrack$y)), min(unlist(lapply(X = cerwList, FUN = function(track){floor(min(track$y))}))))
    maxY <- max(floor(max(origTrack$y)), max(unlist(lapply(X = cerwList, FUN = function(track){floor(max(track$y))}))))+1
    ratio <- (maxY-minY)/(maxX-minX)
    if (is.null(DEM)) {
      DEM <- raster::raster(ncol=min(1000, (maxX-minX)), nrow=min(floor(1000*ratio), (maxY-minY)), xmn=minX, xmx=maxX, ymn=minY, ymx=maxY)
      raster::values(DEM) <- runif(raster::ncell(DEM))
    } else {
      if(!class(DEM)=="RasterLayer") stop("'DEM' is not of type 'RasterLayer'")
      DEM <- raster::crop(DEM, raster::extent(minX, maxX, minY, maxY))
      DEM <- raster::resample(DEM, raster::raster(ncol=min(1000, (maxX-minX)), nrow=min(floor(1000*ratio), (maxY-minY)), xmn=minX, xmx=maxX, ymn=minY, ymx=maxY))
    }
    DEM <- raster::as.matrix(DEM)
    axz <- list(
      title = 'z',
      range = c(min(DEM), maxHeight),
      autoscale = TRUE
    )
  }
  p <- plotly::plot_ly()
  p <- plotly::layout(p, title = titleText)
  p <- plotly::add_trace(p, data = origTrack, x = ~x, y = ~y, z = ~z,
              mode = "lines+markers", type = "scatter3d", name = "Observed",
              line = list(color = "black", width = 3),
              marker = list(size = 2, cmin = -20, cmax = 50), opacity = 0.9, showlegend = (!singleTrack))
  p <- plotly::add_trace(p, data = cerwList[[1]], x = ~x, y = ~y, z = ~z,
              mode = "lines+markers", type = "scatter3d", name = "CERW",
              line = list(color = "blue", width = 3),
              marker = list(size = 2, cmin = -20, cmax = 50), opacity = 0.9, showlegend = (!singleTrack))
  if(multipleTrack){
    for (i in 2:length(cerwList)){
      p <- plotly::add_trace(p, data = cerwList[[i]], x = ~x, y = ~y, z = ~z,
                           mode = "lines+markers", type = "scatter3d", name = "CERW",
                           line = list(color = "rgb(176,196,222)", width = 3),
                           marker = list(size = 2, cmin = -20, cmax = 50),
                           opacity = 0.4,
                           showlegend = FALSE)
    }
  }
  if (surface) {
    p <- plotly::add_surface(p, x = seq(minX, maxX, length.out = ncol(DEM)), y = seq(maxY, minY, length.out = nrow(DEM)), z = DEM, type = "surface",
                opacity=1, colorscale=list(c(0,1,2,3,4,5),terrain.colors(6)), reversescale = TRUE,
                colorbar=list(
                  title='Surface'
                ))
    p <- plotly::layout(p, title = titleText,
           scene = list(xaxis = list(title = 'x', autoscale = FALSE),
                        yaxis = list(title = 'y', autoscale = FALSE),
                        zaxis = axz, bgcolor = "rgb(255, 255, 255)"))}
  print(p)
}
