#' Plot 3D track(s) with a surface
#'
#' @param origTrack a list containing data.frames with x,y,z coordinates or a data.frame
#' @param simTrack a list containing data.frames with x,y,z coordinates or a data.frame
#' @param titleText string with title of the plot
#' @param DEM an object of type 'RasterLayer', needs overlapping extent with the line(s)
#' @param padding adds a pad to the 2-D space in percentage (by default set to 0.05)
#' @param distanceHeightRatio ratio of the height in respect to the ground distance (by default ratio is 10)
#'
#' @return
#' Plots a 3D plotly object
#' @export
#'
#' @examples
#' plot3d(track)
plot3d.new <- function(origTrack, simTrack = NULL, titleText = character(1), DEM = NULL, padding = 0.05, distanceHeightRatio = 10) {
  if (!is.list(origTrack) || (!is.list(simTrack) && !is.null(simTrack))) stop("Track input has to be of type list or data.frame.")
  if (is.list(origTrack) && is.data.frame(origTrack)) {origTrack <- list(origTrack)}
  if (is.list(simTrack) && is.data.frame(simTrack)) {simTrack <-list(simTrack)}
  if (padding > 1 && padding < 0) stop("The variable 'padding' must be a value between 0 and 1.")
  extents <- do.call("rbind", lapply(X = append(origTrack, simTrack), FUN = function(track){
    c(floor(min(track$x)), floor(max(track$x))+1, floor(min(track$y)), floor(max(track$y))+1)
  }))
  minX <- min(extents[,1]); maxX <- max(extents[,2]); minY <- min(extents[,3]); maxY <- max(extents[,4]);
  dx <- maxX-minX; dy <- maxY-minY
  # pad extent
  minX <- minX - round(dx*padding); maxX <- maxX + round(dx*padding);
  minY <- minY - round(dy*padding); maxY <- maxY + round(dy*padding);
  dx <- maxX-minX; dy <- maxY-minY
  ratio <- dy/dx
  # Create minimum square
  if(dy<dx){
    middle <- (minY+maxY)/2
    rMinY <- middle - dx/2; rMaxY <- middle + dx/2
    rMinX <- minX; rMaxX <- maxX
  } else {
    middle <- (minX+maxX)/2
    rMinX <- middle - dy/2; rMaxX <- middle + dy/2
    rMinY <- minY; rMaxY <- maxY
  }
  # Define z axis
  axz <- list(title = 'z', autoscale = TRUE)
  # extend plot area to minimum square
  extentTrack <- as.data.frame(cbind(c(rMinX, rMaxX),c(rMinY, rMaxY),c(0,0)))
  colnames(extentTrack) <- c("x","y","z")
  # add to plotly
  p <- plotly::plot_ly()
  p <- plotly::add_trace(p, data = extentTrack, x = ~x, y = ~y, z = ~z,
                         mode = "lines", type = "scatter3d", opacity = 0, showlegend = FALSE)
  
  for (i in 1:length(origTrack)){
    p <- plotly::add_trace(p, data = origTrack[[i]][1:3], x = ~x, y = ~y, z = ~z,
                           mode = "lines+markers", type = "scatter3d", name = "Original",
                           line = list(color = "black", width = 3),
                           marker = list(size = 2, cmin = -20, cmax = 50),
                           opacity = 0.9, showlegend = (i==1))
  }
  if(!is.null(simTrack)){
    for (i in 1:length(simTrack)){
      p <- plotly::add_trace(p, data = simTrack[[i]][1:3], x = ~x, y = ~y, z = ~z,
                             mode = "lines+markers", type = "scatter3d", name = "Simulation",
                             line = list(color = "rgb(176,196,222)", width = 3),
                             marker = list(size = 2, cmin = -20, cmax = 50),
                             opacity = if(i==1){0.9}else{0.7}, showlegend = (i==1))
    }
  }
  if (!is.null(DEM)) {
    if(!class(DEM)=="RasterLayer") stop("'DEM' is not of type 'RasterLayer'")
    DEM <- raster::crop(DEM, raster::extent(minX, maxX, minY, maxY))
    DEM <- raster::resample(DEM, raster::raster(ncol=min(1000, (maxX-minX)), nrow=min(floor(1000*ratio), (maxY-minY)),
                                                xmn=minX, xmx=maxX, ymn=minY, ymx=maxY))
    DEM <- raster::as.matrix(DEM)
    p <- plotly::add_surface(p,
                             x = seq(minX, maxX, length.out = ncol(DEM)), y = seq(maxY, minY, length.out = nrow(DEM)),
                             z = DEM, type = "surface", opacity=1, colorscale=list(c(0,1,2,3,4,5),terrain.colors(6)),
                             reversescale = TRUE, colorbar=list(title='DEM'))
    minDEM <- min(DEM, na.rm = TRUE)
    axz <- list(title = 'z', autoscale = TRUE, range = c(minDEM, ((max(dx,dy)/distanceHeightRatio)+minDEM)))
  }
  p <- plotly::layout(p, title = titleText,
                      scene = list(xaxis = list(title = 'x', autoscale = TRUE),
                                   yaxis = list(title = 'y', autoscale = TRUE),
                                   zaxis = axz))
  print(p)
}

#' Plot 3D track(s) with a surface
#'
#' @param origTrack a data.frame with x,y,z coordinates
#' @param cerwList a list containing a data.frame with x,y,z coordinates or a data.frame
#' @param titleText string with title of the plot
#' @param DEM an object of type 'RasterLayer', needs overlapping extent with the line(s)
#' @param maxHeight Maximum plot height, default 8000m
#'
#' @return
#' Plots a 2D ggplot2 object
#' @export
#'
#' @examples
#' plot3d(track)
plot3d <- function(origTrack, cerwList=NULL, titleText = character(1), DEM=NULL, maxHeight=8000) {
  origTrack <- origTrack[, 1:3]
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
  if (!is.null(DEM)) {
    if(!class(DEM)=="RasterLayer") stop("'DEM' is not of type 'RasterLayer'")
    minX <- min(floor(min(origTrack$x)), min(unlist(lapply(X = cerwList, FUN = function(track){floor(min(track$x))}))))
    maxX <- max(floor(max(origTrack$x)), max(unlist(lapply(X = cerwList, FUN = function(track){floor(max(track$x))}))))+1
    minY <- min(floor(min(origTrack$y)), min(unlist(lapply(X = cerwList, FUN = function(track){floor(min(track$y))}))))
    maxY <- max(floor(max(origTrack$y)), max(unlist(lapply(X = cerwList, FUN = function(track){floor(max(track$y))}))))+1
    ratio <- (maxY-minY)/(maxX-minX)
    DEM <- raster::crop(DEM, raster::extent(minX, maxX, minY, maxY))
    DEM <- raster::resample(DEM, raster::raster(ncol=min(1000, (maxX-minX)), nrow=min(floor(1000*ratio), (maxY-minY)), xmn=minX, xmx=maxX, ymn=minY, ymx=maxY))
    DEM <- raster::as.matrix(DEM)
    axz <- list(
      title = 'z',
      range = c(min(DEM), maxHeight),
      autoscale = TRUE
    )
  }
  p <- plotly::plot_ly()
  p <- plotly::layout(p, title = titleText)
  p <- plotly::add_trace(p, data = origTrack[, 1:3], x = ~x, y = ~y, z = ~z,
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
  if (!is.null(DEM)) {
    p <- plotly::add_surface(p, x = seq(minX, maxX, length.out = ncol(DEM)), y = seq(maxY, minY, length.out = nrow(DEM)), z = DEM, type = "surface",
                opacity=1, colorscale=list(c(0,1,2,3,4,5),terrain.colors(6)), reversescale = TRUE,
                colorbar=list(
                  title='DEM'
                ))
    p <- plotly::layout(p, title = titleText,
           scene = list(xaxis = list(title = 'x', autoscale = FALSE),
                        yaxis = list(title = 'y', autoscale = FALSE),
                        zaxis = axz, bgcolor = "rgb(255, 255, 255)"))}
  print(p)
}

#' Plot function to plot the 3d tracks in 2d plane
#'
#' @param origTrack a data.frame with x,y,z coordinates
#' @param cerwList a list containing a data.frame with x,y,z coordinates or a data.frame
#' @param titleText string with title of the plot
#' @param DEM an object of type 'RasterLayer', needs overlapping extent with the line(s)
#' @param alpha a number between 0 and 1, to specify the transparency of the cerw line(s)
#'
#' @return A 2D ggplot2 object.
#' @export
#'
#' @examples
#' plot3d(track)
plot2d <- function(origTrack, cerwList = NULL, titleText = character(1), DEM = NULL, BG = NULL, alpha = 0.7)
{
  origTrack <- cbind(origTrack[ ,1:3], group = rep(as.character(1), nrow(origTrack)))
  multipleTrack <- FALSE
  if (is.data.frame(cerwList)) {
    origTrack <- rbind(origTrack, cbind(cerwList[ ,1:3], group = rep(as.character(2), nrow(cerwList))))
  }
  if (!is.data.frame(cerwList) && !is.null(cerwList)){
    cerwList <- cerwList[!unlist(lapply(cerwList, is.null))]
    cerwList <- lapply(X = 1:length(cerwList), FUN = function(X){cbind(cerwList[[X]], group = rep(as.character(X), nrow(cerwList[[X]])))})
    cerwList <- do.call("rbind", cerwList)
    multipleTrack <- TRUE
  }
  if (!is.null(DEM)) {
    if(!class(DEM)=="RasterLayer") stop("'DEM' is not of type 'RasterLayer'")
    if (is.null(cerwList)){
      cerwList <- origTrack
    }
    CRSsave <- raster::projection(DEM)
    minX <- min(floor(min(origTrack$x)), floor(min(cerwList$x)))
    maxX <- max(floor(max(origTrack$x)), floor(max(cerwList$x)))+1
    minY <- min(floor(min(origTrack$y)), floor(min(cerwList$y)))
    maxY <- max(floor(max(origTrack$y)), floor(max(cerwList$y)))+1
    ratio <- (maxY-minY)/(maxX-minX)
    DEM <- raster::crop(DEM, raster::extent(minX, maxX, minY, maxY))
    DEM <- raster::resample(DEM, raster::raster(ncol=min(750, (maxX-minX)), nrow=min(floor(750*ratio), (maxY-minY)), xmn=minX, xmx=maxX, ymn=minY, ymx=maxY))
    raster::projection(DEM) <- CRSsave
    # creating hillshading from DHM:
    terr = raster::terrain(DEM, opt=c("slope", "aspect"))
    hs = raster::hillShade(terr$slope, terr$aspect, angle=70, direction=270)
    # convert rasters to dataframes for plotting with ggplot
    hdf <- data.frame(raster::rasterToPoints(hs));
    colnames(hdf) <- c("X","Y","Hillshade")
    hdf$Hillshade <- 1-hdf$Hillshade
    ddf <- data.frame(raster::rasterToPoints(DEM));
    colnames(ddf) <- c("X","Y","DEM")
  }
  if (!is.null(BG) & is.null(DEM)) {
    if(!class(BG)=="RasterLayer") stop("'BG' is not of type 'RasterLayer'")
    if (is.null(cerwList)){
      cerwList <- origTrack
    }
    CRSsave <- raster::projection(BG)
    minX <- min(floor(min(origTrack$x)), floor(min(cerwList$x)))
    maxX <- max(floor(max(origTrack$x)), floor(max(cerwList$x)))+1
    minY <- min(floor(min(origTrack$y)), floor(min(cerwList$y)))
    maxY <- max(floor(max(origTrack$y)), floor(max(cerwList$y)))+1
    ratio <- (maxY-minY)/(maxX-minX)
    BG <- raster::crop(BG, raster::extent(minX, maxX, minY, maxY))
    BG <- raster::resample(BG, raster::raster(ncol=min(750, (maxX-minX)), nrow=min(floor(750*ratio), (maxY-minY)), xmn=minX, xmx=maxX, ymn=minY, ymx=maxY))
    raster::projection(BG) <- CRSsave
    # convert rasters to dataframes for plotting with ggplot
    BG <- data.frame(raster::rasterToPoints(BG));
    colnames(BG) <- c("X","Y","BG")
  }
  # Plot
  p <- ggplot2::ggplot() +
    ggplot2::theme_classic() +
    ggplot2::coord_fixed(ratio = 1) +
    ggplot2::xlab("Easting") +
    ggplot2::ylab("Northing") +
    ggplot2::ggtitle(titleText)
  if (!is.null(DEM)) {
    p <- p +
      ggplot2::geom_raster(data=ddf, ggplot2::aes(X,Y,fill=DEM), interpolate=TRUE) +
      ggplot2::scale_fill_gradientn(name="Altitude", colours = terrain.colors(4, alpha = 1)) +
      ggplot2::guides(fill = ggplot2::guide_colorbar()) +
      ggplot2::geom_tile(data=hdf, ggplot2::aes(X,Y,alpha=Hillshade), fill = "grey20") +
      ggplot2::scale_alpha(range = c(0, 0.6))
  }
  if (!is.null(BG) & is.null(DEM)) {
    p <- p +
      ggplot2::geom_raster(data=BG, ggplot2::aes(X,Y,fill=BG), interpolate=TRUE) +
      ggplot2::scale_fill_gradientn(name="Thermal Prob", colours = heat.colors(4, alpha = 1)) +
      ggplot2::guides(fill = ggplot2::guide_colorbar())
  }
  if(multipleTrack){
    p <- p + ggplot2::geom_path(data = cerwList, ggplot2::aes(x = x, y = y, color = z, group=group), size = 0.7, alpha = alpha)
  }
  p <- p + ggplot2::geom_path(data = origTrack[origTrack$group==1, ], ggplot2::aes(x = x, y = y, group = group), color="grey", size = 2) +
    ggplot2::geom_path(data = origTrack, ggplot2::aes(x = x, y = y, group = group, color = z), size = 1) +
    ggplot2::geom_point(data=origTrack[1, ], ggplot2::aes(x = x, y = y), size=3.5, shape=7, alpha = 1, color="black") +
    ggplot2::geom_point(data=origTrack[nrow(origTrack), ], ggplot2::aes(x = x, y = y), size=3.5, shape=13, alpha = 1, color="black") +
    ggplot2::labs(colour="Flight height")
  return(p)
}

#' Density plots of turn angle, lift angle and step length
#'
#' The function takes either one track or two tracks.
#' The second track can be a list of tracks (eg. the output of n.sim.cons.3d()),
#' Then the densities of turn angle, lift angle and step length of all the simulations is taken.
#' Additionally the autodifferences parameter can be set to true, then the densities of the autodifferences
#' in turn angle, lift angle and step length are visualized.
#'
#' @param track1 a list containing a data.frame with x,y,z coordinates or a data.frame
#' @param track2 a list containing a data.frame with x,y,z coordinates or a data.frame
#' @param autodifferences logical: Should the densities of the autodifferences in turn angle, lift angle and step length are visualized.
#' @param scaleDensities logical: Should densities be scaled between 0 and 1, then sum of the area under the curve is not 1 anymore!
#'
#' @return A ggplot2 object.
#' @export
#'
#' @examples
#' plot3d.densities(track)
plot3d.densities <- function(track1, track2 = NULL, autodifferences = FALSE, scaleDensities = FALSE)
{
  if (!is.list(track1) || !is.list(track1)) stop("Track input has to be of type list or data.frame.")
  if (is.list(track1) && is.data.frame(track1)) {track1 <- list(track1)}
  if (is.list(track2) && is.data.frame(track2)) {track2 <-list(track2)}
  track1 <- filter.dead.ends(track1); track2 <- filter.dead.ends(track2)
  track1 <- lapply(track1, function(x){track.properties.3d(x)[2:nrow(x), ]})
  track2 <- lapply(track2, function(x){track.properties.3d(x)[2:nrow(x), ]})
  if(autodifferences){
    difftrack1 <- do.call("rbind", lapply(track1, function(x){data.frame(diffT = diff(x$t), diffL = diff(x$l), diffD = diff(x$d))}))
    diffT1 <- difftrack1$diffT; diffL1 <- difftrack1$diffL; diffD1 <- difftrack1$diffD;
    diffTrack2 <- do.call("rbind", lapply(track2, function(x){data.frame(diffT = diff(x$t), diffL = diff(x$l), diffD = diff(x$d))}))
    diffT2 <- diffTrack2$diffT; diffL2 <- diffTrack2$diffL; diffD2 <- diffTrack2$diffD;
    suppressWarnings(plot3d.multiplot(
      .plot3d.density(diffT1, diffT2, titleText = "Turn angle – autodifferences", scaleDensity = scaleDensities),
      .plot3d.density(diffL1, diffL2, titleText = "Lift angle – autodifferences", scaleDensity = scaleDensities),
      .plot3d.density(diffD1, diffD2, titleText = "Step length – autodifferences", scaleDensity = scaleDensities),
      cols = 1))
  } else {
    track1 <- do.call("rbind", track1)
    t1 <- track1$t; l1 <- track1$l; d1 <- track1$d;
    track2 <- do.call("rbind", track2)
    t2 <- track2$t; l2 <- track2$l; d2 <- track2$d;
    suppressWarnings(plot3d.multiplot(
      .plot3d.density(t1, t2, titleText = "Turn angle", scaleDensity = scaleDensities),
      .plot3d.density(l1, l2, titleText = "Lift angle", scaleDensity = scaleDensities),
      .plot3d.density(d1, d2, titleText = "Step length", scaleDensity = scaleDensities),
      cols = 1))
  }
}

#' Density plot of one or two variables
#'
#' @param values1 numeric vector
#' @param values2 numeric vector or NULL
#' @param titleText character of the title text
#' @param scaleDensity logical: should density be scaled between 0 and 1, then sum of the area under the curve is not 1 anymore!
#'
#' @return A ggplot2 object.
#' @export
#'
#' @examples
#' plot3d.density(values)
.plot3d.density <- function(values1, values2=NULL, titleText=character(1), scaleDensity = FALSE)
{
  values1 <- na.omit(values1);
  if(!is.null(values2)){
    values2 <- na.omit(values2)
    dat <- data.frame(values = c(values1, values2),
                      Track = c(rep("Observed", times = length(values1)),
                                rep("Simulated", times = length(values2))))
  } else {
    dat <- data.frame(values = c(values1),
                      Track = c(rep("Track 1", times = length(values1))))
  }
  return(ggplot2::ggplot(dat, ggplot2::aes(x = values, fill = Track)) +
           #ggplot2::theme_classic() +
           (if(scaleDensity){ggplot2::geom_density(ggplot2::aes(x=values, y=..scaled.., fill=Track), alpha=1/2)}
            else{ggplot2::geom_density(alpha = 0.5)}) +
           ggplot2::labs(x = "Values", y = "Density") +
           ggplot2::ggtitle(titleText))
}

#' Multiple plot function for ggplot objects
#'
#' If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
#' then plot 1 will go in the upper left, 2 will go in the upper right, and
#' 3 will go all the way across the bottom.
#'
#' @param ... ggplot objects
#' @param plotlist a list of ggplot objects
#' @param cols number of columns in layout
#' @param layout a matrix specifying the layout. If present, 'cols' is ignored.
#'
#' @return Nothing, plots the ggplot2 objects.
#' @export
#'
#' @examples
#' plot3d.multiplot(p1, p2, p3)
plot3d.multiplot <- function(..., plotlist=NULL, cols=1, layout=NULL)
{
  # make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # if layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    # set up the page
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow(layout), ncol(layout))))
    # make each plot, in the correct location
    for (i in 1:numPlots) {
      # get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = grid::viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
    }
  }
}
