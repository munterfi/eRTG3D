#' Plot function to plot the tracks in 2d
#'
#' @param origTrack a data.frame with x,y,z coordinates
#' @param cerwList a list containing a data.frame with x,y,z coordinates or a data.frame
#' @param titleText string with title of the plot
#' @param DEM an object of type 'RasterLayer', needs overlapping extent with the lines
#'
#' @return Nothing, plots a 2D ggplot2 object.
#' @export
#'
#' @examples
#' plot2d(track)
plot2d <- function(origTrack, cerwList = NULL, titleText = character(1), DEM = NULL)
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
  if(multipleTrack){
    p <- p + ggplot2::geom_path(data = cerwList, ggplot2::aes(x = x, y = y, color = z, group=group), size = 0.7, alpha = 0.7)
  }
  p <- p + ggplot2::geom_path(data = origTrack[origTrack$group==1, ], ggplot2::aes(x = x, y = y, group = group), color="grey", size = 2) +
    ggplot2::geom_path(data = origTrack, ggplot2::aes(x = x, y = y, group = group, color = z), size = 1) +
    ggplot2::geom_point(data=origTrack[1, ], ggplot2::aes(x = x, y = y), size=3.5, shape=7, alpha = 1, color="black") +
    ggplot2::geom_point(data=origTrack[nrow(origTrack), ], ggplot2::aes(x = x, y = y), size=3.5, shape=13, alpha = 1, color="black") +
    ggplot2::labs(colour="Flight height")
  print(p)
}

#' Density plots of turn angle, lift angle and step length
#'
#' The function takes either one track or two tracks.
#' The second track can be a list of tracks (eg. the output of n.sim.cons.3d()),
#' Then the densities of turn angle, lift angle and step length of all the simulations is taken.
#' Additionally the autodifferences parameter can be set to true, then the densities of the autodifferences
#' in turn angle, lift angle and step length are visualized.
#'
#' @param track1 a data.frame with x,y,z coordinates
#' @param track2 a list containing a data.frame with x,y,z coordinates or a data.frame
#' @param autodifferences logical: Should the densities of the autodifferences in turn angle, lift angle and step length are visualized.
#' @param scaleDensities logical: Should densities be scaled between 0 and 1, then sum of the area under the curve is not 1 anymore!
#'
#' @return A ggplot2 object.
#' @export
#'
#' @examples
#' plot2d.densities(track)
plot2d.densities <- function(track1, track2 = NULL, autodifferences = FALSE, scaleDensities = FALSE)
{
  track1 <- track.properties.3d(track1)[2:nrow(track1), ]
  if(!autodifferences) {
    t1 <- track1$t; l1 <- track1$l; d1 <- track1$d3d;
  } else {
    diffT1 <- diff(track1$t); diffL1 <- diff(track1$l); diffD1 <- diff(track1$d3d);
  }
  if(is.null(track2)) {
    t2 <- l2 <- d2 <- diffT2 <- diffL2 <- diffD2 <- NULL
  } else {
    if(is.data.frame(track2)){
      track2 <- track.properties.3d(track2)[2:nrow(track2), ]
      if(!autodifferences) {
        t2 <- track2$t; l2 <- track2$l; d2 <- track2$d3d;
      } else {
        diffT2 <- diff(track2$t); diffL2 <- diff(track2$l); diffD2 <- diff(track2$d3d);}
    }
    if(!is.data.frame(track2) && is.list(track2)) {
      track2 <- filter.dead.ends(track2)
      track2 <- lapply(track2, function(x){track.properties.3d(x)[2:nrow(x), ]})
      if(autodifferences) {
        diffTrack2 <- lapply(track2, function(x){data.frame(diffT = diff(x$t), diffL = diff(x$l), diffD = diff(x$d3d))})
        diffTrack2 <- do.call("rbind", diffTrack2)
        diffT2 <- diffTrack2$diffT; diffL2 <- diffTrack2$diffL; diffD2 <- diffTrack2$diffD;
      } else {
        track2 <- do.call("rbind", track2)
        t2 <- track2$t; l2 <- track2$l; d2 <- track2$d3d;
      }
    }
  }
  if(autodifferences){
    suppressWarnings(plot2d.multiplot(
      .plot2d.density(diffT1, diffT2, titleText = "Turn angle – autodifferences", scaleDensity = scaleDensities),
      .plot2d.density(diffL1, diffL2, titleText = "Lift angle – autodifferences", scaleDensity = scaleDensities),
      .plot2d.density(diffD1, diffD2, titleText = "Step length – autodifferences", scaleDensity = scaleDensities),
      cols = 1))
  } else {
    suppressWarnings(plot2d.multiplot(
      .plot2d.density(t1, t2, titleText = "Turn angle", scaleDensity = scaleDensities),
      .plot2d.density(l1, l2, titleText = "Lift angle", scaleDensity = scaleDensities),
      .plot2d.density(d1, d2, titleText = "Step length", scaleDensity = scaleDensities),
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
#' plot2d.density(values)
.plot2d.density <- function(values1, values2=NULL, titleText=character(1), scaleDensity = FALSE)
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
#' plot2d.multiplot(p1, p2, p3)
plot2d.multiplot <- function(..., plotlist=NULL, cols=1, layout=NULL)
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
