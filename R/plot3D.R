#' Plot track(s) with a surface of a digital elevation model in three dimensions
#'
#' @param origTrack a list containing data.frames with x,y,z coordinates or a data.frame
#' @param simTrack a list containing data.frames with x,y,z coordinates or a data.frame
#' @param titleText string with title of the plot
#' @param DEM an object of type \code{RasterLayer}, needs overlapping extent with the line(s)
#' @param padding adds a pad to the 2-D space in percentage (by default set to 0.1)
#' @param timesHeight multiply the height scale by a scalar (by default set to 10)
#'
#' @return
#' Plots a plotly object
#' @export
#'
#' @examples
#' plot3d(niclas)
plot3d <- function(origTrack, simTrack = NULL, titleText = character(1), DEM = NULL, padding = 0.1, timesHeight = 10) {
  if (!is.list(origTrack) || (!is.list(simTrack) && !is.null(simTrack))) stop("Track input has to be of type list or data.frame.")
  if (is.list(origTrack) && is.data.frame(origTrack)) {
    origTrack <- list(origTrack)
  }
  if (is.list(simTrack) && is.data.frame(simTrack)) {
    simTrack <- list(simTrack)
  }
  if (!is.null(simTrack)) {
    simTrack <- simTrack[!unlist(lapply(simTrack, is.null))]
  }
  if (padding > 1 && padding < 0) stop("The variable 'padding' must be a value between 0 and 1.")
  extents <- do.call("rbind", lapply(X = append(origTrack, simTrack), FUN = function(track) {
    c(
      floor(min(track$x)), floor(max(track$x)) + 1,
      floor(min(track$y)), floor(max(track$y)) + 1,
      floor(min(track$z)), floor(max(track$z)) + 1
    )
  }))
  minX <- min(extents[, 1])
  maxX <- max(extents[, 2])
  minY <- min(extents[, 3])
  maxY <- max(extents[, 4])
  minZ <- min(extents[, 5])
  maxZ <- max(extents[, 6])
  dx <- maxX - minX
  dy <- maxY - minY
  dz <- maxZ - minZ
  # pad extent
  minX <- minX - round(dx * padding)
  maxX <- maxX + round(dx * padding)
  minY <- minY - round(dy * padding)
  maxY <- maxY + round(dy * padding)
  minZ <- minZ - round(dz * padding)
  maxZ <- maxZ + round(dz * padding)
  dx <- maxX - minX
  dy <- maxY - minY
  dz <- maxZ - minZ
  ratio <- dy / dx
  # Create minimum cube
  d <- max(c(dx, dy, dz))
  middle <- (minX + maxX) / 2
  rMinX <- middle - d / 2
  rMaxX <- middle + d / 2
  middle <- (minY + maxY) / 2
  rMinY <- middle - d / 2
  rMaxY <- middle + d / 2
  middle <- (minZ + maxZ) / 2
  rMinZ <- middle - (d * 1 / timesHeight) / 2
  rMaxZ <- middle + (d * 1 / timesHeight) / 2
  # Define z axis
  axz <- list(title = "z", autoscale = TRUE, range = c(min(minZ, rMinZ), max(maxZ, rMaxZ)))
  p <- plotly::plot_ly()
  for (i in seq_len(length(origTrack))) {
    p <- plotly::add_trace(p,
      data = origTrack[[i]][1:3], x = ~x, y = ~y, z = ~z,
      mode = "lines+markers", type = "scatter3d", name = "Original",
      line = list(color = "black", width = 3),
      marker = list(color = "black", size = 2, cmin = -20, cmax = 50),
      opacity = 0.9, showlegend = (i == 1)
    )
  }
  if (!is.null(simTrack)) {
    for (i in seq_len(length(simTrack))) {
      p <- plotly::add_trace(p,
        data = simTrack[[i]][1:3], x = ~x, y = ~y, z = ~z,
        mode = "lines+markers", type = "scatter3d", name = "Simulation",
        line = list(color = "rgb(176,196,222)", width = 3),
        marker = list(color = "rgb(176,196,222)", size = 2, cmin = -20, cmax = 50),
        opacity = if (i == 1) {
          0.9
        } else {
          0.7
        }, showlegend = (i == 1)
      )
    }
  }
  if (!is.null(DEM)) {
    if (!class(DEM) == "RasterLayer") stop("'DEM' is not of type 'RasterLayer'")
    DEM <- raster::crop(DEM, raster::extent(minX, maxX, minY, maxY))
    DEM <- raster::resample(DEM, raster::raster(
      ncol = min(1000, (maxX - minX)), nrow = min(floor(1000 * ratio), (maxY - minY)),
      xmn = minX, xmx = maxX, ymn = minY, ymx = maxY
    ))
    DEM <- raster::as.matrix(DEM)
    p <- plotly::add_surface(p,
      x = seq(minX, maxX, length.out = ncol(DEM)), y = seq(maxY, minY, length.out = nrow(DEM)),
      z = DEM, type = "surface", opacity = 1, colorscale = list(c(0, 1, 2, 3, 4, 5), grDevices::terrain.colors(6)),
      reversescale = TRUE, colorbar = list(title = "DEM")
    )
    minDEM <- min(DEM, na.rm = TRUE)
    axz <- list(title = "z", autoscale = TRUE, range = c(minDEM, max(maxZ, ((rMaxZ - rMinZ) + minDEM))))
  }
  # extend plot area to minimum square with scaled z axis and add title
  p <- plotly::layout(p,
    title = titleText,
    scene = list(
      xaxis = list(title = "x", autoscale = TRUE, range = c(rMinX, rMaxX)),
      yaxis = list(title = "y", autoscale = TRUE, range = c(rMinY, rMaxY)),
      zaxis = axz
    )
  )
  print(p)
}

#' Plot function to plot the 3-D tracks in 2-D plane
#'
#' @param origTrack a list containing data.frames with x,y,z coordinates or a data.frame
#' @param simTrack a list containing data.frames with x,y,z coordinates or a data.frame
#' @param titleText string with title of the plot
#' @param DEM an object of type \code{RasterLayer}, needs overlapping extent with the line(s)
#' @param BG an object of type \code{RasterLayer}, needs overlapping extent with the line(s)
#' @param padding adds a pad to the 2-D space in percentage (by default set to 0.1)
#' @param alpha a number between 0 and 1, to specify the transparency of the simulated line(s)
#' @param resolution number of pixels the rasters are downsampled to (by default set to 500 pixels)
#'
#' @return A ggplot2 object.
#' @export
#'
#' @examples
#' plot2d(niclas)
plot2d <- function(origTrack, simTrack = NULL, titleText = character(1), DEM = NULL, BG = NULL, padding = 0.1, alpha = 0.7, resolution = 500) {
  if (!is.list(origTrack) || (!is.list(simTrack) && !is.null(simTrack))) stop("Track input has to be of type list or data.frame.")
  if (is.list(origTrack) && is.data.frame(origTrack)) {
    origTrack <- list(origTrack)
  }
  if (is.list(simTrack) && is.data.frame(simTrack)) {
    simTrack <- list(simTrack)
  }
  if (!is.null(simTrack)) {
    simTrack <- simTrack[!unlist(lapply(simTrack, is.null))]
  }
  if (padding > 1 && padding < 0) stop("The variable 'padding' must be a value between 0 and 1.")
  extents <- do.call("rbind", lapply(X = append(origTrack, simTrack), FUN = function(track) {
    c(
      floor(min(track$x)), floor(max(track$x)) + 1,
      floor(min(track$y)), floor(max(track$y)) + 1
    )
  }))
  minX <- min(extents[, 1])
  maxX <- max(extents[, 2])
  minY <- min(extents[, 3])
  maxY <- max(extents[, 4])
  dx <- maxX - minX
  dy <- maxY - minY
  # pad extent
  minX <- minX - round(dx * padding)
  maxX <- maxX + round(dx * padding)
  minY <- minY - round(dy * padding)
  maxY <- maxY + round(dy * padding)
  # set up plot
  p <- ggplot2::ggplot() +
    ggplot2::theme_classic() +
    ggplot2::coord_fixed(ratio = 1) +
    ggplot2::xlab("Easting") +
    ggplot2::ylab("Northing") +
    ggplot2::ggtitle(titleText)
  # prepare and add DEM to the plot
  if (!is.null(DEM)) {
    if (!class(DEM) == "RasterLayer") stop("'DEM' is not of type 'RasterLayer'")
    ratio <- (maxY - minY) / (maxX - minX)
    CRSsave <- raster::projection(DEM)
    DEM <- raster::crop(DEM, raster::extent(minX, maxX, minY, maxY))
    DEM <- raster::resample(DEM, raster::raster(
      ncol = min(resolution, (maxX - minX)), nrow = min(floor(resolution * ratio), (maxY - minY)),
      xmn = minX, xmx = maxX, ymn = minY, ymx = maxY
    ))
    raster::projection(DEM) <- CRSsave
    # creating hillshading from DHM:
    terr <- raster::terrain(DEM, opt = c("slope", "aspect"))
    hs <- raster::hillShade(terr$slope, terr$aspect, angle = 70, direction = 270)
    # convert rasters to dataframes for plotting with ggplot
    hdf <- data.frame(raster::rasterToPoints(hs))
    colnames(hdf) <- c("X", "Y", "Hillshade")
    hdf$Hillshade <- 1 - hdf$Hillshade
    ddf <- data.frame(raster::rasterToPoints(DEM))
    colnames(ddf) <- c("X", "Y", "DEM")
    p <- p +
      ggplot2::geom_raster(data = ddf, ggplot2::aes_string("X", "Y", fill = "DEM"), interpolate = TRUE) +
      ggplot2::scale_fill_gradientn(name = "Altitude", colours = grDevices::terrain.colors(4, alpha = 1)) +
      ggplot2::guides(fill = ggplot2::guide_colorbar()) +
      ggplot2::geom_tile(data = hdf, ggplot2::aes_string("X", "Y", alpha = "Hillshade"), fill = "grey20") +
      ggplot2::scale_alpha(range = c(0, 0.6))
  }
  # prepare and add BG to the plot
  if (!is.null(BG) & is.null(DEM)) {
    if (!class(BG) == "RasterLayer") stop("'DEM' is not of type 'RasterLayer'")
    ratio <- (maxY - minY) / (maxX - minX)
    CRSsave <- raster::projection(BG)
    BG <- raster::crop(BG, raster::extent(minX, maxX, minY, maxY))
    BG <- raster::resample(BG, raster::raster(ncol = min(resolution, (maxX - minX)), nrow = min(floor(resolution * ratio), (maxY - minY)), xmn = minX, xmx = maxX, ymn = minY, ymx = maxY))
    raster::projection(BG) <- CRSsave
    # convert rasters to dataframes for plotting with ggplot
    BG <- data.frame(raster::rasterToPoints(BG))
    colnames(BG) <- c("X", "Y", "BG")
    p <- p +
      ggplot2::geom_raster(data = BG, ggplot2::aes_string("X", "Y", fill = "BG"), interpolate = TRUE) +
      ggplot2::scale_fill_gradientn(name = "Uplift", colours = grDevices::heat.colors(4, alpha = 1)) +
      ggplot2::guides(fill = ggplot2::guide_colorbar())
  }
  # prepare tracks and add to plot
  if (!is.null(simTrack)) {
    simTrack <- simTrack[!unlist(lapply(simTrack, is.null))]
    simTrack <- lapply(X = 1:length(simTrack), FUN = function(X) {
      cbind(simTrack[[X]][, 1:3], group = rep(as.character(X), nrow(simTrack[[X]])))
    })
    simTrack <- do.call("rbind", simTrack)
    p <- p + ggplot2::geom_path(data = simTrack, ggplot2::aes_string(x = "x", y = "y", color = "z", group = "group"), size = 0.7, alpha = alpha)
  }
  origTrack <- origTrack[!unlist(lapply(origTrack, is.null))]
  origTrack <- lapply(X = seq_len(length(origTrack)), FUN = function(X) {
    cbind(origTrack[[X]][, 1:3], group = rep(as.character(X), nrow(origTrack[[X]])))
  })
  origTrack <- do.call("rbind", origTrack)
  p <- p + ggplot2::geom_path(data = origTrack, ggplot2::aes_string(x = "x", y = "y", group = "group"), color = "white", size = 2.5) +
    ggplot2::geom_path(data = origTrack, ggplot2::aes_string(x = "x", y = "y", group = "group", color = "z"), size = 1.5) +
    ggplot2::geom_point(data = origTrack[1, ], ggplot2::aes_string(x = "x", y = "y"), size = 3.5, shape = 7, alpha = 1, color = "black") +
    ggplot2::geom_text(data = origTrack[1, ], ggplot2::aes_string(x = "x", y = "y"), label = "START", fontface = "bold", size = 3, vjust = -1.4, color = "white") +
    ggplot2::geom_text(data = origTrack[1, ], ggplot2::aes_string(x = "x", y = "y"), label = "START", size = 3, vjust = -1.4, color = "black") +
    ggplot2::geom_point(data = origTrack[nrow(origTrack), ], ggplot2::aes_string(x = "x", y = "y"), size = 3.5, shape = 13, alpha = 1, color = "black") +
    ggplot2::geom_text(data = origTrack[nrow(origTrack), ], ggplot2::aes_string(x = "x", y = "y"), label = "END", fontface = "bold", size = 3, vjust = -1.4, color = "white") +
    ggplot2::geom_text(data = origTrack[nrow(origTrack), ], ggplot2::aes_string(x = "x", y = "y"), label = "END", size = 3, vjust = -1.4, color = "black") +
    ggplot2::labs(colour = "Height")
  return(p)
}

#' Density plots of turn angle, lift angle and step length
#'
#' The function takes either one track or two tracks.
#' The second track can be a list of tracks (eg. the output of \link[eRTG3D]{n.sim.cond.3d}),
#' Then the densities of turn angle, lift angle and step length of all the simulations is taken.
#' Additionally the autodifferences parameter can be set to true, then the densities of the autodifferences
#' in turn angle, lift angle and step length are visualized.
#'
#' @param track1 a list containing a data.frame with x,y,z coordinates or a data.frame
#' @param track2 a list containing a data.frame with x,y,z coordinates or a data.frame
#' @param autodifferences logical: should the densities of the autodifferences in turn angle, lift angle and step length are visualized.
#' @param scaleDensities logical: should densities be scaled between 0 and 1, then sum of the area under the curve is not 1 anymore!
#'
#' @return A ggplot2 object.
#' @export
#'
#' @examples
#' plot3d.densities(niclas)
plot3d.densities <- function(track1, track2 = NULL, autodifferences = FALSE, scaleDensities = FALSE) {
  if (!is.list(track1) || !is.list(track1)) stop("Track input has to be of type list or data.frame.")
  if (is.list(track1) && is.data.frame(track1)) {
    track1 <- list(track1)
    track1 <- filter.dead.ends(track1)
  }
  if (is.list(track2) && is.data.frame(track2)) {
    track2 <- list(track2)
    track2 <- filter.dead.ends(track2)
  }
  track1 <- lapply(track1, function(x) {
    track.properties.3d(x)[2:nrow(x), ]
  })
  track2 <- lapply(track2, function(x) {
    track.properties.3d(x)[2:nrow(x), ]
  })
  if (autodifferences) {
    difftrack1 <- do.call("rbind", lapply(track1, function(x) {
      data.frame(diffT = diff(x$t), diffL = diff(x$l), diffD = diff(x$d))
    }))
    diffT1 <- difftrack1$diffT
    diffL1 <- difftrack1$diffL
    diffD1 <- difftrack1$diffD
    diffTrack2 <- do.call("rbind", lapply(track2, function(x) {
      data.frame(diffT = diff(x$t), diffL = diff(x$l), diffD = diff(x$d))
    }))
    diffT2 <- diffTrack2$diffT
    diffL2 <- diffTrack2$diffL
    diffD2 <- diffTrack2$diffD
    suppressWarnings(plot3d.multiplot(
      .plot3d.density(diffT1, diffT2, titleText = "Turn angle - autodifferences", xlab = expression(paste(Delta, "t")), ylab = expression(paste("P(", Delta, "t)")), scaleDensity = scaleDensities),
      .plot3d.density(diffL1, diffL2, titleText = "Lift angle - autodifferences", xlab = expression(paste(Delta, "l")), ylab = expression(paste("P(", Delta, "l)")), scaleDensity = scaleDensities),
      .plot3d.density(diffD1, diffD2, titleText = "Step length - autodifferences", xlab = expression(paste(Delta, "d")), ylab = expression(paste("P(", Delta, "d)")), scaleDensity = scaleDensities),
      cols = 1
    ))
  } else {
    track1 <- do.call("rbind", track1)
    t1 <- track1$t
    l1 <- track1$l
    d1 <- track1$d
    track2 <- do.call("rbind", track2)
    t2 <- track2$t
    l2 <- track2$l
    d2 <- track2$d
    suppressWarnings(plot3d.multiplot(
      .plot3d.density(t1, t2, titleText = "Turn angle", xlab = "t", ylab = "P(t)", scaleDensity = scaleDensities),
      .plot3d.density(l1, l2, titleText = "Lift angle", xlab = "l", ylab = "P(l)", scaleDensity = scaleDensities),
      .plot3d.density(d1, d2, titleText = "Step length", xlab = "d", ylab = "P(d)", scaleDensity = scaleDensities),
      cols = 1
    ))
  }
}

#' Density plot of one or two variables
#'
#' @param values1 numeric vector
#' @param values2 numeric vector or \code{NULL}
#' @param titleText character of the title
#' @param xlab character of x label
#' @param ylab character of y label
#' @param scaleDensity logical: should density be scaled between 0 and 1, then sum of the area under the curve is not 1 anymore!
#'
#' @return A ggplot2 object.
#' @export
#'
#' @examples
#' .plot3d.density(niclas$x)
#' @noRd
.plot3d.density <- function(values1, values2 = NULL, titleText = character(1), xlab = "x", ylab = "P(x)", scaleDensity = FALSE) {
  values1 <- stats::na.omit(values1)
  if (!is.null(values2)) {
    values2 <- stats::na.omit(values2)
    dat <- data.frame(
      values = c(values1, values2),
      PDF = c(
        rep("Observed", times = length(values1)),
        rep("Simulated", times = length(values2))
      )
    )
  } else {
    dat <- data.frame(
      values = c(values1),
      PDF = c(rep("Trajectory", times = length(values1)))
    )
  }
  return(ggplot2::ggplot(dat, ggplot2::aes_string(x = "values", fill = "PDF")) +
    ggplot2::theme_classic() +
    (if (scaleDensity) {
      ggplot2::geom_density(ggplot2::aes_string(x = "values", y = "..scaled..", fill = "PDF"), alpha = 1 / 2)
    } else {
      ggplot2::geom_density(alpha = 0.5)
    }) +
    ggplot2::labs(x = xlab, y = ylab) +
    ggplot2::ggtitle(titleText))
}

#' Multiple plot function for ggplot objects
#'
#' If the layout is something like \code{matrix(c(1,2,3,3), nrow=2, byrow=TRUE)},
#' then plot 1 will go in the upper left, 2 will go in the upper right, and
#' 3 will go all the way across the bottom.
#'
#' @param ... ggplot objects
#' @param plotlist a list of ggplot objects
#' @param cols number of columns in layout
#' @param layout a matrix specifying the layout. If present, \code{cols} is ignored.
#'
#' @return Nothing, plots the ggplot2 objects.
#' @export
#'
#' @examples
#' plot3d.multiplot(plot2d(niclas), plot2d(niclas), plot2d(niclas))
plot3d.multiplot <- function(..., plotlist = NULL, cols = 1, layout = NULL) {
  # make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots <- length(plots)
  # if layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots / cols)),
      ncol = cols, nrow = ceiling(numPlots / cols)
    )
  }
  if (numPlots == 1) {
    print(plots[[1]])
  } else {
    # set up the page
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow(layout), ncol(layout))))
    # make each plot, in the correct location
    for (i in 1:numPlots) {
      # get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = grid::viewport(
        layout.pos.row = matchidx$row,
        layout.pos.col = matchidx$col
      ))
    }
  }
}

#' Visualize turn-lift-step histogram
#'
#' Creates a three dimensional scatterplot of the possibles next steps,
#' based on the tldCube, which was extracted from a track.
#'
#' @param tldCube tldCube; the ouptut from \link[eRTG3D]{turnLiftStepHist} or \link[eRTG3D]{get.densities.3d}
#'
#' @return Plots a plotly object
#' @export
#'
#' @examples
#' P <- get.track.densities.3d(niclas)
#' suppressWarnings(plot3d.tldCube(P$tldCube))
plot3d.tldCube <- function(tldCube) {
  if (!all(names(tldCube) == c("values", "tRes", "lRes", "dRes"))) stop("Input must be a tldCube, the ouptut from 'turnLiftStepHist()' and 'get.densities.3d()'.")
  # get coordinates of the tldCube
  t <- tldCube$values$turn
  l <- pi / 2 + tldCube$values$lift
  d <- tldCube$values$step
  # convert the coordinates from step length turning angle dimension
  df <- data.frame(
    x = d * sin(l) * cos(t),
    y = d * sin(l) * sin(t),
    z = d * cos(l),
    prob = tldCube$values$prob
  ) # get probs for each combination
  # intial direction
  dist <- max(df$x) * 0.1
  dfAx <- data.frame(x = c(0, max(df$x)), y = c(0, 0), z = c(0, 0))
  dfAy <- data.frame(x = c(0, 0), y = c(-dist, dist), z = c(0, 0))
  dfAz <- data.frame(x = c(0, 0), y = c(0, 0), z = c(-dist, dist))
  # 3-D scatterplot
  p <- plotly::plot_ly()
  if (length(unique(df$prob)) == 1) {
    p <- plotly::add_markers(p,
      data = df, x = ~x, y = ~y, z = ~z,
      marker = list(
        cmin = 0, cmax = unique(df$prob), opacity = 0.9, colorbar = list(title = "Probability"),
        color = ~prob, colorscale = "Bluered", showscale = TRUE
      )
    )
  } else {
    p <- plotly::add_markers(p,
      data = df, x = ~x, y = ~y, z = ~z, size = ~prob,
      marker = list(
        opacity = 0.9, colorbar = list(title = "Probability"),
        color = ~prob, colorscale = "Bluered", showscale = TRUE
      )
    )
  }
  p <- plotly::add_trace(p,
    data = dfAx, x = ~x, y = ~y, z = ~z,
    mode = "lines+markers", type = "scatter3d", name = "DirectionX",
    line = list(color = "black", width = 3),
    marker = list(color = "black", size = 0),
    opacity = 0.9, showlegend = FALSE
  )
  p <- plotly::add_trace(p,
    data = dfAy, x = ~x, y = ~y, z = ~z,
    mode = "lines+markers", type = "scatter3d", name = "DirectionY",
    line = list(color = "black", width = 3),
    marker = list(color = "black", size = 0),
    opacity = 0.9, showlegend = FALSE
  )
  p <- plotly::add_trace(p,
    data = dfAz, x = ~x, y = ~y, z = ~z,
    mode = "lines+markers", type = "scatter3d", name = "DirectionZ",
    line = list(color = "black", width = 3),
    marker = list(color = "black", size = 0),
    opacity = 0.9, showlegend = FALSE
  )
  p <- plotly::layout(p, title = paste("Bin width - t:", round(tldCube$tRes, 3),
    ", l:", round(tldCube$lRes, 3),
    ", d:", round(tldCube$dRes, 3),
    sep = ""
  ))
  suppressWarnings(print(p))
}
