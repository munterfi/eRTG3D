#' Transform coordinates reference system (CRS) of a 3-D track
#'
#' Attention: Please use this function for CRS transformations,
#' since it is based on the \link[sf]{st_transform} from the sf package and therefore
#' supports CRS transformations in 3-D. Note: \code{spTransform} from the \code{sp} package
#' only supports transformations in the 2D plane, which will cause distortions
#' in the third dimension.
#'
#' @param track data.frame with x,y,z coordinates
#' @param fromCRS numeric, EPSG code of the current CRS
#' @param toCRS numeric, EPSG code of the CRS to be converted in
#'
#' @return A data.frame containing x,y,z and variables.
#' @export
#'
#' @examples
#' transformCRS.3d(niclas, fromCRS = 2056, toCRS = 4326)
transformCRS.3d <- function(track, fromCRS, toCRS) {
  track <- track2sf.3d(track = track, CRS = fromCRS)
  track <- sf::st_transform(track, toCRS)
  track <- sf2df.3d(track)
  track <- track.properties.3d(track)
  return(track)
}

#' Tests if the object is a simple feature collection (class: \code{'sf, data.frame'})
#'
#' @param track any object to test
#'
#' @return A logical: \code{TRUE} if is a simple feature collection (class: \code{'sf, data.frame'}) of the sf package, \code{FALSE} otherwise.
#' @export
#'
#' @examples
#' is.sf.3d(niclas)
#' is.sf.3d(track2sf.3d(track = niclas, CRS = 2056))
is.sf.3d <- function(track) {
  if (inherits(track, "sf")) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' Converts a sf data.frame to a normal dataframe
#'
#' @param track An object of type \code{'sf, data.frame'}
#'
#' @return A data.frame.
#' @export
#'
#' @examples
#' sf2df.3d(track2sf.3d(niclas, CRS = 4326))
sf2df.3d <- function(track) {
  track <- cbind(sf::st_coordinates(track), as.data.frame(track))
  track$geometry <- NULL
  colnames(track)[1:3] <- c("x", "y", "z")
  return(track)
}

#' Converts a track to a \code{'sf, data.frame'}
#'
#' @param track eRTG3D track data.frame or a matrix
#' @param CRS numeric, EPSG code of the CRS
#'
#' @return A track of type \code{'sf, data.frame'}.
#' @export
#'
#' @examples
#' track2sf.3d(niclas, 4326)
track2sf.3d <- function(track, CRS = NA) {
  if (is.data.frame(track)) {
    return(.df2sf.3d(track, CRS = CRS))
  }
  if (is.matrix(track)) {
    return(.matrix2sf.3d(track, CRS = CRS))
  }
  if (inherits(track, "Move")) {
    return(.move2sf.3d(track))
  }
}

#' Converts a track data.frame to a \code{'sf, data.frame'}
#'
#' @param track eRTG3D track data.frame
#' @param CRS numeric, EPSG code of the CRS
#'
#' @return A track of type \code{'sf, data.frame'}.
#' @export
#'
#' @examples
#' .df2sf.3d(track, CRS = NA)
#' @noRd
.df2sf.3d <- function(track, CRS = NA) {
  if (any(is.na(track[, 1:3]))) stop("Track 'data.frame' contains NA values.")
  track <- track.properties.3d(track)
  return(sf::st_as_sf(track, coords = c(1, 2, 3), crs = CRS))
}

#' Converts a track matrix to a sf data.frame
#'
#' @param track matrix with x, y and z coordinates
#' @param CRS numeric, EPSG code of the CRS
#'
#' @return A track of type \code{'sf, data.frame'}.
#' @export
#'
#' @examples
#' .matrix2sf.3d(track, CRS = NA)
#' @noRd
.matrix2sf.3d <- function(track, CRS = NA) {
  if (any(is.na(track[, 1:3]))) stop("Track 'matrix' contains NA values.")
  track <- track.properties.3d(as.data.frame(track))
  return(sf::st_as_sf(track, coords = c(1, 2, 3), crs = CRS))
}

#' Converts a move object to a sf data.frame
#'
#' @param track move object with a single track
#'
#' @return A track of type \code{'sf, data.frame'}.
#' @export
#'
#' @examples
#' .move2sf.3d(track)
#' @noRd
.move2sf.3d <- function(track) {
  CRS <- as.character(track@proj4string)
  track <- data.frame(
    x = track$location_long,
    y = track$location_lat,
    z = track$height_above_ellipsoid
  )
  na.ind <- which(is.na(track$z))
  if (any(na.ind)) {
    warning("NA values in z: Filled by linear interpolation")
    track$z[na.ind] <- (track$z[na.ind + 1] + track$z[na.ind - 1]) / 2
  }
  track <- track.properties.3d(track)
  return(sf::st_as_sf(track, coords = c(1, 2, 3), crs = CRS))
}
