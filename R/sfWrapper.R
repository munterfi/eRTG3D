sf2df <- function(track)
{
  track <- cbind(sf::st_coordinates(track), as.data.frame(track))
  track$geometry <- NULL
  return(track)
}

track2sf.3d <- function(track, CRS = NA)
{
  if(is.data.frame(track)) {return(.df2sf.3d(track, CRS = CRS))}
  if(is.matrix(track)) {return(.matrix2sf.3d(track, CRS = CRS))}
}

.df2sf.3d <- function(track, CRS = NA)
{
  cols <- colnames(track)
  return(sf::st_as_sf(track, coords = cols[match(c("x", "y", "z"), cols)], crs = CRS))
}

.matrix2sf.3d <- function(track, CRS = NA)
{
  return(sf::st_as_sf(as.data.frame(track), coords = c(1,2,3), crs = CRS))
}
