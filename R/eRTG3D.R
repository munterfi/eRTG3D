#' eRTG3D: Empirically Informed Random Trajectory Generator in 3-D
#'
#' The empirically informed random trajectory generator in three dimensions (eRTG3D)
#' is an algorithm to generate realistic random trajectories in a 3-D space
#' between two given fix points in space, so-called Conditional Empirical Random Walks.
#' The trajectory generation is based on empirical distribution functions extracted from
#' observed trajectories (training data) and thus reflects the geometrical movement
#' characteristics of the mover. A digital elevation model (DEM),
#' representing the Earth's surface, and a background layer of probabilities
#' (e.g. food sources, uplift potential, waterbodies, etc.) can be used to
#' influence the trajectories.
#'
#' See the packages site on \href{https://munterfi.github.io/eRTG3D/}{GitHub},
#' detailed information about the algorithm in this \href{https://www.geo.uzh.ch/dam/jcr:6194e41e-055c-4635-9807-53c5a54a3be7/MasterThesis_Unterfinger_2018.pdf}{Master’s Thesis},
#' or test the algorithm online in the \href{https://mufi.shinyapps.io/ertg3d-simulator}{eRTG3D Simulator}.
#'
#' @docType package
#' @name eRTG3D
NULL

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))
