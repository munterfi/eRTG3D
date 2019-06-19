# eRTG3D 0.6.0
Initial submission to CRAN. New github page created with `pkgdown`.
Added continuous integration by Travis CI, unit testing with `testthat` and `codecovr` for test coverage assessment. Restructured and updated package documentation. 

## New features
* Package now available on CRAN: `install.packages("eRTG3D")`
* GitHub page for the package. Created with `pkgdown`
* Documentation for conducting 2-D simulations and voxel counting in order to extract utilization distributions.
* Continuous integration by Travis CI
* Unit testing with `testthat`
* Test coverage assessment using `codecovr`
* New eRTG3D logo
* New index/home area with interactive 3-D plots
* Shiny app - Online eRTG3D simulator: https://mufi.shinyapps.io/ertg3d-simulator/

## R CMD Check
0 errors ✔ | 0 warnings ✔ | 0 notes ✔

R CMD check succeeded

# eRTG3D 0.5.9
The 2-D case should also be handled by a 3-D algorithm, therefore the handling of 2-D trajectories is now enabled. To use simulations in 2-D, just pass in a data.frame where the z dimension is zero or has a constant value. The algorithm handles 2-D input correctly and the resulting simulations are valid. Since the third dimension is still part of the 2-D trajectory (either 0 or a constant value), the combination of 3-D and 2-D simulations is straight forward. Furthermore, the wrapper, plotting and testing functions are now also capable of handling 2-D input.

## New features
* Support for simulations in 2-D
* Wrapper, plotting and testing functions also adjusted to 2-D

## R CMD Check
0 errors ✔ | 0 warnings ✔ | 3 notes ✖

R CMD check succeeded

# eRTG3D 0.5.8
Replacement of `.fd.bw` with the R base function `grDevices::nclass.FD`, closes `#29`.
Update required packages, documentation and vignettes.

## New features
* `.fd.bw` with the R base function `grDevices::nclass.FD`

## R CMD Check
0 errors ✔ | 0 warnings ✔ | 3 notes ✖
R CMD check succeeded

# eRTG3D 0.5.6
Updated `roxygen2` function documentations with examples that run properly. Declared all import statements properly. `R CMD check --as-cran eRTG3D_0.5.6.tar.gz` runs properly, no errors and warnings only 3 notes.

## New features

* Running examples in roxygen function documentations

## R CMD Check
0 errors ✔ | 0 warnings ✔ | 3 notes ✖

R CMD check succeeded

# eRTG3D 0.5.5
Parallel calculations are now also enabled for Windows, just set `multicore=TRUE` to use it. The 3-D plots are now adjusted to `plotly 4.8.0`, since the coloring is handled differently for markers and lines.

The package now is based on `R 3.5.1` and all packages have been updated.

## New features
* Parallel processing enabled also on Windows OS
* Better coloring in 3-D plots

# eRTG3D 0.5.4
Dependencies updated, documentation updated and variable names in plots adjusted.

## New features
* Consistent variable names in plots

# eRTG3D 0.5.3
Public release of the R package as it was at the submission of the master's thesis:
Unterfinger, M (2018). “3-D Trajectory Simulation in Movement Ecology: Conditional Empirical Random Walk”. Master’s thesis. University of Zurich.

## New features
* Package is available on GitHub
