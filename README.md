# Empirically Informed Random Trajectory Generation in 3-D

The empirically informed random trajectory generator in three dimensions (eRTG3D)
is a R package to generate realistic random trajectories in a 3-D space
between two given fix points. The trajectory generation is based on
empirical distribution functions extracted from observed trajectories (training data)
and thus reflects the geometrical movement characteristics of the mover.

The **eRTG3D** package contains functions to:

* calculate **properties of 3-D GPS tracking data**, turn angle, lift angle and step length
* **extract distributions** from properties; (1) P probability - the movers behavior from its perspective, and (2) Q probability - the pull towards the target
* simulate **Unconditioned Empirical Random Walks**
* simulate **Conditioned Empirical Random Walks**
* simulate conditioned **gliding and soaring behavior** of birds between to given points
* **verify the results** statistically
* **visualize** tracks, simulations and distributions in 3-D and 2-D
* extract **3-D Space Utilization Distributions** from observed or simulated tracking data
* project 3-D tracking data into different **Coordinate Reference Systems**
* export data to **sf package objects**; 'sf, data.frames'
* manipulate **extent of raster layers**

## Installing
### Prerequisites

* [R](https://www.r-project.org/) - R is a free software environment for statistical computing and graphics.
* [RStudio](https://www.rstudio.com/) - Open source and enterprise-ready professional software for R.

### Package

Install either from CRAN with:
```r
install.packages("eRTG3D")
```
this will install binary packages on Windows and MacOS, unless you configured R such that it tries to install source packages; in that case, see below.

Install development versions from github with
```r
library(devtools)
install_github("munterfinger/eRTG3D")
```

## Authors

* **Merlin Unterfinger** - *3-D version* - [munterfinger](https://github.com/munterfinger)
* **George Technitis** - *2-D eRTG* - [nnneogeorge](https://github.com/nnneogeorge)
* **Kamran Safi** - *2-D eRTG*
* **Robert Weibel** - *2-D eRTG*

## License

This R package is licensed under the GPL (>= 3) License - see the [LICENSE](LICENSE) file for details.
