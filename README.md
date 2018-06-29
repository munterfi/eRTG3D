# Empirically Informed Random Trajectory Generation in 3-D

The empirically informed random trajectory generator in three dimensions (eRTG3D)
is an algorithm to generate realistic random trajectories in a 3-D space
between two given fix points. The trajectory generation is based on
empirical distribution functions extracted from observed trajectories (training data)
and thus reflects the geometrical movement characteristics of the mover.

The eRTG3D algorithm was developed and implemented as an R package within the scope of a master's thesis (Unterfinger, 2018) at the Department of Geography, University of Zurich. The development startet from a 2-D version of the eRTG algorithm by Technitis et al. ([2016](https://doi.org/10.5167/uzh-130652)).

**Functionality** - The **eRTG3D** package contains functions to:

* calculate **movement parameters of 3-D GPS tracking data**, turning angle, lift angle and step length
* **extract distributions** from movement parameters;
    1. **P probability** - The mover's behavior from its perspective
    2. **Q probability** - The pull towards the target
* simulate **Unconditional Empirical Random Walks (UERW)**
* simulate **Conditional Empirical Random Walks (CERW)**
* simulate conditional **gliding and soaring behavior** of birds between two given points
* **statistically test** the simulated tracks against the original input
* **visualize** tracks, simulations and distributions in 3-D and 2-D
* extract **3-D Utilization Distributions (UDs)** from observed or simulated tracking data by means of voxel counting
* project 3-D tracking data into different **Coordinate Reference Systems (CRSs)**
* export data to **sf package objects**; 'sf, data.frames'
* manipulate **extent of raster layers**

## Installing
**Prerequisites** - Software needed:

* [R](https://www.r-project.org/) - R is a free software environment for statistical computing and graphics.
* [RStudio](https://www.rstudio.com/) - Open source and enterprise-ready professional software for R.

**Install Package** - Get development version from github:

```r
library(devtools)
install_github("munterfinger/eRTG3D")
```

## Authors

* **Merlin Unterfinger** - *eRTG3D and R Package* - [munterfinger](http://www.munterfinger.ch)
* **George Technitis** - *2-D eRTG* - [nnneogeorge](https://github.com/nnneogeorge)
* **Dr. Kamran Safi** - *2-D eRTG* - [MPIO](https://www.orn.mpg.de/person/26381/2168)
* **Prof. Dr. Robert Weibel** - *2-D eRTG* - [GIUZ](https://www.geo.uzh.ch/en/studying/spez_master/msc_spez_giscience/People/weibel.html)

## License

This R package is licensed under the GPL (>= 3) License - see the [LICENSE](LICENSE) file for details.

## References

Technitis, G., Weibel, R., Kranstauber, B., and Safi, K. ([2016](https://doi.org/10.5167/uzh-130652)). “An algorithm for empirically informed random trajectory generation between two endpoints”. In: GIScience 2016: Ninth International Conference on Geographic Information Science. Vol. 9. s.n., online. DOI: [10.5167/uzh-130652](https://doi.org/10.5167/uzh-130652).

Unterfinger, M (2018). “3-D Trajectory Simulation in Movement Ecology: Conditional Empirical Random Walk”. Master's thesis. University of Zurich.
