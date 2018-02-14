# eRTG3D

Empirically Informed Random Trajectory Generation in 3-D

## About

The empirically informed random trajectory generator in three dimensions (eRTG3D)
is an algorithm to generate realistic random trajectories in a 3-D space
between two given fix points in space. The trajectory generation is based on
empirical distribution functions extracted from observed trajectories (training data)
and thus reflects the geometrical movement characteristics of the mover.

### Prerequisites

* [R](https://www.r-project.org/) - R is a free software environment for statistical computing and graphics.
* [RStudio](https://www.rstudio.com/) - Open source and enterprise-ready professional software for R.

### R Packages needed

```
install.packages("doParallel")
install.packages("ggplot2")
install.packages("parallel")
install.packages("plotly")
install.packages("raster")
install.packages("sf")
```

### Test if eRTG3D is running properly

Singlecore:
```
test.eRTG3D(multicore = FALSE, returnResult = FALSE)
```

And multicore
```
test.eRTG3D(multicore = TRUE, returnResult = FALSE)
```

## Authors

* **Merlin Unterfinger** - *3-D version* - [munterfinger](https://github.com/munterfinger)
* **George Technitis** - *original eRTG* - [nnneogeorge](https://github.com/nnneogeorge)
* **Kamran Safi** - *original eRTG*
* **Robert Weibel** - *original eRTG*

## License

This R package is licensed under the GPL (>= 3) License - see the [LICENSE](LICENSE) file for details.
