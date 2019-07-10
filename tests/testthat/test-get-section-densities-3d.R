test_that("get.section.densities.3d works", {
  # Split track and delete short sections (nrow = 1)
  timeLag <- movingMedian(rnorm(nrow(niclas)-1), window = 3)
  trackList <- track.split.3d(niclas, timeLag = timeLag)
  trackList <- trackList[sapply(trackList, nrow) > 1]
  # Section densities
  P <- get.section.densities.3d(trackSections = trackList, gradientDensity = TRUE, heightDistEllipsoid = TRUE, DEM = dem)
  # Test P probability
  expect_equal(names(P$tldCube), c("values", "tRes", "lRes", "dRes"))
  expect_is(P$autoT, "function")
  expect_is(P$autoL, "function")
  expect_is(P$autoD, "function")
  expect_is(P$gDens, "function")
  expect_is(P$hDistEllipsoid, "function")
  expect_is(P$hDistTopo, "function")
  # Section densities with DEM and height distributions
  P <- get.section.densities.3d(trackSections = trackList, gradientDensity = FALSE, heightDistEllipsoid = FALSE, DEM = NULL)
  # Gradient needs to be between x > 0 & x < pi
  expect_equal(P$gDens(pi), 0)
  expect_equal(P$gDens(0), 0)
  expect_equal(P$gDens(pi/2), 1)
  # If no input DEM output should ALWAYS be 1.
  expect_equal(P$hDistEllipsoid(rnorm(1)), 1)
  expect_equal(P$hDistTopo(rnorm(1)), 1)
})
