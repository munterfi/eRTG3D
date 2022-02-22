test_that("reproduce.track.3d works", {
  grDevices::pdf(NULL)
  set.seed(123)
  # Singlecore
  invisible(capture.output(
    repr_singlecore <- reproduce.track.3d(
      n.sim = 5, track = niclas, parallel = FALSE, error = TRUE,
      DEM = dem, BG = (dem < 1000), filterDeadEnds = TRUE, plot2d = TRUE,
      plot3d = TRUE, maxBin = 25, gradientDensity = TRUE
    )
  ))
  expect_is(repr_singlecore, "list")
  expect_equal(any(sapply(repr_singlecore, class) != "data.frame"), FALSE)
  expect_equal(any(sapply(repr_singlecore, nrow) != nrow(niclas)), FALSE)

  # Singlecore
  invisible(capture.output(
    repr_singlecore_minimum_opts <- reproduce.track.3d(
      n.sim = 5, track = niclas, parallel = FALSE, error = FALSE,
      DEM = NULL, BG = NULL, filterDeadEnds = FALSE, plot2d = FALSE,
      plot3d = FALSE, maxBin = 25, gradientDensity = FALSE
    )
  ))
  expect_is(repr_singlecore_minimum_opts, "list")

  # Parallel
  invisible(capture.output(
    repr_parallel <- reproduce.track.3d(n.sim = 5, track = niclas, DEM = dem, parallel = 2, gradientDensity = TRUE, maxBin = 25)
  ))
  expect_is(repr_parallel, "list")
  expect_equal(any(sapply(repr_parallel, class) != "data.frame"), FALSE)
  expect_equal(any(sapply(repr_parallel, nrow) != nrow(niclas)), FALSE)
})
