test_that("reproduce.track.3d works", {
  set.seed(123)
  # Singlecore
  invisible(capture.output(
    repr_singlecore <- reproduce.track.3d(n.sim = 5, track = niclas, DEM = dem, gradientDensity = TRUE, maxBin = 25)
  ))
  expect_is(repr_singlecore, "list")
  expect_equal(any(sapply(repr_singlecore, class) != "data.frame"), FALSE)
  expect_equal(any(sapply(repr_singlecore, nrow) != nrow(niclas)), FALSE)
  
  # Multicore
  # invisible(capture.output(
  #   repr_multicore <- reproduce.track.3d(n.sim = 5, track = niclas, DEM = dem, multicore = TRUE, gradientDensity = TRUE, maxBin = 25)
  # ))
  # expect_is(repr_multicore, "list")
  # expect_equal(any(sapply(repr_multicore, class) != "data.frame"), FALSE)
  # expect_equal(any(sapply(repr_multicore, nrow) != nrow(niclas)), FALSE)
})
