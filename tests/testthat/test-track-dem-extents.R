test_that("track.extent and dem works", {
  e <- track.extent(niclas, zAxis = FALSE)
  expect_is(e, "Extent")
  e.z <- track.extent(niclas, zAxis = TRUE)
  expect_is(e.z, "matrix")
  dem.cropped <- dem2track.extent(dem, niclas, buffer = 100)
  expect_is(dem.cropped, "RasterLayer")
})
