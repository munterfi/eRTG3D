test_that("point-cloud-analysis works", {
  crws <- lapply(X = seq(1:100), FUN = function(X) {
    sim.crw.3d(nStep = 100, rTurn = 0.99, rLift = 0.99, meanStep = 0.1)
  })
  points <- do.call("rbind", crws)
  extent <- raster::extent(c(-10, 10, -10, 10))
  # Voxel count standartize
  invisible(capture.output(
    ud <- voxelCount(points, extent, xyRes = 5, zMin = -10, zMax = 10, standartize = TRUE)
  ))
  expect_is(ud, "RasterStack")
  # Voxel count
  invisible(capture.output(
    ud <- voxelCount(points, extent, xyRes = 5, zMin = -10, zMax = 10)
  ))
  expect_is(ud, "RasterStack")
  # Chi maps
  invisible(capture.output(
    chi <- chiMaps(ud)
  ))
  expect_is(chi, "RasterStack")
  # Log raster stack
  expect_is(logRasterStack(abs(chi)), "RasterStack")
  # Raster stack plots
  expect_is(plotRaster(chi, centerColorBar = TRUE), c("gtable", "gTree", "grob", "gDesc"))
  expect_is(plotRaster(chi, centerColorBar = FALSE), c("gtable", "gTree", "grob", "gDesc"))
  expect_is(plotRaster(chi, ncol = 1), c("gtable", "gTree", "grob", "gDesc"))
})
