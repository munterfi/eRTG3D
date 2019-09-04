test_that("plot3d.densities works", {
  grDevices::pdf(NULL)
  expect_equal(plot3d.densities(track1 = niclas, track2 = niclas, autodifferences = FALSE), NULL)
  expect_equal(plot3d.densities(track1 = niclas, track2 = niclas, autodifferences = TRUE), NULL)
})
