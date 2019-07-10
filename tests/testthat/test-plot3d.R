test_that("plot3d works", {
  expect_is(plot3d(origTrack = niclas), c("plotly", "htmlwidget"))
  expect_is(plot3d(origTrack = niclas, simTrack = niclas, DEM = dem, titleText = "Test"), c("plotly", "htmlwidget"))
})
