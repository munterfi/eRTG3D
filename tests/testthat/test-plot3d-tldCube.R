test_that("plot3d.tldCube works", {
  P <- get.track.densities.3d(niclas, maxBin = Inf, gradientDensity = FALSE)
  expect_is(plot3d.tldCube(P$tldCube), c("plotly", "htmlwidget"))
})
