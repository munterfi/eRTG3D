test_that("plot2d works", {
  grDevices::pdf(NULL)
  expect_is(plot2d(origTrack = niclas), c("gg", "ggplot"))
  expect_is(plot2d(origTrack = niclas, simTrack = niclas, DEM = dem, titleText = "Test"), c("gg", "ggplot"))
  expect_is(plot2d(origTrack = niclas, simTrack = niclas, BG = dem, titleText = "Test"), c("gg", "ggplot"))
})
