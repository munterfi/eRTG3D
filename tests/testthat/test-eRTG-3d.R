test_that("test.eRTG.3d() works", {
  grDevices::pdf(NULL)
  invisible(capture.output(
    t <- test.eRTG.3d(plot2d = TRUE, plot3d = TRUE, plotDensities = TRUE, returnResult = TRUE)
  ))
  expect_is(t, "list")
  expect_equal(nrow(t$crw), nrow(t$cerw))
  expect_equal(any(apply(t$crw, class, MARGIN = 2) != "numeric"), FALSE)
  expect_equal(any(apply(t$cerw, class, MARGIN = 2) != "numeric"), FALSE)
  expect_equal(colnames(t$crw), c("x", "y", "z", "a", "g", "t", "l", "d"))
  expect_equal(colnames(t$cerw), c("x", "y", "z", "a", "g", "t", "l", "d", "p"))
})
