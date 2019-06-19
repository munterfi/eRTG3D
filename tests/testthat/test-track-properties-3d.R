test_that("track.properties.3d works", {
  properties <- track.properties.3d(niclas)
  expect_is(properties, "data.frame")
  expect_equal(any(apply(properties, class, MARGIN = 2) != "numeric"), FALSE)
  expect_equal(colnames(properties), c("x", "y", "z", "a", "g", "t", "l", "d"))
  expect_equal(nrow(properties), nrow(niclas))
})
