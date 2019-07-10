test_that("track points calculations in respect to target work", {
  # Distance to point
  d_ground <- dist2point.3d(niclas, c(median(niclas$x), median(niclas$y), median(niclas$z)), groundDistance = TRUE)
  d_air <- dist2point.3d(niclas, c(median(niclas$x), median(niclas$y), median(niclas$z)), groundDistance = FALSE)
  expect_equal(length(d_ground), length(d_air))
  expect_equal(length(d_ground), nrow(niclas))
  expect_is(d_ground, "numeric")
  expect_is(d_air, "numeric")
  # Distance to target
  d <- dist2target.3d(niclas)
  expect_equal(length(d), nrow(niclas))
  expect_is(d, "numeric")
  # Lift angle to target
  l <- lift2target.3d(niclas)
  expect_equal(length(l), nrow(niclas))
  expect_is(l, "numeric")
  # Turn angle to target
  t <- turn2target.3d(niclas)
  expect_equal(length(t), nrow(niclas))
  expect_is(t, "numeric")
})
