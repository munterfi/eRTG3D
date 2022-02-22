test_that("get.densities.3d() returns the correct structure", {
  test_vec <- rep(0, 5)
  test_delta_vec <- rep(0, 4)
  P <- get.densities.3d(
    turnAngle = test_vec, liftAngle = test_vec, stepLength = test_vec,
    deltaLift = test_delta_vec, deltaTurn = test_delta_vec, deltaStep = test_delta_vec,
    gradientAngle = NULL, heightEllipsoid = NULL, heightTopo = NULL, maxBin = Inf
  )
  # Test P probability
  expect_equal(names(P$tldCube), c("values", "tRes", "lRes", "dRes"))
  expect_is(P$autoT, "function")
  expect_is(P$autoL, "function")
  expect_is(P$autoD, "function")
  expect_is(P$gDens, "function")
  expect_is(P$hDistEllipsoid, "function")
  expect_is(P$hDistTopo, "function")
  # Gradient needs to be between x > 0 & x < pi
  expect_equal(P$gDens(pi), 0)
  expect_equal(P$gDens(0), 0)
  expect_equal(P$gDens(pi / 2), 1)
  # If no input DEM output should ALWAYS be 1.
  expect_equal(P$hDistEllipsoid(rnorm(1)), 1)
  expect_equal(P$hDistTopo(rnorm(1)), 1)
})
