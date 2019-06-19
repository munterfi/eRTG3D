test_that("sim.cond.3d() works", {
  set.seed(123)
  n <- 5
  crw <- sim.crw.3d(nStep = 10000, start = c(0, 0, 0), rTurn = 0.5, rLift = 0.5, meanStep = 1)
  P <- get.track.densities.3d(track = crw)
  invisible(capture.output(Q <- qProb.3d(n.locs = n, crw)))
  crw <- track.properties.3d(crw[1:n, ])
  invisible(capture.output(sim <- sim.cond.3d(n.locs = n,
              start = Reduce(c, crw[1, 1:3]),
              end = Reduce(c, crw[n, 1:3]),
              a0 = crw$a[1],
              g0 = crw$g[1],
              densities = P,
              qProbs = Q)))
  expect_is(sim, "data.frame")
  expect_equal(any(apply(sim, class, MARGIN = 2) != "numeric"), FALSE)
  expect_equal(colnames(sim), c("x", "y", "z", "a", "g", "t", "l", "d", "p"))
  expect_equal(nrow(sim), n)
})
