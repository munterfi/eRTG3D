test_that("sim.uncond.3d() works", {
  n <- 10
  invisible(capture.output(
    sim <- sim.uncond.3d(n,
      start = c(0, 0, 0), a0 = pi / 2, g0 = pi / 2,
      densities = get.track.densities.3d(niclas)
    )
  ))
  expect_is(sim, "data.frame")
  expect_equal(any(apply(sim, class, MARGIN = 2) != "numeric"), FALSE)
  expect_equal(colnames(sim), c("x", "y", "z", "a", "g", "t", "l", "d", "p"))
  expect_equal(nrow(sim), n)
})
