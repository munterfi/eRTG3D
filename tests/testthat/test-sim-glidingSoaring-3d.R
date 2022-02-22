test_that("sim.glidingSoaring.3d works", {
  set.seed(123)
  # Soaring
  t <- seq(0, 10 * pi, pi / 6)
  skew <- 10
  mu <- 0
  sd <- 0.1
  soaring <- data.frame(
    x = sin(t) * skew + cumsum(rnorm(length(t), mean = mu, sd = 0.3)),
    y = cos(t) * skew + cumsum(rnorm(length(t), mean = mu, sd = 0.5)),
    z = t + cumsum(rnorm(length(t), mean = mu, sd = 0.1))
  )
  soaring <- track.properties.3d(soaring)
  # Gliding
  n <- nrow(soaring)
  gliding <- sim.crw.3d(nStep = n, rTurn = 0.99, rLift = 0.99, meanStep = mean(soaring$d, na.rm = TRUE))
  gliding$z <- abs(gliding$z) * -0.1
  gliding <- track.properties.3d(gliding)
  # Soaring2
  soaring2 <- data.frame(
    x = sin(-t) * skew + cumsum(rnorm(length(t), mean = mu, sd = 0.3)),
    y = cos(-t) * skew + cumsum(rnorm(length(t), mean = mu, sd = 0.5)),
    z = t + cumsum(rnorm(length(t), mean = mu, sd = 0.1))
  )
  # P Probs
  P.gliding <- get.track.densities.3d(gliding)
  P.soaring <- get.section.densities.3d(list(soaring, soaring2, gliding))
  invisible(capture.output(
    uerw.soaring <- sim.uncond.3d(n * 1500, densities = P.soaring, a0 = soaring$a[1], g0 = soaring$g[1])
  ))
  P.soaring <- get.track.densities.3d(uerw.soaring)
  # Q probs
  invisible(capture.output(
    uerw <- sim.uncond.3d(n * 1500, densities = P.gliding, a0 = gliding$a[1], g0 = gliding$g[1])
  ))
  P.gliding <- get.track.densities.3d(uerw)
  invisible(capture.output(
    Q.gliding <- qProb.3d(n.locs = round(n * 1), uerw)
  ))
  # Modelayer
  r1 <- r2 <- raster::raster(xmn = -250, xmx = 250, ymn = -250, ymx = 250, resolution = 10, vals = NULL)
  raster::values(r1) <- 1:raster::ncell(r1)
  raster::values(r2) <- raster::ncell(r1):1
  MODE <- ((r1 * r2) > 1559000) + 1
  # Simulate soaring and gliding
  glideRatio <- get.glideRatio.3d(gliding)
  invisible(capture.output(
    soaringGliding <- suppressWarnings(
      n.sim.glidingSoaring.3d(
        n.sim = 2, MODE = MODE, dGliding = P.gliding, dSoaring = P.soaring, qGliding = Q.gliding,
        start = c(25, 25, 0), end = Reduce(c, c(tail(gliding, 1)[, 1:2] + 25, 0)),
        a0 = soaring$a[1], g0 = soaring$g[1], error = TRUE, smoothTransition = TRUE,
        glideRatio = glideRatio, DEM = NULL, BG = NULL, verbose = TRUE
      )
    )
  ))
  expect_is(soaringGliding, "list")
})
