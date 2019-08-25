test_that("parpbapply functions work", {
  # Set up
  n <- 1000
  nNodes <- 2
  df <- data.frame(
    x = seq(1, n, 1),
    y = -seq(1, n, 1)
  )
  #square <- function(x){x*x}
  # Run parallel
  invisible(capture.output(
    s <- parpbsapply(X = df$x, FUN = sum, nNodes = nNodes)
  ))
  invisible(capture.output(
    l <- parpblapply(X = df$x, FUN = sum, nNodes = nNodes)
  ))
  invisible(capture.output(
    a <- parpbapply(X = df, FUN = sum, MARGIN = 1, nNodes = nNodes)
  ))
  # Test
  expect_is(s, "numeric")
  expect_is(l, "list")
  expect_is(a, "numeric")
})
