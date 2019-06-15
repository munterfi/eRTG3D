test_that("qProb.3d() works", {
  n <- 3
  test_vec <- rep(0, n)
  invisible(capture.output(Q <- qProb.3d(niclas, n)))
  expect_is(Q, "list")
  expect_is(Q[[1]], "list")
  expect_equal(length(Q), n-2)
  expect_equal(names(Q[[1]]),
               names(turnLiftStepHist(turn=test_vec,
                                      lift=test_vec,
                                      step=test_vec,
                                      printDims = TRUE,
                                      rm.zeros = TRUE,
                                      maxBin = Inf)))
})
