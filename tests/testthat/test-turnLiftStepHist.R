test_that("turnLiftStepHist() has the correct structure", {
turnLiftStep <- turnLiftStepHist(turn = c(0, 0, 0), lif = c(0, 0, 0), step = c(0, 0, 0),
                                 printDims = FALSE, rm.zeros = TRUE, maxBin = Inf)
expect_equal(turnLiftStep$values$turn, 0)
expect_equal(turnLiftStep$values$lift, 0)
expect_equal(turnLiftStep$values$step, 0)
expect_equal(turnLiftStep$values$prob, 1)
expect_equal(turnLiftStep$tRes, 0)
expect_equal(turnLiftStep$lRes, 0)
expect_equal(turnLiftStep$dRes, 0)
})
