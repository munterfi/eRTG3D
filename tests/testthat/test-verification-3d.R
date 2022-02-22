test_that("test.verification.3d works", {
  ks <- test.verification.3d(niclas, niclas, alpha = 0.05, plot = TRUE, test = "ks")
  ttest <- test.verification.3d(niclas, niclas, alpha = 0.05, plot = TRUE, test = "ttest")
  expect_is(ks, "list")
  expect_is(ttest, "list")
})
