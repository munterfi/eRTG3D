test_that("sf linkage functions work", {
  # To sf
  niclas.sf <- track2sf.3d(niclas, CRS = "+init=epsg:2056")
  expect_equal(is.sf.3d(niclas.sf), TRUE)
  expect_equal(nrow(niclas), nrow(niclas.sf))
  expect_is(niclas.sf$geometry, c("sfc_POINT", "sfc"))
  niclas.sf.m <- track2sf.3d(as.matrix(niclas), CRS = "+init=epsg:2056")
  expect_equal(is.sf.3d(niclas.sf.m), TRUE)
  expect_equal(nrow(niclas), nrow(niclas.sf.m))
  expect_is(niclas.sf.m$geometry, c("sfc_POINT", "sfc"))
  # from sf
  niclas.df <- sf2df.3d(niclas.sf)
  expect_equal(is.sf.3d(niclas.df), FALSE)
  niclas.df.m <- sf2df.3d(niclas.sf.m)
  expect_equal(is.sf.3d(niclas.df.m), FALSE)
  # transform CRS
  row.names(niclas) <- NULL
  expect_equal(niclas,
               transformCRS.3d(
                 transformCRS.3d(niclas, fromCRS="+init=epsg:2056", toCRS="+init=epsg:4326"),
                 toCRS="+init=epsg:2056", fromCRS="+init=epsg:4326")[, 1:3])
})
