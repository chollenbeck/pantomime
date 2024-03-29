test_that("ind_stats works", {
  library(adegenet)
  data(nancycats)

  df <- get_ind_stats(nancycats)

  expect_equal(dim(df)[1], 237)
  expect_equal(dim(df)[2], 5)
})
