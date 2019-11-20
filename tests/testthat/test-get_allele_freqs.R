test_that("Basic usage works", {
  library(adegenet)
  data(nancycats)

  df <- get_allele_freqs(nancycats)

  expect_equal(dim(df)[1], 1836)
  expect_equal(dim(df)[2], 4)
})
