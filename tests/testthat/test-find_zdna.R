test_that("find_zdna() errors until C interface is implemented", {
  result <- find_zdna(example_fasta)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), nrow(expected$Z))
  expect_equal(result$start, expected$Z$Start)
  expect_equal(result$end,   expected$Z$Stop)
})
