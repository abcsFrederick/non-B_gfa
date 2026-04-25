test_that("find_mr() errors until C interface is implemented", {
  result <- find_mr(example_fasta)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), nrow(expected$MR))
  expect_equal(result$start, expected$MR$Start)
  expect_equal(result$end,   expected$MR$Stop)
})
