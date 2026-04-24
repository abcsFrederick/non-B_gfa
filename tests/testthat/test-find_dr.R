test_that("find_dr() errors until C interface is implemented", {
  result <- find_dr(example_fasta)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), nrow(expected$DR))
  expect_equal(result$start, expected$DR$Start)
  expect_equal(result$end,   expected$DR$Stop)
})
