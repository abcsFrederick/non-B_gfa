test_that("find_str() errors until C interface is implemented", {
  result <- find_str(example_fasta)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), nrow(expected$STR))
  expect_equal(result$start, expected$STR$Start)
  expect_equal(result$end,   expected$STR$Stop)
})
