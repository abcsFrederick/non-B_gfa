test_that("find_ir() errors until C interface is implemented", {
  result <- find_ir(example_fasta)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("seq_name", "start", "end", "strand",
                          "length", "spacer", "num_repeats", "remainder", "subset"))
  expect_equal(nrow(result), nrow(expected$IR))
  expect_equal(result$start, expected$IR$Start)
  expect_equal(result$end,   expected$IR$Stop)
})

test_that("find_ir() returns GRanges when format = 'GRanges'", {
  skip_if_not_installed("GenomicRanges")
  result <- find_ir(example_fasta, format = "GRanges")
  expect_s4_class(result, "GRanges")
})
