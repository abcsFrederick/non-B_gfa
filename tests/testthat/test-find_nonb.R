test_that("find_nonb() returns a named list with all seven motif types", {
  result <- find_nonb(example_fasta)
  expect_type(result, "list")
  expect_named(result, c("IR", "MR", "DR", "GQ", "Z", "STR", "APR"))
})

test_that("find_nonb() respects the skip argument", {
  result <- find_nonb(example_fasta, skip = c("MR", "DR", "STR"))
  expect_named(result, c("IR", "GQ", "Z", "APR"))
})
