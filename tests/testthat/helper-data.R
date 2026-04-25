library(nonbgfa)

# Paths to reference data bundled with the package
extdata_dir <- system.file("extdata", package = "nonbgfa")
example_fasta <- file.path(extdata_dir, "gfa_test.fasta")

# Load expected TSV outputs produced by the original gfa binary
expected <- list(
  IR  = read.table(file.path(extdata_dir, "gfa_test_IR.tsv"),
                   header = TRUE, sep = "\t", comment.char = "#"),
  MR  = read.table(file.path(extdata_dir, "gfa_test_MR.tsv"),
                   header = TRUE, sep = "\t", comment.char = "#"),
  DR  = read.table(file.path(extdata_dir, "gfa_test_DR.tsv"),
                   header = TRUE, sep = "\t", comment.char = "#"),
  GQ  = read.table(file.path(extdata_dir, "gfa_test_GQ.tsv"),
                   header = TRUE, sep = "\t", comment.char = "#"),
  Z   = read.table(file.path(extdata_dir, "gfa_test_Z.tsv"),
                   header = TRUE, sep = "\t", comment.char = "#"),
  STR = read.table(file.path(extdata_dir, "gfa_test_STR.tsv"),
                   header = TRUE, sep = "\t", comment.char = "#"),
  APR = read.table(file.path(extdata_dir, "gfa_test_APR.tsv"),
                   header = TRUE, sep = "\t", comment.char = "#")
)
