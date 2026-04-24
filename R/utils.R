# Utility functions shared across motif finders.

# Default parameter values matching the original gfa CLI defaults.
.gfa_defaults <- list(
  # G-Quadruplex
  minGQrep      = 3L,
  maxGQspacer   = 7L,
  # Mirror Repeat
  minMRrep      = 10L,
  maxMRspacer   = 100L,
  # Inverted Repeat
  minIRrep      = 6L,
  maxIRspacer   = 100L,
  shortIRcut    = 9L,
  shortIRspacer = 4L,
  # Direct Repeat
  minDRrep      = 10L,
  maxDRrep      = 300L,
  maxDRspacer   = 10L,
  # A-Phased Repeat
  minATracts    = 3L,
  minATractSep  = 10L,
  maxATractSep  = 11L,
  maxAPRlen     = 9L,
  minAPRlen     = 3L,
  # Z-DNA
  minZlen       = 10L,
  # Short Tandem Repeat
  minSTR        = 1L,
  maxSTR        = 9L,
  minSTRbp      = 10L,
  # Subset classification thresholds
  minCruciformRep       = 6L,
  maxCruciformSpacer    = 4L,
  minTriplexYRpercent   = 10L,
  maxTriplexSpacer      = 8L,
  maxSlippedSpacer      = 0L
)

# NULL-coalescing operator used in finder functions
`%||%` <- function(x, y) if (!is.null(x) && length(x) > 0L) x else y

# Internal: accept either a file path or a raw DNA character vector.
# Returns a named character vector (names = seq identifiers).
.resolve_seq <- function(seq) {
  if (length(seq) == 1L &&
      grepl("\\.(fa|fasta|fna)(\\.gz)?$", seq, ignore.case = TRUE) &&
      file.exists(seq)) {
    read_fasta(seq)
  } else {
    if (is.null(names(seq))) names(seq) <- paste0("seq", seq_along(seq))
    seq
  }
}


#' Read a FASTA file
#'
#' Parses a FASTA file (single or multi-sequence) into a named character vector.
#'
#' @param path Path to a FASTA file.
#' @return A named character vector.  Names are the first word of each `>`
#'   header line; values are the concatenated sequence strings.
#' @examples
#' fa <- system.file("extdata", "gfa_test.fasta", package = "nonbgfa")
#' seq <- read_fasta(fa)
#' nchar(seq)
#' @export
read_fasta <- function(path) {
  lines <- readLines(path)
  header_idx <- grep("^>", lines)
  if (length(header_idx) == 0L)
    stop("Not a valid FASTA file: no '>' header lines found in ", path)

  starts <- header_idx
  ends   <- c(header_idx[-1L] - 1L, length(lines))
  names  <- sub("^>([^ \t]+).*", "\\1", lines[header_idx])

  seqs <- mapply(function(s, e) {
    paste(lines[(s + 1L):e], collapse = "")
  }, starts, ends, SIMPLIFY = TRUE)

  stats::setNames(seqs, names)
}

# Internal: convert the named list returned from C into a data.frame,
# prepending the seq_name column.
.rep_to_df <- function(rep_list, seq_name) {
  df <- as.data.frame(rep_list, stringsAsFactors = FALSE)
  df$subset <- as.logical(df$subset)
  if (nrow(df) == 0L) {
    cbind(seq_name = character(0), df, stringsAsFactors = FALSE)
  } else {
    cbind(seq_name = seq_name, df, stringsAsFactors = FALSE)
  }
}

#' Convert a nonbgfa data.frame to a GRanges object
#'
#' Requires the \pkg{GenomicRanges} package (Bioconductor).
#'
#' @param df A `data.frame` returned by one of the `find_*()` functions.
#' @return A [GenomicRanges::GRanges] object.
#' @export
to_granges <- function(df) {
  rlang::check_installed("GenomicRanges",
                         reason = "to convert results to GRanges format")
  rlang::check_installed("IRanges",
                         reason = "to convert results to GRanges format")
  GenomicRanges::GRanges(
    seqnames = df$seq_name,
    ranges   = IRanges::IRanges(start = df$start, end = df$end),
    strand   = df$strand,
    df[, setdiff(names(df), c("seq_name", "start", "end", "strand"))]
  )
}
