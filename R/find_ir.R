#' Find inverted repeats (cruciform DNA)
#'
#' Searches a DNA sequence for inverted repeats — complementary palindromic
#' sequences on the same strand separated by a spacer — which can form
#' cruciform structures.
#'
#' @param seq A single DNA sequence as a character string, or a path to a
#'   FASTA file (detected automatically if the string ends in `.fa`, `.fasta`,
#'   `.fna`, or `.fa.gz`).
#' @param minIRrep Minimum length of each repeat arm. Default: `6`.
#' @param maxIRspacer Maximum spacer length between arms. Default: `100`.
#' @param shortIRcut Maximum arm length to be classified as a "short" IR.
#'   Default: `9`.
#' @param shortIRspacer Maximum spacer for short IRs. Default: `4`.
#' @param minCruciformRep Minimum arm length for cruciform subset flag.
#'   Default: `6`.
#' @param maxCruciformSpacer Maximum spacer for cruciform subset flag.
#'   Default: `4`.
#' @param format Output format: `"data.frame"` (default) or `"GRanges"`.
#' @return A `data.frame` (or `GRanges`) with columns: `seq_name`, `start`,
#'   `end`, `strand`, `length`, `spacer`, `num_repeats`, `remainder`,
#'   `subset`.
#' @seealso [find_nonb()] for running all motif finders at once.
#' @references Cer et al. (2013) Non-B DB v2.0. *Nucleic Acids Research*,
#'   41(D1):D94–D100. \doi{10.1093/nar/gks955}
#' @export
find_ir <- function(seq,
                    minIRrep          = 6L,
                    maxIRspacer       = 100L,
                    shortIRcut        = 9L,
                    shortIRspacer     = 4L,
                    minCruciformRep   = 6L,
                    maxCruciformSpacer = 4L,
                    format            = c("data.frame", "GRanges")) {
  format <- match.arg(format)
  seq <- .resolve_seq(seq)
  results <- lapply(seq, function(s) {
    raw <- .Call("gfa_find_ir", s,
                 as.integer(minIRrep), as.integer(maxIRspacer),
                 as.integer(shortIRcut), as.integer(shortIRspacer),
                 as.integer(minCruciformRep), as.integer(maxCruciformSpacer),
                 PACKAGE = "nonbgfa")
    .rep_to_df(raw, names(s) %||% "seq1")
  })
  df <- do.call(rbind, results)
  if (format == "GRanges") return(to_granges(df))
  df
}
