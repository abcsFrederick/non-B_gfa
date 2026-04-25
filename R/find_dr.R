#' Find direct repeats (slipped DNA)
#'
#' Searches a DNA sequence for direct repeats — identical sequence runs in
#' tandem — which can form slipped-strand (hairpin) structures.
#'
#' @param seq A single DNA sequence as a character string, or a path to a
#'   FASTA file.
#' @param minDRrep Minimum length of each repeat unit. Default: `10`.
#' @param maxDRrep Maximum length of each repeat unit. Default: `300`.
#' @param maxDRspacer Maximum spacer length between repeats. Default: `10`.
#' @param maxSlippedSpacer Maximum spacer for slipped subset flag. Default: `0`.
#' @param format Output format: `"data.frame"` (default) or `"GRanges"`.
#' @return A `data.frame` (or `GRanges`) with columns: `seq_name`, `start`,
#'   `end`, `strand`, `length`, `spacer`, `num_repeats`, `remainder`,
#'   `subset`.
#' @seealso [find_nonb()]
#' @references Cer et al. (2013) \doi{10.1093/nar/gks955}
#' @export
find_dr <- function(seq,
                    minDRrep         = 10L,
                    maxDRrep         = 300L,
                    maxDRspacer      = 10L,
                    maxSlippedSpacer = 0L,
                    format           = c("data.frame", "GRanges")) {
  format <- match.arg(format)
  seq <- .resolve_seq(seq)
  results <- lapply(seq, function(s) {
    raw <- .Call("gfa_find_dr", s,
                 as.integer(minDRrep), as.integer(maxDRrep),
                 as.integer(maxDRspacer), as.integer(maxSlippedSpacer),
                 PACKAGE = "nonbgfa")
    .rep_to_df(raw, names(s) %||% "seq1")
  })
  df <- do.call(rbind, results)
  if (format == "GRanges") return(to_granges(df))
  df
}
