#' Find short tandem repeats (STR)
#'
#' Searches a DNA sequence for short tandem repeats — repeating elements of
#' 1–9 bp occurring three or more times in tandem (microsatellites).
#'
#' @param seq A single DNA sequence as a character string, or a path to a
#'   FASTA file.
#' @param minSTR Minimum length of the repeating unit (bp). Default: `1`.
#' @param maxSTR Maximum length of the repeating unit (bp). Default: `9`.
#' @param minSTRbp Minimum total length of the STR locus (bp). Default: `10`.
#' @param format Output format: `"data.frame"` (default) or `"GRanges"`.
#' @return A `data.frame` (or `GRanges`) with columns: `seq_name`, `start`,
#'   `end`, `strand`, `length`, `spacer`, `num_repeats`, `remainder`,
#'   `subset`.
#' @seealso [find_nonb()]
#' @references Cer et al. (2013) \doi{10.1093/nar/gks955}
#' @export
find_str <- function(seq,
                     minSTR   = 1L,
                     maxSTR   = 9L,
                     minSTRbp = 10L,
                     format   = c("data.frame", "GRanges")) {
  format <- match.arg(format)
  seq <- .resolve_seq(seq)
  results <- lapply(seq, function(s) {
    raw <- .Call("gfa_find_str", s,
                 as.integer(minSTR), as.integer(maxSTR),
                 as.integer(minSTRbp), 3L,
                 PACKAGE = "nonbgfa")
    .rep_to_df(raw, names(s) %||% "seq1")
  })
  df <- do.call(rbind, results)
  if (format == "GRanges") return(to_granges(df))
  df
}
