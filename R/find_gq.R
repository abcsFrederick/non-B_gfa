#' Find G-quadruplex motifs
#'
#' Searches a DNA sequence for G-quadruplex (G4) motifs — runs of guanines
#' separated by short spacers — on both strands. G4 structures are formed by
#' four G-runs stacked in a planar tetrad arrangement.
#'
#' @param seq A single DNA sequence as a character string, or a path to a
#'   FASTA file.
#' @param minGQrep Minimum number of consecutive G's to form a G-run.
#'   Default: `3`.
#' @param maxGQspacer Maximum spacer between G-runs. Default: `7`.
#' @param format Output format: `"data.frame"` (default) or `"GRanges"`.
#' @return A `data.frame` (or `GRanges`) with columns: `seq_name`, `start`,
#'   `end`, `strand`, `length`, `spacer`, `num_repeats`, `remainder`,
#'   `subset`.
#' @seealso [find_nonb()]
#' @references Cer et al. (2013) \doi{10.1093/nar/gks955}
#' @export
find_gq <- function(seq,
                    minGQrep    = 3L,
                    maxGQspacer = 7L,
                    format      = c("data.frame", "GRanges")) {
  format <- match.arg(format)
  seq <- .resolve_seq(seq)
  results <- lapply(seq, function(s) {
    raw <- .Call("gfa_find_gq", s,
                 as.integer(minGQrep), as.integer(maxGQspacer),
                 PACKAGE = "nonbgfa")
    .rep_to_df(raw, names(s) %||% "seq1")
  })
  df <- do.call(rbind, results)
  if (format == "GRanges") return(to_granges(df))
  df
}
