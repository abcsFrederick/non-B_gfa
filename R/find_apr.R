#' Find A-phased repeats (bent DNA)
#'
#' Searches a DNA sequence for A-phased repeats — A/T tracts recurrently
#' spaced at ~10–11 bp intervals (one helical turn) — which cause
#' macroscopic DNA bending.
#'
#' @param seq A single DNA sequence as a character string, or a path to a
#'   FASTA file.
#' @param minATracts Minimum number of consecutive A-tracts. Default: `3`.
#' @param minATractSep Minimum center-to-center separation between A-tracts
#'   (bp). Default: `10`.
#' @param maxATractSep Maximum center-to-center separation between A-tracts
#'   (bp). Default: `11`.
#' @param minAPRlen Minimum number of consecutive A's in an A-tract.
#'   Default: `3`.
#' @param maxAPRlen Maximum number of consecutive A's in an A-tract.
#'   Default: `9`.
#' @param format Output format: `"data.frame"` (default) or `"GRanges"`.
#' @return A `data.frame` (or `GRanges`) with columns: `seq_name`, `start`,
#'   `end`, `strand`, `length`, `spacer`, `num_repeats`, `remainder`,
#'   `subset`.
#' @seealso [find_nonb()]
#' @references Cer et al. (2013) \doi{10.1093/nar/gks955}
#' @export
find_apr <- function(seq,
                     minATracts   = 3L,
                     minATractSep = 10L,
                     maxATractSep = 11L,
                     minAPRlen    = 3L,
                     maxAPRlen    = 9L,
                     format       = c("data.frame", "GRanges")) {
  format <- match.arg(format)
  seq <- .resolve_seq(seq)
  results <- lapply(seq, function(s) {
    raw <- .Call("gfa_find_apr", s,
                 as.integer(minAPRlen), as.integer(maxAPRlen),
                 as.integer(minATracts),
                 PACKAGE = "nonbgfa")
    .rep_to_df(raw, names(s) %||% "seq1")
  })
  df <- do.call(rbind, results)
  if (format == "GRanges") return(to_granges(df))
  df
}
