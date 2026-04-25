#' Find all non-B DNA forming motifs
#'
#' A unified wrapper that runs any combination of the seven non-B DNA motif
#' finders on a single DNA sequence or FASTA file and returns a named list of
#' results.
#'
#' @param seq A single DNA sequence as a character string, or a path to a
#'   FASTA file.
#' @param format Output format for each element of the returned list:
#'   `"data.frame"` (default) or `"GRanges"`.
#' @param skip A character vector of motif types to skip. Valid values:
#'   `"IR"`, `"MR"`, `"DR"`, `"GQ"`, `"Z"`, `"STR"`, `"APR"`.
#' @param ir_args,mr_args,dr_args,gq_args,zdna_args,str_args,apr_args
#'   Named lists of additional arguments forwarded to the corresponding
#'   `find_*()` function. See each function's documentation for available
#'   parameters.
#' @return A named list with up to seven elements (`IR`, `MR`, `DR`, `GQ`,
#'   `Z`, `STR`, `APR`). Skipped motif types are absent from the list.
#'   Each element is a `data.frame` or `GRanges` as controlled by `format`.
#' @examples
#' \dontrun{
#' fasta <- system.file("extdata", "gfa_test.fasta", package = "nonbgfa")
#' results <- find_nonb(fasta)
#' results$IR   # inverted repeats
#' results$GQ   # G-quadruplexes
#' }
#' @seealso [find_ir()], [find_mr()], [find_dr()], [find_gq()],
#'   [find_zdna()], [find_str()], [find_apr()]
#' @references Cer et al. (2013) Non-B DB v2.0. *Nucleic Acids Research*,
#'   41(D1):D94–D100. \doi{10.1093/nar/gks955}
#' @export
find_nonb <- function(seq,
                      format    = c("data.frame", "GRanges"),
                      skip      = character(0),
                      ir_args   = list(),
                      mr_args   = list(),
                      dr_args   = list(),
                      gq_args   = list(),
                      zdna_args = list(),
                      str_args  = list(),
                      apr_args  = list()) {
  format <- match.arg(format)
  results <- list()

  if (!"IR"  %in% skip) results$IR  <- do.call(find_ir,   c(list(seq = seq, format = format), ir_args))
  if (!"MR"  %in% skip) results$MR  <- do.call(find_mr,   c(list(seq = seq, format = format), mr_args))
  if (!"DR"  %in% skip) results$DR  <- do.call(find_dr,   c(list(seq = seq, format = format), dr_args))
  if (!"GQ"  %in% skip) results$GQ  <- do.call(find_gq,   c(list(seq = seq, format = format), gq_args))
  if (!"Z"   %in% skip) results$Z   <- do.call(find_zdna, c(list(seq = seq, format = format), zdna_args))
  if (!"STR" %in% skip) results$STR <- do.call(find_str,  c(list(seq = seq, format = format), str_args))
  if (!"APR" %in% skip) results$APR <- do.call(find_apr,  c(list(seq = seq, format = format), apr_args))

  results
}
