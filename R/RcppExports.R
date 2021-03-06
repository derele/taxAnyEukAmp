# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Checks a vector of character sequences for whether they are entirely ACGT.
#'
#' @param seqs A \code{character} of candidate DNA sequences.
#' @return A \code{logical}. Whether or not each input character was ACGT only.
#' @export
isACGT <- function(seqs) {
    .Call(`_taxAnyEukAmp_isACGT`, seqs)
}

