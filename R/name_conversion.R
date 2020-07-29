##' Convert ENA names to (NCBI) accession names
##'
##' Extract the NCBI accession name for sequences from the name
##' European Nucleotide Archive gives to downloaded (fasta) sequences
##' @title getAccession4ENAname
##' @param ENAnames the full name of the (fasta) sequence in the ENA
##'     download
##' @return A string with the (NCBI) accession for a sequence
##' @export
##' @author Emanuel Heitlinger
getAccession4ENAname <- function(ENAnames){
    gsub("ENA\\|.*?\\|(.*?\\.\\d+):.*", "\\1", ENAnames)
}
