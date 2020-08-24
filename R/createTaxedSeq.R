

##' Create an object of class taxedSeq
##'
##' An object of class taxedSeq can be compiled from sequences,
##' matching taxiIDs and a taxonomy data-frame, usually obtained from
##' the taxnomomizr package. This object class makes sure that the
##' sequence and contents can be subsetted in a correct synchronous
##' manner.
##' @title createTaxedSeq
##' @param sequences DNAStringSet object containing the DNA
##' @param taxIDs NCBI taxonomy IDs (usually from the taxonomizr
##'     package)
##' @param taxonomy taxonomy datafraome (usually from the taxonomizr
##'     package)
##' @return a taxedSeq object
##' @author Emanuel Heitlinger
createTaxedSeq <- function(sequences, taxIDs, taxonomy) {
    if(!class(sequences)%in%"DNAStringSet"){
        stop("sequence must be provided as DNAStringSet object")
    }
    if(!length(taxIDs)==nrow(taxonomy)){
        stop("taxIDs and taxonomy are incongruent")
    }
    taxedSeq <- list(sequences=sequences, taxIDs=taxIDs, taxonomy=taxonomy)
    class(taxedSeq) <- append(class(taxedSeq), "taxedSeq")
    return(taxedSeq)
}

##' Print the contents of a taxedSeq object
##'
##' Prints the DNAStingSet using the default and head an tail of the taxonomy
##' @title print.taxedSeq
##' @param taxedSeq an object of class taxedSeq to be printed
##' @return printing to the console
##' @export
##' @author Emanuel Heitlinger
print.taxedSeq <- function (taxedSeq){
    print (taxedSeq[["sequences"]])
    print(head(taxedSeq[["taxonomy"]], n=4))
    cat( "... ", nrow(taxedSeq[["taxonomy"]]) - 8, " more rows ... \n")
    print(tail(taxedSeq[["taxonomy"]], n=4))
}
