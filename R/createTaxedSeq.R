

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
##' @export
##' @author Emanuel Heitlinger
createTaxedSeq <- function(sequences, taxonomy) {
    if(!class(sequences)%in%"DNAStringSet"){
        stop("sequence must be provided as DNAStringSet object")
    }
    if(!length(sequences)==nrow(taxonomy)){
        stop("sequences and taxonomy are incongruent")
    }
    taxedSeq <- list(sequences=sequences, taxonomy=taxonomy)
    class(taxedSeq) <- "taxedSeq"
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
    cat("class", class(taxedSeq), "object, consiting of:\n")
    print (taxedSeq[["sequences"]])
    cat("taxonomy data frame:\n")
    if (nrow(taxedSeq[["taxonomy"]])>10){
        print(head(taxedSeq[["taxonomy"]], n=5))
        cat( "... ", nrow(taxedSeq[["taxonomy"]]) - 10, " more rows ... \n")
        print(tail(taxedSeq[["taxonomy"]], n=5))
    } else{print(taxedSeq[["taxonomy"]])}
}


##' Subset taxedSeq objects
##'
##' Subsetting an object of class taxedSeq
##' @title 
##' @param x object of class taxedSeq
##' @param i index usedf for subsetting
##' @param taxedSeq an object of class taxedSeq to be subsetted
##' @return subsetted object of class taxedSeq
##' @author Emanuel Heitlinger
subselectTaxedSeq <- function (x, i){
    createTaxedSeq(sequences = x[["sequences"]][i],
                   taxonomy = x[["taxonomy"]][i, ])
}


