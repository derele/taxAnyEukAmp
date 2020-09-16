

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


##' Check for a DNAstringSet potentially within a taxedSeq object
##'
##' This internal function checks for a DNAstringSet or taxedSeq
##' object. If hte latter, it etextracts this DNAStringSet
##'
##' @title .checkSequences
##' @param DNAdata taxedSeq object containing a DNAStringSet or only
##'     the latter
##' @return A DNAStringSet
##' @author Emanuel Heitlinger
.checkSequences <- function(DNAdata){
    if("taxedSeq" %in% class(DNAdata)) {
        DNAdata <- DNAdata[["sequences"]]
    }
    if(!"DNAStringSet" %in% class(DNAdata)){
        stop("please provide a DNAStrinSet or a taxedSeq object containing one")
    }
    DNAdata
}


##' Extract and check the content of a taxonomy matrix from a taxedSeq
##' object
##' 
##' This function extracts the taxonomy table from a taxedSeq object
##' and checks whether the rquired species column is given in a
##' taxonomy table.
##' @title .checkTaxonomy
##' @param taxonomy from taxonomizr::getTaxonomy or a taxedSeq object
##'     containing this taxonomy
##' @return a taxonomy matrix
##' @author Emanuel Heitlinger
.checkTaxonomy <- function(taxonomy){
    if("taxedSeq" %in% class(taxonomy)) {
        taxonomy <- taxedSeq[["taxonomy"]]
    }
    if(!"matrix" %in% class(taxonomy)  |
       ! "species" %in% colnames(taxonomy) ){
        stop("please provide a matrix of taxonomy annotations or a taxedSeq object containing such a matrix")
    }
    taxonomy
}
