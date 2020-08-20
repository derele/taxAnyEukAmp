##' Recursively train a IdTaxa classifier
##'
##' This function is a wrapper around DECIPER's Learn Taxa. It
##' recursively trains the classifier removing problem sequences in
##' each iteration.
##' @title idTaxaTrainAndClean
##' @param DNAdata DNAStringSet of sequences to be used for training
##' @param TAXdata the taxonomic annotation in a dataframe as obtained
##'     from the taxonomizr package
##' @param normalizeSize integer indicating to use only this number of
##'     sequences in case of mutliple sequnces for the same taxon
##' @param badRetain postiv integer indicating to repeat classifier
##'     training until the number of "problem sequences" is lower than
##'     this number
##' @return A list with three elements: [1] the sequence the latest
##'     classifier is based on [2] the taxonomy for the sequences in
##'     the latest classifier [3] the latest classifier itself
##' @importFrom DECIPHER LearnTaxa
##' @export
##' @author Emanuel Heitlinger
idTaxaTrainAndClean <- function (DNAdata, TAXdata,
                                 normalizeSize=10,
                                 badRetain=100
                                 ) {


    ## specific formatting for training data
    TAXdata[, "species"] <- gsub(".*? (.*)", "\\1",
                                          TAXdata[, "species"])
    TaxonomyDNAStringSet <- apply(TAXdata, 1, paste, collapse="; ")
    TaxonomyDNAStringSet <- paste("Root", TaxonomyDNAStringSet, sep="; ")

    ## first normalize group size for each taxon
    normalizeBigGroups <- function (Taxonomy, maxGroupSize) {
        groupCounts <- table(Taxonomy)        
        remove <- logical(length(Taxonomy))
        for (i in which(groupCounts > maxGroupSize)) {
            index <- which(Taxonomy==names(groupCounts)[i])
            keep <- sample(length(index),
                           maxGroupSize)
            remove[index[-keep]] <- TRUE
        }
        remove
    }
    
    removeNormData <- normalizeBigGroups(TaxonomyDNAStringSet,
                                         maxGroupSize = normalizeSize)

## remove taxa which cause problems... see ClassifySequences vignette
## of DECIPHER. We do this recursively but also exhaustively!

    rmSeqs <- list()
    rmSeqs[[1]] <- which(removeNormData) ## first to remove are the
    ## results of normilizstion

    continue <- TRUE
    tSetList <- list()
    i <- 1

    while(continue){
        allrm <- unlist(rmSeqs)
        tSetList[[i]] <- LearnTaxa(DNAdata[-allrm], 
                                   TaxonomyDNAStringSet[-allrm])
        ## make the new indices relative to the overall dataset
        orig.idx.before <- seq_len(length(DNAdata))[-allrm]
        orig.idx.to.rm <- orig.idx.before[tSetList[[i]]$problemSequences$Index]
        rmSeqs[[i+1]] <- orig.idx.to.rm
        if(length(rmSeqs[[i+1]]) < badRetain) {
            continue <- FALSE
        }
        i <- i+1
        cat("\nrep", i, "removing", length(rmSeqs[[i]]), "sequences", 
            "\nhead:", head(rmSeqs[[i]]), "\ntail:", tail(rmSeqs[[i]]), "\n")
    }

    list(DNAdata[-allrm], ## the sequence the latest classifier is based on
         TAXdata[-allrm], ## the taxonomy for the sequences in the latest classifier
         ## the latest classifier itself
         tSetList[[length(tSetList)]])
}
    
