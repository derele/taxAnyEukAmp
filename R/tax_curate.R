##' Identify uninformative annotations ("bad" taxa) based on a taxonomy
##' table from taxonomizr
##'
##' Some sequences from databases have uninformative taxonomical
##' annotation. This function can idetnify those sequences and returns
##' a logical indicating whether the taxonomy for an entry was
##' uniformative.
##' 
##' @title getBadTaxa
##' @param taxonomy matrix of annotations from taxonomizr::getTaxonomy
##'     including a "species" column
##' @param fromN perform screening only for taxa with more abundant in
##'     the dataset than fromN (set to 0 to disable).
##' @param allowedNA the number of meaningless annotations allowed
##'     over all levels of the taxonomy
##' @param badstring a vector of uniformative taxonomical
##'     annotations. These will be treated as if no annotation was
##'     available at all.
##' @return a logical vector indicating which of the input taxa
##'     belongs to a "bad" taxon.
##' @export
##' @author Emanuel Heitlinger
getBadTaxa <- function(taxonomy, fromN=10, allowedNA=3,
                       badstring = c("uncultured", "sp\\.", "environmental", "sample",
                                     "unidentified")){
    taxonomy <- .checkTaxonomy(taxonomy)
    nSp <- table(taxonomy[, "species"])
    fromSp <- nSp[nSp > fromN]
    tnmy <- taxonomy[taxonomy[, "species"]%in%names(fromSp), ]
    ## fist set NAs
    badregex <- paste(badstring, collapse="|")
    tnmy[grepl(badregex, tnmy)] <-NA
    is.badNA <- rowSums(is.na(tnmy)) > allowedNA
    badSp <- tnmy[is.badNA, "species"]
    taxonomy[, "species"]%in%badSp
}

##' Identify uniformative species annotations based on a taxonomy
##' table from taxonomizr.
##'
##' Some sequences from databases have uninformative taxonomical
##' annotation. This function can be used to idetnify those sequences,
##' it returns a logical indicating whether the species taxonomy for
##' an entry was uniformative.
##' 
##' @title getBadSpecies
##' @param taxonomy from taxonomizr::getTaxonomy
##' @param badstring a vector of uniformative taxonomical
##'     annotations. These will be treated as if no annotation was
##'     available at all.
##' @return
##' @export
##' @author Emanuel Heitlinger
getBadSpecies <- function(taxonomy,
                       badstring = c("uncultured", "sp\\.",
                                     "environmental", "sample",
                                     "unidentified")){
    taxonomy <- .checkTaxonomy(taxonomy)
    badregex <- paste(badstring, collapse="|")
    sP <- taxonomy[, "species"]
    sP[grepl(badregex, sP)] <-NA
    is.na(sP)
}
        
