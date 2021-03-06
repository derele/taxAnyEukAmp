##' Get distance matrices from sequence aligments of each species in a dataset
##'
##' Alignments are created for each species seperately, if the species
##' has multiple sequences in the dataset
##' @title getMatrices
##' @param taxedSeq A taxedSeq object containing a DNAstingSet and an
##'     associated taxonomy
##' @param mc.cores number of threads/processor used. This works only
##'     on UNIX systems via the package parallel.
##' @return a list of distance matrices
##' @importFrom DECIPHER AlignSeqs
##' @importFrom DECIPHER DistanceMatrix
##' @export
##' @author Emanuel Heitlinger
getMatrices <- function (taxedSeq, mc.cores=1) { 
    taxonomy <- .checkTaxonomy(taxedSeq)
    DNAdata <- .checkSequences(taxedSeq)
    nSeqTax <- table(taxonomy[, "species"])
    repSeqTax <- names(nSeqTax[nSeqTax>1])
    DmatList <- mclapply(repSeqTax, function (spec){
        specSeq <- DNAdata[taxonomy[, "species"]%in%spec]
        alnSpecSeq <- AlignSeqs(specSeq)
        DistanceMatrix(alnSpecSeq)
    }, mc.cores=mc.cores)

    names(DmatList) <- repSeqTax
    DmatList
}

##' Get outlier sequences based on distance matrices for each species
##'
##' Based on distance matrices for each species (see getMatrices)
##' @title getOutliers
##' @param matrices A list of distance matrices from getMatrices
##' @param DNAdata A DNAstrinSet or a taxedSeq object containg such
##'     DNA sequence data
##' @param cutoff percent difference in the distance matrix above
##'     which a sequence is considered an outlier
##' @param taxonomy from taxonomizr::getTaxonomy
##' @return a vector strings representing names of outlier sequences
##'     in the original dataset further than cutoff distance from the
##'     mayority of sequences of the same species
##' @author Emanuel Heitlinger
getOutliers <- function (matrices, DNAdata, cutoff=0.1){
    DNAdata <- .checkSequences(DNAdata)
    outliersList <- lapply(matrices, clstutils::findOutliers, cutoff=cutoff)
    ## for a quartile distance based cutoff
    ## outliersList <- lapply(DmatList, findOutliers, quant=quant) 
    ## this many outliers exist for this many taxa:
    with.outliers <- unlist(lapply(outliersList, any))
    ## but not if are only 2x2 matrices with only two sequences not clear
    ## which one is the outlier
    twos <- unlist(lapply(matrices, function (x) dim(x)[[1]]==2))
    ## we can't do anything against the outliers in 2x2 matrices as we
    ## don't know which is the correct version
    outliers <- unlist(lapply(outliersList[!twos], function (x) names(x)[x]))
    names(DNAdata)%in%outliers
}

##' Get (shorter) subsequences relative to longer sequences for each
##' species in a taxnomically annotated dataset
##'
##' Based on distance matrices for each species (see getMatrices) this
##' function identifies the longest unique sequence and returns
##' shorter subsequencs with (perfect, when cutoff=0) matches (and
##' hence terminal gaps) in distance matrices based on alignments of
##' each species in a dataset
##' @title getSubsequeces
##' @param matrices A list of distance matrices from getMatrices
##' @param DNAdata A XString set of DNA or RNA sequences or an object
##'     of class taxedSeq
##' @param taxonomy from taxonomizr::getTaxonomy
##' @return a vector strings representing names of subsequences in the
##'     original dataset shorter than other sequences for the same
##'     species
##' @export
##' @importFrom DECIPHER IdClusters
##' @author Emanuel Heitlinger
getSubsequeces <- function(matrices, DNAdata, cutoff=0, mc.cores=1){ 
    DNAdata <- .checkSequences(DNAdata)
    subsequencesList <- lapply(matrices, function (mat) { 
        clusters <- IdClusters(mat, cutoff=cutoff)
        clustersList <- by(clusters, clusters$cluster, rownames)        
        cSubS <- mclapply(clustersList, function(x) {
            targetSeq <- DNAdata[x]
            seqlength <- width(targetSeq)
            ## select only the first longest sequence (in case of multiple longest)
            good <- which(seqlength == max(seqlength))[1]
            one.good <- names(targetSeq)[good]
            ## excludelist of the shorter subsequnces
            names(targetSeq)[!names(targetSeq)%in%one.good]
        }, mc.cores=mc.cores)
        unlist(cSubS)
    })
    subsequences <- unlist(subsequencesList)
    names(DNAdata)%in%subsequences
}



