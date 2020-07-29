library(parallel)
library(DECIPHER)
library(taxonomizr)
library(clstutils)


## Remove sequences with bad characters after removing gap characters. 
removeBadSeq <- function (DNAStringSet, min.length=0){
    DNAStringSet <- RemoveGaps(DNAStringSet)
    DNAStringSet <- DNAStringSet[width(DNAStringSet)>min.length]
    goodSeq <- dada2:::C_isACGT(as.character(DNAStringSet))
    DNAStringSet[goodSeq]
}

## also make them "full length"
X18SFullLength <- removeBadSeq(X18Seq, 1700)

## or not
X18SALL <- removeBadSeq(X18Seq, 300)


## Use the accession as name in the  DNAStringSet
getAccessions <- function(ENAnames){
    gsub("ENA\\|.*?\\|(.*?\\.\\d+):.*", "\\1", ENAnames)
}

## rename with the accessions
names(X18SFullLength) <- getAccessions(names(X18SFullLength))
names(X18SALL) <- getAccessions(names(X18SALL))

## exclude duplicates
X18SFullLength <- X18SFullLength[!duplicated(names(X18SFullLength))]
X18SALL <- X18SALL[!duplicated(names(X18SALL))]

TaxIDS <- accessionToTaxa(names(X18SFullLength), "/SAN/db/taxonomy/taxonomizr.sql")
TaxIDSALL <- accessionToTaxa(names(X18SALL), "/SAN/db/taxonomy/taxonomizr.sql")

## exclude wonky taxa (uncultured and other "badly" annotated)
getTaxaBlackList <- function(taxids, maxN=10, allowedNA=3, allowSpeciesNA=FALSE){
    ntaxa <- table(taxids)
    badTax <- getTaxonomy(names(ntaxa[ntaxa > maxN]), "/SAN/db/taxonomy/taxonomizr.sql")
    badTax[grepl("uncultured", badTax)] <- NA
    badTax[grepl("sp\\.|environmental sample", badTax[, "species"]), "species"] <- NA
    is.badNA <- rowSums(is.na(badTax)) > allowedNA
    is.badSpecies <- is.na(badTax[, "species"])
    gsub(" ", "", rownames(badTax)[is.badNA|is.badSpecies])
}

## a vector of "bad" taxids 
badTaxa <- getTaxaBlackList(TaxIDS)
badTaxaSeq <- names(X18SFullLength[TaxIDS%in%badTaxa])

badTaxaALL <- getTaxaBlackList(TaxIDSALL)
badTaxaSeqALL <- names(X18SALL[TaxIDSALL%in%badTaxaALL])

## now working on taxa with multiple sequences ...
## to see which ones can be excluded

getMatrices <- function (Tids, DNAdata) { 
    nSeqTax <- table(Tids)
    repSeqTax <- names(nSeqTax[nSeqTax>1])
    DmatList <- mclapply(repSeqTax, function (spec){
        specSeq <- DNAdata[Tids%in%spec]
        alnSpecSeq <- AlignSeqs(specSeq)
        DistanceMatrix(alnSpecSeq)
    }, mc.cores=20)

    names(DmatList) <- repSeqTax
    DmatList
}

DmatList <- getMatrices(TaxIDS[!TaxIDS%in%badTaxa],
                        X18SFullLength)

DmatListALL <- getMatrices(TaxIDSALL[!TaxIDSALL%in%badTaxaALL],
                           X18SALL)

## exclude outlier sequences above a certain distance ## outliers can
## be found in matrices (above 2, 2 dimension) by a static cutoff


getOutliers <- function (matices, cutoff=0.1){
    outliersList <- lapply(DmatList, clsutils::findOutliers, cutoff=cutoff)
    ## for a quartile distance based cutoff
    ## outliersList <- lapply(DmatList, findOutliers, quant=quant) 
    ## this many outliers exist for this many taxa:
    with.outliers <- unlist(lapply(outliersList, any))
    ## but not if are only 2x2 matrices with only two sequences not clear
    ## which one is the outlier
    twos <- unlist(lapply(DmatList, function (x) dim(x)[[1]]==2))
    ## we can't do anything against the outlers in 2x2 matrices
    unlist(lapply(outliersList[!twos], function (x) names(x)[x]))
}

outliers <- getOutliers(DmatList)

outliersALL <- getOutliers(DmatListALL)

## identify sorther sequences which are a subset fo longer ones

getSubsequeces <- function(matrices, DNAdata){ 
    subsequencesList <- mclapply(matrices, function (mat) { 
        clusters <- IdClusters(mat, cutoff=0)
        clustersList <- by(clusters, clusters$cluster, rownames)
        cSubS <- lapply(clustersList, function(x) {
            targetSeq <- DNAdata[x]
            seqlength <- width(targetSeq)
            ## select only the first longest sequence (in case of multiple longest)
            good <- which(seqlength == max(seqlength))[1]
            one.good <- names(targetSeq)[good]
            ## blacklist of the "bad" subsequnces
            names(targetSeq)[!names(targetSeq)%in%one.good]
        })
        unlist(cSubS)
    }, mc.cores=20)
    unlist(subsequencesList)
}
      
subsequences <-  getSubsequeces(DmatList, X18SFullLength)

subsequencesALL <-  getSubsequeces(DmatListALL, X18SALL)




