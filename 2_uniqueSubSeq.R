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


## Use the accession as name in the  DNAStringSet
getAccessions <- function(ENAnames){
    gsub("ENA\\|.*?\\|(.*?\\.\\d+):.*", "\\1", ENAnames)
}

## rename with the accessions
names(X18SFullLength) <- getAccessions(names(X18SFullLength))

## exclude duplicates
X18SFullLength <- X18SFullLength[!duplicated(names(X18SFullLength))]

TaxIDS <- accessionToTaxa(names(X18SFullLength), "/SAN/db/taxonomy/taxonomizr.sql")

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

## now working on taxa with multiple sequences ...
## to see which ones can be excluded
nSeqTax <- table(TaxIDS[!TaxIDS%in%badTaxa])

repSeqTax <- names(nSeqTax[nSeqTax>1])

DmatList <- mclapply(repSeqTax, function (spec){
    specSeq <- X18SFullLength[TaxIDS%in%spec]
    alnSpecSeq <- AlignSeqs(specSeq)
    DistanceMatrix(alnSpecSeq)
}, mc.cores=20)

names(DmatList) <- repSeqTax

## exclude outlier sequences above a certain distance ## outliers can
## be found in matrices (above 2, 2 dimension) by a static cutoff
cutoff <- 0.1
## quant <- 0.5

outliersList <- lapply(DmatList, findOutliers, cutoff=cutoff)

## for a quartile distance based cutoff
### outliersList <- lapply(DmatList, findOutliers, quant=quant) 

## this many outliers exist for this many taxa:
table(unlist(outliersList))
with.outliers <- unlist(lapply(outliersList, any))
table(with.outliers)
## but not if are only 2x2 matrices with only two sequences not clear
## which one is the outlier
twos <- unlist(lapply(DmatList, function (x) dim(x)[[1]]==2))
table(twos, with.outliers)

### we can't do anything against the outlers in 2x2 matrices
table(twos, with.outliers)

outliers <- unlist(lapply(outliersList[!twos], function (x) names(x)[x]))

## Looking for off diagonal zeros! And then mark all but the one of
## them for removal


#### Sometimes it just takes some understanding...
## IdClusters: Cluster Sequences By Distance or Sequence
## In DECIPHER: Tools for curating, analyzing, and manipulating biological sequences 
subsequencesList <- mclapply(DmatList, function (mat) { 
    clusters <- IdClusters(mat, cutoff=0)
    clustersList <- by(clusters, clusters$cluster, rownames)
    cSubS <- lapply(clustersList, function(x) {
        targetSeq <- X18SFullLength[x]
        seqlength <- width(targetSeq)
        ## select only the first longest sequence (in case of multiple longest)
        good <- which(seqlength == max(seqlength))[1]
        one.good <- names(targetSeq)[good]
        ## blacklist of the "bad" subsequnces
        names(targetSeq)[!names(targetSeq)%in%one.good]
    })
    unlist(cSubS)
}, mc.cores=20)
      
subsequences <- unlist(subsequencesList)


foo <- data.frame(badTax=names(X18SFullLength)%in%badTaxaSeq, 
                  outlier=names(X18SFullLength)%in%outliers, 
                  subsequence=names(X18SFullLength)%in%subsequences)

bar <- as.data.frame(apply(foo, 2, as.numeric))
rownames(bar) <- 1:nrow(bar)

upset(bar)

Cleaned18S <- X18SFullLength[
    !names(X18SFullLength)%in%badTaxaSeq &
    !names(X18SFullLength)%in%outliers &
    !names(X18SFullLength)%in%subsequences
] 


