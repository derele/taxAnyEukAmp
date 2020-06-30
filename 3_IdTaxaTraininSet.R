
Cleaned18S <- X18SFullLength[
    !names(X18SFullLength)%in%badTaxaSeq &
    !names(X18SFullLength)%in%outliers &
    !names(X18SFullLength)%in%subsequences
] 

CleanedTaxIDs <- TaxIDS[
    !names(X18SFullLength)%in%badTaxaSeq &
    !names(X18SFullLength)%in%outliers &
    !names(X18SFullLength)%in%subsequences
] 



FTaxDNAStringSet <- getTaxonomy(CleanedTaxIDs, "/SAN/db/taxonomy/taxonomizr.sql")

## specific formatting for training data
FTaxDNAStringSet[, "species"] <- gsub(".*? (.*)", "\\1",
                                      FTaxDNAStringSet[, "species"])
TaxonomyDNAStringSet <- apply(FTaxDNAStringSet, 1, paste, collapse="; ")
TaxonomyDNAStringSet <- paste("Root", TaxonomyDNAStringSet, sep="; ")

## first normalize group size for each taxon
normalizeBigGroups <- function (Taxonomy, maxGroupSize = 10) {
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
    
removeNormData <- normalizeBigGroups(TaxonomyDNAStringSet)

## remove taxa which cause problems... see ClassifySequences vignette
## of DECIPHER. We do this recursively but also exhaustively!

rmSeqs <- list()
rmSeqs[[1]] <- which(removeNormData) ## first are the results of
                                     ## normilizstion

continue <- TRUE
tSetList <- list()
i <- 1

while(continue){
    allrm <- unlist(rmSeqs)
    tSetList[[i]] <- LearnTaxa(Cleaned18S[-allrm], 
                               TaxonomyDNAStringSet[-allrm])
    rmSeqs[[i+1]] <- tSetList[[i]]$problemSequences$Index
    if(length(rmSeqs[[i+1]]) < 100) {
        continue <- FALSE
    }
    i <- i+1
}

