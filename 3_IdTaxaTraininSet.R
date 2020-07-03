library(ShortRead)

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
rmSeqs[[1]] <- which(removeNormData) ## first to remove are the
                                     ## results of normilizstion

continue <- TRUE
tSetList <- list()
i <- 1

while(continue){
    allrm <- unlist(rmSeqs)
    tSetList[[i]] <- LearnTaxa(Cleaned18S[-allrm], 
                               TaxonomyDNAStringSet[-allrm])
    ## make the new indices relative to the overall dataset
    orig.idx.before <- seq_len(length(Cleaned18S))[-allrm]
    orig.idx.to.rm <- orig.idx.before[tSetList[[i]]$problemSequences$Index]
    rmSeqs[[i+1]] <- orig.idx.to.rm
    if(length(rmSeqs[[i+1]]) < 10) {
        continue <- FALSE
    }
    i <- i+1
    cat("\nrep", i, "removing", length(rmSeqs[[i]]), "sequences", 
        "\nhead:", head(rmSeqs[[i]]), "\ntail:", tail(rmSeqs[[i]]), "\n")
}

## that's the latest classifier
tSetList[[16]]

length(allrm)

## these last 9 we didn't add to the removal list
## rmSeqs[[17]]

## allrmfinal <- unlist(rmSeqs)

finalFL18S <- Cleaned18S[-allrm]

writeFasta(finalFL18S, "/SAN/db/blastdb/18S_ENA/Full_length_1700.fasta")



