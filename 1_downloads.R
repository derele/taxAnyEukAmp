library(Biostrings)


## ITS1 search
## http://itsonedb.cloud.ba.infn.it/

## ## and the download from the link below for ITS1
## wget -O markers_ITS1_ITSoneDB.fasta.gz -X GET http://itsonedb.cloud.ba.infn.it/ExportController?dbexport=ENA

## SHOULD a new download be performed?
download <- c(
    ##    "18S",
    ##    "COI",
    ##    "28S",
    ##    "12S",
    ##    "16S",
    FALSE
    )

## We do this programmatically via THE NEW INTERFACE to ENA (2020)
## marker search FASTA, e.g.:
## "https://www.ebi.ac.uk/ena/browser/api/fasta/search?query=marker=%2218S%22&result=noncoding_release&fields=accession,marker,rna_class&limit=100"



## ##### 18S: how many released sequences on 05/27/2020

## curl -X POST -H "Content-Type: application/x-www-form-urlencoded" -d 'result=noncoding_release&query=marker%3D%2218S%22&limit=1000000000&format=tsv' "https://www.ebi.ac.uk/ena/portal/api/search" | wc -l

## 1170464
## because of paging we have to get chunks of 100,000 sequences
pages18S <- seq(0, 1170464, by=100000)
pages18S <- format(pages18S, digits=16, trim=TRUE, scientific=FALSE)

X18SURLs <- paste0("https://www.ebi.ac.uk/ena/browser/api/fasta/search?query=marker=%2218S%22&result=noncoding_release&limit=100000&offset=", pages18S)

X18SFiles <- paste0('/SAN/db/ENA_marker/18S/', pages18S, '.fasta')

commands18S <- paste0("wget -O ", X18SFiles, " -X GET ", "\'", X18SURLs, "\'")

if("18S"%in%download){
    resultsList18S <- list()
    for(i in seq_along(commands18S)){
        resultsList18S[[i]] <- system(commands18S[[i]])
    }
    names(resultsList18S) <- pages18S
    table(unlist(resultsList18S)==0)
}

X18Seq <- Biostrings::readDNAStringSet(X18SFiles)
## "lost" only about 13,000 records 

## and many sequences are recoverd
summary(width(X18Seq)>100)


## ##### 28S: how many released sequences on 05/27/2020

## curl -X POST -H "Content-Type: application/x-www-form-urlencoded" -d 'result=noncoding_release&query=marker%3D%2228S%22&limit=1000000000&format=tsv' "https://www.ebi.ac.uk/ena/portal/api/search" | wc -l

## Only 120474
## because of paging we have to get chunks of 100,000 sequences
pages28S <- seq(0, 120474, by=100000)
pages28S <- format(pages28S, digits=16, trim=TRUE, scientific=FALSE)

X28SURLs <- paste0("https://www.ebi.ac.uk/ena/browser/api/fasta/search?query=marker=%2228S%22&result=noncoding_release&limit=100000&offset=", pages28S)

X28SFiles <- paste0('/SAN/db/ENA_marker/28S/', pages28S, '.fasta')

commands28S <- paste0("wget -O ", X28SFiles, " -X GET ", "\'", X28SURLs, "\'")


if("28S"%in%download){
    resultsList28S <- list()
    for(i in seq_along(commands28S)){
        resultsList28S[[i]] <- system(commands28S[[i]])
    }
    names(resultsList28S) <- pages28S
    table(unlist(resultsList28S)==0)
}

X28Seq <- Biostrings::readDNAStringSet(X28SFiles)
## "lost" only about 13,000 records 

## and many sequences are recoverd
summary(width(X28Seq)>100)

## ##### COI: how many released sequences on 05/27/2020

## curl -X POST -H "Content-Type: application/x-www-form-urlencoded" -d 'result=coding_release&query=marker%3D%22COX1%22&limit=1000000000&format=tsv' "https://www.ebi.ac.uk/ena/portal/api/search" | wc -l

## 3481122

pagesCOI <- seq(0, 3481122, by=100000)
pagesCOI <- format(pagesCOI, digits=16, trim=TRUE)

COIURLs <- paste0("https://www.ebi.ac.uk/ena/browser/api/fasta/search?query=marker=%22COX1%22&result=coding_release&limit=100000&offset=", pagesCOI)

COIFiles <- paste0('/SAN/db/ENA_marker/COI/', pagesCOI, '.fasta')

commandsCOI <- paste0("wget -O ", COIFiles, " -X GET ", "\'", COIURLs, "\'")

if("COI"%in%download){
    resultsListCOI <- list()
    for(i in seq_along(commandsCOI)){
        resultsListCOI[[i]] <- system(commandsCOI[[i]])
    }
    table(unlist(resultsListCOI)==0)
}

COISeq <- Biostrings::readDNAStringSet(COIFiles)
## 3 sequences more than records??!!!

## and many sequences are recoverd
summary(width(COISeq)>100)
summary(width(COISeq)>500)
## and most are "proper"


## ##### 12S: how many released sequences on 05/28/2020

## curl -X POST -H "Content-Type: application/x-www-form-urlencoded" -d 'result=noncoding_release&query=marker%3D%2212S%22&limit=1000000000&format=tsv' "https://www.ebi.ac.uk/ena/portal/api/search" | wc -l

## 263377

pages12S <- seq(0, 263377, by=100000)
pages12S <- format(pages12S, digits=16, scientific=FALSE, trim=TRUE)

X12SURLs <- paste0("https://www.ebi.ac.uk/ena/browser/api/fasta/search?query=marker=%2212S%22&result=noncoding_release&limit=100000&offset=", pages12S)

X12SFiles <- paste0('/SAN/db/ENA_marker/12S/', pages12S, '.fasta')

commands12S <- paste0("wget -O ", X12SFiles, " -X GET ", "\'", X12SURLs, "\'")


if("12S"%in%download){
    resultsList12S <- list()
    for(i in seq_along(commands12S)){
        resultsList12S[[i]] <- system(commands12S[[i]])
    }

    table(unlist(resultsList12S)==0)
}

X12Seq <- Biostrings::readDNAStringSet(X12SFiles)
## 3 sequences more than records??!!!

## and only lost one sequence!
## 263376- 263377
table(width(X12Seq)>100)
table(width(X12Seq)>500)
table(width(X12Seq)<1200)
## and most are "proper"


## Eukaryote (mitochondrial) 16S

## curl -X POST -H "Content-Type: application/x-www-form-urlencoded" -d 'result=noncoding_release&query=tax_tree(2759)%20AND%20marker%3D%2216S%22&limit=1000000000&format=tsv' "https://www.ebi.ac.uk/ena/portal/api/search" | wc -l

## 489567

pages16S <- seq(0, 489567, by=100000)
pages16S <- format(pages16S, digits=16, scientific=FALSE, trim=TRUE)

X16SURLs <- paste0("https://www.ebi.ac.uk/ena/browser/api/fasta/search?query=tax_tree(2759)%20AND%20marker=%2216S%22&result=noncoding_release&limit=100000&offset=", pages16S)

X16SFiles <- paste0('/SAN/db/ENA_marker/16S_Euk/', pages16S, '.fasta')

commands16S <- paste0("wget -O ", X16SFiles, " -X GET ", "\'", X16SURLs, "\'")


if("16S"%in%download){
    resultsList16S <- list()
    for(i in seq_along(commands16S)){
        resultsList16S[[i]] <- system(commands16S[[i]])
    }
    table(unlist(resultsList16S)==0)
}

X16Seq <- Biostrings::readDNAStringSet(X16SFiles)
## 3 sequences more than records??!!!

## and only lost one sequence!
## 263376- 263377
table(width(X16Seq)>100)
table(width(X16Seq)>500)
table(width(X16Seq)<1200)


