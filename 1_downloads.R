library(Biostrings)
library(RCurl)

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


## general function to get number of ENA released sequences
getENAsequenceCount <- function(marker, data="noncoding_release", freequery="") { 
    countURL <- "https://www.ebi.ac.uk/ena/portal/api/count"
    Copts <- curlOptions(verbose = TRUE, # to debug
                         header = "Content-Type: application/x-www-form-urlencoded")
    query=paste0(freequery, "marker%3D%22", marker, "%22")
    countRes <- RCurl::postForm(uri=countURL, 
                                .opts=Copts,
                                query=query,
                                result=data,
                                style= "POST")
    as.numeric(countRes)
}



getENAdownloads <- function(marker,
                            downloadDir,
                            data="noncoding_release",
                            paging="100000",
                            freequery=""){
    count <- getENAsequenceCount(marker, data=data, freequery=freequery)
    pages <- seq(0, count, by=100000)
    pages <- format(pages, digits=16, trim=TRUE, scientific=FALSE)
    URLs <- paste0("https://www.ebi.ac.uk/ena/browser/api/fasta/search?",
                   "query=",freequery ,
                   "marker=%22", marker,
                   "%22&result=", data,
                   "&limit=", paging,
                   "&offset=", pages)
    files <- paste0(downloadDir, marker, "_", pages, '.fasta')
    if(all(file.exists(files))){
        message("files in ", downloadDir, " exist!\n",
                "Returning names of existing files,",
                "but new URLs and new counts for the number of sequences in ENA, ",
                "in case of discrepanicies you may want to", 
                "delete files or give different downloadDir to repeat the download!\n"
                )
        return(list(files=files, URLs=URLs, should_be=count))
    }
    if (sum(file.exists(files))>0){
        stop("some files in ", downloadDir, " exist, ",
             "but don't match the number or names of files expected from a new download, ",
             "delete files or give different downloadDir to repeat the download!\n"
             )
    } else {
        for(i in seq_along(URLs)){
            download.file(URLs[[i]], files[[i]])
        }
    }
    list(files=files, URLs=URLs, should_be=count)
}

X18SDownloads <- getENAdownloads("18S", "/SAN/db/ENA_marker/18S/")

X18Seq <- Biostrings::readDNAStringSet(X18SDownloads[[1]])

## how many are recovered
length(X18Seq) - X18SDownloads[[3]] 
## we are missing 12k approximately
summary(width(X18Seq)>100)


## ##### 28S: how many released sequences on 05/27/2020

X28SDownloads <- getENAdownloads("28S", "/SAN/db/ENA_marker/28S/")
X28Seq <- Biostrings::readDNAStringSet(X28SDownloads[[1]])

## and how many sequences are recoverd
length(X28Seq) - X28SDownloads[[3]] 
summary(width(X28Seq)>100)

### COI or COX1 in ENA terms

COIDownloads <- getENAdownloads("COX1", "/SAN/db/ENA_marker/COI/", "coding_release")
COISeq <- Biostrings::readDNAStringSet(COIDownloads[[1]])

## and how many sequences are recoverd
length(COISeq) - COIDownloads[[3]] 
summary(width(COISeq)>100)

### 12S 
X12SDownloads <- getENAdownloads("12S", "/SAN/db/ENA_marker/12S/")
X12Seq <- Biostrings::readDNAStringSet(X12SDownloads[["files"]])

## how many recoverd
length(X12Seq) - X12SDownloads[[3]] 

table(width(X12Seq)>100)
table(width(X12Seq)>500)
table(width(X12Seq)<1200)
## and most are "proper"


## Eukaryote (mitochondrial) 16S, providing the eukaryote taxonomy as
## limit for the search
X16SDownloads <- getENAdownloads("16S", "/SAN/db/ENA_marker/16S_Euk/",
                                 freequery="tax_tree(2759)%20AND%20")


X16Seq <- Biostrings::readDNAStringSet(X16SDownloads[["files"]])

length(X16Seq) - X16SDownloads[["should_be"]]
## 3 sequences more than records??!!!

## and only lost one sequence!
## 263376- 263377
table(width(X16Seq)>100)
table(width(X16Seq)>500)
table(width(X16Seq)<1200)


