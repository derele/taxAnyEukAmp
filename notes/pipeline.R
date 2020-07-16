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


