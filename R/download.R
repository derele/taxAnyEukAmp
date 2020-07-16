##' Get the number of ENA released sequences for a marker
##'
##' Post an http request to ENA to count the number of available
##' sequences for a marker
##' @title .getENAsequenceCount
##' @param marker The marker to be queried see https://www.ebi.ac.uk/ena/browse/marker-portal-rest
##' @param data the database to be queried usually "noncoding_release", "coding_release" see https://www.ebi.ac.uk/ena/browser/advanced-search
##' @param freequery terms in addtion to marker to be included included in the query (e.g. to limit the scope of the search to particular nodes in the taxonomy)
##' @return the number of sequences found for this marker search
##' @author Emanuel Heitlinger
##' @importFrom RCurl curlOptions
.getENAsequenceCount <- function(marker, data="noncoding_release", freequery="") { 
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


##' Download marker sequences from ENA. 
##'
##' ENA has indexed their databases using hidden markov models. This allows the retrieval of marker sequences beyond what is given in the sequence description (e.g. also subsequences from complete genomes, etc.). This function downloads these marker sequences
##' @title getENAdownloads
##' @param marker The marker to be queried see https://www.ebi.ac.uk/ena/browse/marker-portal-rest
##' @param downloadDir Directory into which the downloaded sequences are saved 
##' @param data the database to be queried usually "noncoding_release", "coding_release" see https://www.ebi.ac.uk/ena/browser/advanced-search
##' @param paging Downloads are broken down in into chunks. This is the chunk-size for this pageing.
##' @param freequery terms in addtion to marker to be included included in the query (e.g. to limit the scope of the search to particular nodes in the taxonomy)
##' @return Returns a list containing filenames (files), the URLs the files were downloaded from (URLs) and the number of expected sequences in the database that should be downloaded (should_be)
##' @author Emanuel Heitlinger
##' @importFrom utils download.file
##' @export
getENAdownloads <- function(marker,
                            downloadDir,
                            data="noncoding_release",
                            paging="100000",
                            freequery=""){
    count <- .getENAsequenceCount(marker, data=data, freequery=freequery)
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

