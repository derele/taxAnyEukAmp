## https://mothur.org/wiki/silva_reference_files/

silvaSeed <- Biostrings::readDNAStringSet("/SAN/db/RDP_Silva/Silva_138/silva.seed_v138_Euk.align")

silvaSeedAl <- DECIPHER::AdjustAlignment(silvaSeed)

silRNAAl <- Biostrings::RNAStringSet(silvaSeedAl)

X18SRNA <- Biostrings::RNAStringSet(X18SCurated)

gT <- lapply(order(width(X18SRNA), decreasing=TRUE),
             function(x) {
                 attr(x, "height") <- 0
                 attr(x, "label") <- names(X18SRNA)[x]
                 attr(x, "members") <- 1L
                 attr(x, "leaf") <- TRUE
                 x
             })

attr(gT, "height") <- 0.5
attr(gT, "members") <- length(X18SRNA)
class(gT) <- "dendrogram"
# use the guide tree as inp

X18SRNAl <- DECIPHER::AlignSeqs(X18SRNA,
                                ## ## if including shorter sequences
                                ## restrict=-1e10, normPower=0,
                                iterations=0, refinements=0,
                                guideTree=gT
                                )




DECIPHER::Seqs2DB(X18SRNA, "RNAStringSet", identifier=names(X18SRNA),
                  "/SAN/db/DEC18S.sqlite",
                 tblName="input")

DECIPHER::Seqs2DB(silRNAAl, "RNAStringSet", identifier=names(silRNAAl),
                  "/SAN/db/DEC18S.sqlite",
                  tblName="silvaSeed")

DECIPHER::AlignDB("/SAN/db/DEC18S.sqlite", tblName=c("input"),
                  add2tbl="Aligned18S")


