# taxAnyEukAmp
## Taxonomy annotation for any Eukaryote amplicon

Marker sequences are downloaded from ENA marker search and downloaded
files read into Biostring "DNAStringSet instances" for subsequent
work.

```r
X18SDownloads <- getENAdownloads("18S", "/SAN/db/ENA_marker/18S/")

X18Seq <- Biostrings::readDNAStringSet(X18SDownloads[[1]])

## how many are recovered
length(X18Seq) - X18SDownloads[[3]] 
## we are missing 12k approximately
summary(width(X18Seq)>100)

```





- "2_uniqueSubSeq.R":

First a length threshold is specified allowing to either length filter
or even create full length databases. A taxonomic annotation for each
accession number is obtained via the taxnomizr R package and the NCBI
taxonomy database.

A list of taxa to exclude is created by tabulating the most abundant
"species", setting "environmental", "uncultured" and "sp." annotations
to NA and removing any taxa with eigher NA as species annotation or a
defined number of (currently 3) NAs at any level of the taxonmy.

Based on alignemts of (if available) multiple sequences for each
speces: a) outlier sequences are removed if they have more than 10%
nucleotide differences with the majority cluster of sequences. Then,
b) full length sequences are made non-redundant by removing non unique
sub sequences (including potential sub-sequences) for each
species. This is performed only for identity clusters but could be
changed to remove simlilar sub-sequences at a clustering threshold

- "3_IdTaxaTraininSet.R"

The DECIPHER pageage is used to train datsaets as classifiers. This
allows the identification of "problem sequences".


- ... run pipeline... to add here... 