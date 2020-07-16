# taxAnyEukAmp

Create a taxnomy database for any (Eukaryote) amplicon

## Taxonomy annotation for any Eukaryote amplicon

This package can build a datbase from genetic marker sequences for
taxomic annotation. It focuses on markers for Eukaryote taxa but can
in principle also be used for Prokaryotes. Sequence databases for to
annotate Eukaryote marker seqeunces with taxonomy are currently not
available in a comprehensive form. This package can be used to create
them. Typical markers employed are ribosomal RNA subunits (small and
large; nuclear and mitochondrial 18S, 28S, 12S and 16S) and Cytochrome
C oxidase subunit I (COI), for those we provide workflows and readiely
constructed databases. For others (e.g. other mitochondrially encoded
genes) those workflows are directly transferrable.

The package provides very high level functions to make this quite
easy. So let's jump right into it: 

## Install
```r
require(devtools)
devtools::install_github("derele/taxAnyEukAmp")
```

### Download of maker seqeunces for ENA

Marker sequences are downloaded from European nucleotide archive (ENA)
marker search. ENA has indexed their databases using hidden markov
models. This allows the retrieval of marker sequences beyond the
identification given in the sequence description (e.g. also
subsequences from complete genomes, etc.). This function downloads
these marker sequences. 

```r
X18SDownloads <- getENAdownloads("18S", "/SAN/db/ENA_marker/18S/")
```
Downloaded files can be read into Biostring "DNAStringSet instances"
for subsequent work.

```r
X18Seq <- Biostrings::readDNAStringSet(X18SDownloads[[1]])
```

Before download the function also queries the database for the number
of sequends available and returns this number, so we can compare
whether all were downloaded sucessfully

```r
length(X18Seq) - X18SDownloads[[3]] ## we
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