# taxAnyEukAmp

Create a taxonomy database for any (Eukaryote) amplicon

## Taxonomy annotation for any Eukaryote amplicon

This package can build a database from genetic marker sequences for
taxonomic annotation. It focuses on markers for Eukaryote taxa but can
in principle also be used for Prokaryotes. Sequence databases for to
annotate Eukaryote marker sequences with taxonomy are currently not
available in a comprehensive form. This package can be used to create
them. Typical markers employed are ribosomal RNA subunits (small and
large; nuclear and mitochondrial 18S, 28S, 12S and 16S) and Cytochrome
C oxidase subunit I (COI), for those we provide workflows and readily
constructed databases. For others (e.g. other mitochondrially encoded
genes) those workflows are directly transferable.

The package provides very high level functions to make this quite
easy. So let's jump right into it: 

## Install
```S
require(devtools)
devtools::install_github("derele/taxAnyEukAmp")
```

### Download of marker sequences for ENA

Marker sequences are downloaded from European nucleotide archive (ENA)
marker search. ENA has indexed their databases using hidden Markov
models. This allows the retrieval of marker sequences beyond the
identification given in the sequence description (e.g. also
sub-sequences from complete genomes, etc.). This function downloads
these marker sequences. 


```S
library(taxAnyEukAmp)
X18SDownloads <- getENAdownloads("18S", "/SAN/db/ENA_marker/18S/")
```
Downloaded files can be read into Biostring "DNAStringSet instances"
for subsequent work.

```S
X18Seq <- Biostrings::readDNAStringSet(X18SDownloads[[1]])
```

Before download the function also queries the database for the number
of sequences available and returns this number, so we can compare
whether all were downloaded successfully.

In case all the expected files were already in the folder given for
the download the function will only query the ENA for the number of
sequences available. Especially in this case we should check how many
sequences are available. 

```S
length(X18Seq) - X18SDownloads[["should_be"]] 
```

### Curate the database

As a basic design principle, we collect names (accession numbers) of
sequences exclution when creating the database.

#### Cleaning

Sequences with bad charactersalignemnt)

```S
isACGT(as.character(X18Seq))
```


#### Taxonomy curration

A taxonomic annotation for each accession number is obtained via the
taxnomizr R package and the NCBI taxonomy database.





A list of taxa to exclude is created by tabulating the most abundant
"species" (or pseudo-species in case of "bad annotations", setting
"environmental", "uncultured" and "sp." annotations to NA and removing
any taxa with either NA as species annotations (if the arument
speciesNArm is set to TRUE; the default), or a defined number of (by
default 3) NAs at any level of the taxonmy.

#### Sequence similarity curration

Based on alignemts of (if available) multiple sequences for each
speces: a) outlier sequences are removed if they have more than 10%
nucleotide differences with the majority cluster of sequences. Then,
b) full length sequences are made non-redundant by removing non unique
sub sequences (including potential sub-sequences) for each
species. This is performed only for identity clusters but could be
changed to remove simlilar sub-sequences at a clustering threshold

### IdTaxa (DECIPHER) curration and training set

The DECIPHER pageage is used to train datsaets as classifiers. This
allows the identification of "problem sequences".
