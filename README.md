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

Sequences with bad characters should be marked for removal at this
point.

```S
nonACGT <- !isACGT(as.character(X18Seq))
```

#### Taxonomy curration

A taxonomic annotation for each accession number is obtained via the
taxnomizr R package from the NCBI taxonomy database. We first extract
the accessions numbers from the names of the sequence object and then
query a [taxonomizr](https://github.com/sherrillmix/taxonomizr)
database (see the link on how to construct this on you rown system) to
obtain the taxonomy idntifiers. This then allows us to highlight
duplicated sequences for the same taxon.

```S
accessions <- getAccession4ENAname(names(X18Seq))
TaxIDs <- taxonomizr::accessionToTaxa(accessions, "/SAN/db/taxonomy/taxonomizr.sql")
duplicates <- duplicated(X18Seq)&duplicated(TaxIDs)           
```

Then we can obtain the full taxonomy path for all those taxids (again
through the taxonomizr package). By default this comprises the ranks
"superkingdom", "phylum", "class", "order", "family", "genus" and
"species" levels. We identify "bad taxa" based on this.


```S
taxonomy <- taxonomizr::getTaxonomy(TaxIDs, "/SAN/db/taxonomy/taxonomizr.sql")
rownames(taxonomy) <- names(X18Seq)
         
badTaxa <- getBadTaxa(taxonomy, fromN=2)
badSpecies <- getBadSpecies(taxonomy)         
table(badTaxa, badSpecies)
```

Taxa to (potentialy) exclude are compiled created: 1) Allowing only
for a maximal number of undefined ranks, optionally after declaring
"environmental", "uncultured" and "sp."  (pseudo-) annotations as
undfined (setting them to NA) and optionally only for the most
abundant "species" (or pseudo-species in case of "bad
annotations"). This is done by ```getBadTaxa```. 2) Highlighting taxa
with undefined species annotations, again optionall after removing
pseudo-annotations with ```getBadSpecies```.

Which taxa are in this way excluded depends on the purpose of the
created database. I generally tend to be more inclusive, i.e. to
exclude only the most uniformative and most abundant pseudo-taxa.


After all this somewhat harder tasks we can in the next steps also
simply decide to only use sequneces of a certain length (i.e. when
building a [close to] full length database).

```S
X18SClean <- X18Seq[!nonACGT & !duplicates & !badTaxa & !badSpecies &
                    width(X18Seq)>1500] ## this excludes below 1500 bases

taxonomyClean <- taxonomy[!nonACGT & !duplicates & !badTaxa & !badSpecies
                          & width(X18Seq)>1500, ]  
```

#### Sequence similarity curration

Taxonomic annotation databases are more easiely usable if they have a
relatively low redundancy.  Based on alignemts of (if available)
multiple sequences for each species, thus: a) outlier sequences are
removed if they have more than a particular threshold of nucleotide
differences (default 10%) with the majority cluster of
sequences. Then, b) full length sequences are made non-redundant by
removing non unique (sub-)sequences for each species. This is
performed only for (100%) identity clusters but could be changed to
remove simlilar sub-sequences at a clustering threshold


```S
X18SMatrices <- getMatrices(X18SClean, taxonomyClean, mc.cores=20)

outliers <- getOutliers(X18SMatrices)

subsequences <- getSubsequeces(X18SMatrices, X18SClean, mc.cores=20)

table(outliers=names(X18SClean)%in%outliers,
      subseq=names(X18SClean)%in%subsequences)

X18SCurated <- X18SClean[!names(X18SClean)%in%outliers &
                         !names(X18SClean)%in%subsequences]

taxonomyCurated <- taxonomyClean[!rownames(taxonomyClean)%in%outliers &
                                 !rownames(taxonomyClean)%in%subsequences, ]

                         
```

Not that (in contrast to the procedure above, in which we just
obtained an index for sequences to remove) we obtain the sequence
names here. This is to make sure that those sequences can be removed
irrespective of the order of sequences (i.e. from the original and/or
the otherwise cleand datasets).


### IdTaxa (DECIPHER) curration and training set

The DECIPHER pageage is used to train datsaets as classifiers. This
allows the identification of "problem sequences". The wrapper function
below recursivel removes such problematic sequences until only cerain
number ("badRetain") are left in the dataset.

The function also "normalizes" the number of sequences for taxa with a
large number of sequences to a certain number (normalizeSize; default
10). If many such redundant sequenes for the same taxon have to be
reomved consider screening your input data more stringently for
subsequences, "bad taxa" or short sequences. 


```
IdTaxaResults <- idTaxaTrainAndClean(X18SCurated, 
                                     taxonomyCurated,
                                     badRetain=10)
```

Your final training data set is now found in ```IdTaxaResults[[3]]```
the corresponding sequences in ```IdTaxaResults[[3]]``` can be writen
to a fasta file to use them in other classification methods
(e.g. BLAST, RDP, etc...).

```
writeFasta(IdTaxaResults[[1]], "/SAN/db/blastdb/18S_ENA/Full_length_1700.fasta")
```
