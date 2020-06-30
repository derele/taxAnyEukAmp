# taxAnyEukAmp
## Taxonomy annotation for any Eukaryote amplicon

- "1_downloads.R": Marker sequences are downloaded from ENA marker
   search and downloaded files read into Biostring "DNAStringSet
   instances" for subsequent work.

- "2_uniqueSubSeq.R": Sequences are made non-redundant by removing non
  unique sequences (including potential sub-sequences)