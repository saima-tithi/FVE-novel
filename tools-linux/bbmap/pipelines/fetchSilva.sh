#!/bin/bash

#Fetches and sketches Silva.
#Be sure taxonomy is updated first!
#To use this script outside of NERSC,
#add "taxpath=/path/to/taxonomy_directory/" to each BBTools command
#Also, check Silva's archive (https://www.arb-silva.de/no_cache/download/archive/) for a newer version.


#Fetch files.
#This script only uses SSURef, but it fetches a couple others also.
#"_tax_silva_trunc" means it has taxonomic information, and it is truncated to retain only alignable 16S bases, not additional sequence.
wget -nv https://www.arb-silva.de/no_cache/download/archive/release_128/Exports/SILVA_128_SSURef_tax_silva_trunc.fasta.gz
wget -nv https://www.arb-silva.de/no_cache/download/archive/release_128/Exports/SILVA_128_LSURef_tax_silva_trunc.fasta.gz
wget -nv https://www.arb-silva.de/no_cache/download/archive/release_128/Exports/SILVA_128_SSURef_Nr99_tax_silva.fasta.gz

#This transforms U to T, trims leading or trailing Ns, discards sequences under 80bp, and changes the line wrap to 4000 to save space.
reformat.sh in=SILVA_128_SSURef_tax_silva_trunc.fasta.gz out=utot_ssu.fa.gz fastawrap=4000 utot zl=8 qtrim=rl trimq=1 ow minlen=80

#Removes all identical or fully-contained sequences.
dedupe.sh in=utot_ssu.fa.gz out=deduped_ssu.fa.gz zl=8 pigz unpigz fastawrap=4000

#Places similar sequences near each other to dramatically increase compression.
#May be possible to increase compression even more by using a blacklist, or doing a taxonomic sort, or choosing a different hash seed.
clumpify.sh in=deduped_ssu.fa.gz out=clumped_ssu.fa.gz zl=9 pigz unpigz reorder fastawrap=4000

#Make a blacklist of kmers occuring in at least 500 different species.
#A blacklist is HIGHLY recommended for Silva or any ribo database.
sketchblacklist.sh -Xmx31g in=both_clumped.fa.gz prefilter=f tree=auto taxa silva taxlevel=species ow out=blacklist_silva_species_500.sketch mincount=500

#Generate one sketch per sequence (most accurate and recommended)
time sketch.sh files=31 out=bl_ssu_seq#.sketch in=clumped_ssu.fa.gz size=200 maxgenomefraction=0.1 -Xmx8g tree=auto mode=sequence ow silva blacklist=blacklist_silva_species_500.sketch autosize ow

#...or one sketch per taxonomic unit at the subspecies level (smaller output)
time sketch.sh files=31 out=bl_ssu_taxa#.sketch in=clumped_ssu.fa.gz size=200 maxgenomefraction=0.1 -Xmx8g tree=auto mode=taxa ow silva blacklist=blacklist_silva_species_500.sketch autosize ow

