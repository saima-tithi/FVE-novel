#!/bin/bash

#Fetches and sketches RefSeq.
#Be sure taxonomy is updated first!
#To use this script outside of NERSC,
#add "taxpath=/path/to/taxonomy_directory/" to each BBTools command

#Ensure necessary executables are in your path
#module load bbtools
module load pigz

#Fetch RefSeq
time wget -nv ftp://ftp.ncbi.nlm.nih.gov/refseq/release/complete/*genomic.fna.gz

#Concatenate into a single file
time cat *genomic.fna.gz > all.fa.gz

#Optionally, delete the old files
#rm *genomic.fna.gz

#Rename by taxID by looking up gi numbers or accessions
time gi2taxid.sh -Xmx63g in=all.fa.gz out=renamed.fa.gz tree=auto table=auto accession=auto zl=6

#Sort by taxonomy.
#This makes sketching by taxa use much less memory because sketches can be written to disk as soon as they are finished.
time sortbyname.sh -Xmx63g in=renamed.fa.gz out=sorted.fa.gz zl=6 pigz=32 taxa tree=auto gi=ignore fastawrap=255 minlen=60

#Make a blacklist of kmers occuring in at least 300 different species.
time sketchblacklist.sh -Xmx63g in=sorted.fa.gz prepasses=1 tree=auto taxa taxlevel=species ow out=blacklist_refseq_species_300.sketch mincount=300 k=31,24

#Generate 31 sketch files, with one sketch per species.
time bbsketch.sh -Xmx63g in=sorted.fa.gz out=taxa#.sketch mode=taxa tree=auto accession=null gi=null files=31 ow unpigz minsize=400 prefilter autosize blacklist=blacklist_refseq_species_300.sketch k=31,24 depth

#A query such as contigs.fa can now be compared to the new reference sketches like this:
#comparesketch.sh in=contigs.fa k=31,24 tree=auto taxa*.sketch blacklist=blacklist_refseq_species_300.sketch

#On NERSC systems, you can then set the default path to nt by pointing /global/projectb/sandbox/gaag/bbtools/refseq/current at the path to the new sketches.
#Then you can use the default set of nt sketches like this:
#comparesketch.sh in=contigs.fa refseq tree=auto
#That command automatically adds the default path to the sketches, the blacklist, and the correct values for K.
