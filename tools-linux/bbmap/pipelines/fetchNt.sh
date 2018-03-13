#!/bin/bash

#Fetches and sketches nt.
#Be sure taxonomy is updated first!
#To use this script outside of NERSC,
#add "taxpath=/path/to/taxonomy_directory/" to each BBTools command

#Ensure necessary executables are in your path
#module load bbtools
module load pigz

#Fetch nt.
wget -nv ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nt.gz

#Rename it because NCBI named it incorrectly.
mv nt.gz nt.fa.gz

#Rename by taxID by looking up gi numbers or accessions
time gi2taxid.sh -Xmx63g in=nt.fa.gz out=renamed.fa.gz pigz=32 zl=8 tree=auto accession=auto table=auto ow shrinknames

#Sort by taxonomy.
#This makes sketching by taxa use much less memory because sketches can be written to disk as soon as they are finished.
time sortbyname.sh -Xmx63g in=renamed.fa.gz out=sorted.fa.gz ow taxa tree=auto gi=ignore fastawrap=255 zl=8 pigz=32 minlen=60

#Make a blacklist of kmers occuring in at least 1000 different species.
time sketchblacklist.sh -Xmx63g in=sorted.fa.gz prepasses=1 tree=auto taxa taxlevel=species ow out=blacklist_nt_species_1000.sketch mincount=1000 k=31,24

#Generate 31 sketch files, with one sketch per species.
#Multiple files allow sketches to load faster on multicore systems. 
time bbsketch.sh -Xmx63g in=sorted.fa.gz out=taxa#.sketch mode=taxa tree=auto accession=null gi=null files=31 ow unpigz minsize=300 prefilter autosize blacklist=blacklist_nt_species_1000.sketch k=31,24 depth

#A query such as contigs.fa can now be compared to the new reference sketches like this:
#comparesketch.sh in=contigs.fa k=31,24 tree=auto taxa*.sketch blacklist=blacklist_nt_species_1000.sketch

#On NERSC systems, you can then set the default path to nt by pointing /global/projectb/sandbox/gaag/bbtools/nt/current at the path to the new sketches.
#Then you can use the default set of nt sketches like this:
#comparesketch.sh in=contigs.fa nt tree=auto
#That command automatically adds the default path to the sketches, the blacklist, and the correct values for K.
