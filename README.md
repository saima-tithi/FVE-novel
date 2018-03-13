# FastViromeExplorer-novel
Recover draft genomes of viruses and phages through reference based mapping using FastViromeExplorer and then iterative assembly using Spades.

# Installation
FastViromeExplorer requires the following tools installed in the user's machine.
1. FastViromeExplorer (JAVA (JDK) 1.8 or later, Samtools 1.4 or later, and Kallisto 0.43.0 or later installed)
2. bbmap
3. Spades
4. Salmon
5. CD-HIT
6. Bedtools
7. Mummer

# Run
Download the reference database and kallisto index file by downloading "GOV-viral-populations" folder from http://bench.cs.vt.edu/FastViromeExplorer/. Then download 2 read files containing ocean metagenome from NCBI SRA under the accession number SRR5677466 (https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR5677466).

Run FastViromeExplorer using the following commands:
```bash
mkdir FVE-outputDirectory
java -cp bin FastViromeExplorer -1 $read1File -2 $read2File -i /path-to-referencedb-folder/GOV-viral-populations/GOV_viral_populations.idx -l /path-to-referencedb-folder/GOV-viral-populations/gov_viral_populations-length.txt -o FVE-outputDirectory
```

Then run FastViromeExplorer-novel using the following commands:
```bash
mkdir FVENovel-outputDirectory
java -cp bin FVENovel -1 $read1File -2 $read2File -o FVENovel-outputDirectory -fveres /path-to-FVE-res/FVE-outputDirectory -dbType gov -dbDir /path-to-referencedb-folder/GOV-viral-populations -topBins 2
```

# Output
The ouput folders and files will be generated in the `FVENovel-outputDirectory` folder. The outputs will be stored in `final-results` folder inside `FVENovel-outputDirectory` folder.

The output files are:
1. *final-scaffolds.fasta* : all scaffolds assembled from the reads, they are draft viral genomes
2. *bin-info.tsv* : For each scaffold, gives the bin information of that scaffold.
3. *ANI-info.tsv* : For each scaffold, gives the ANI(Average Nucleotide Identity) of this scaffold with the reference genome (and all other reference genomes present in the same bin) from where this scaffold got assembled.  