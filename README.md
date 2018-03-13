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

## Install FastViromeExplorer
Install FastViromeExplorer from here: https://code.vt.edu/saima5/FastViromeExplorer
 
## Download FastViromeExplorer-novel
You can download FastViromeExplorer-novel directly from github and extract it. You can also download it using the following command:
```bash
git clone https://github.com/saima-tithi/FVE-novel.git
```
From now on, we will refer the FastViromeExplorer-novel directory in the user's local machine as `project directory`. The `project directory` will contain 3 folders: src, bin, and tools-linux.

Install bbmap:
```bash 
export PATH=$PATH:/path-to-project-directory/tools-linux/bbmap
```
Install Spades:
```bash
export PATH=$PATH:/path-to-project-directory/tools-linux/SPAdes-3.10.1-Linux/bin
```
Install Salmon:
```bash
export PATH=$PATH:/path-to-project-directory/tools-linux/Salmon-latest_linux_x86_64/bin
```
Install CD-HIT:
```bash
export PATH=$PATH:/path-to-project-directory/tools-linux/cd-hit-v4.6.8-2017-0621
```
Install Bedtools:
```bash
export PATH=$PATH:/path-to-project-directory/tools-linux/bedtools2/bin
```
Install Mummer:
```bash
export PATH=$PATH:/path-to-project-directory/tools-linux/MUMmer3.23
```

# Run FastViromeExplorer-novel
Download the reference database and kallisto index file by downloading "GOV-viral-populations" folder from http://bench.cs.vt.edu/FastViromeExplorer/. Then download 2 read files containing ocean metagenome from NCBI SRA under the accession number SRR5677466 (https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR5677466). You can download the read files using the following command:
```bash
fastq-dump --skip-technical  --readids --read-filter pass --dumpbase --split-files --clip SRR5677466
```
Two paired-end read files, SRR5677466_pass_1.fastq and SRR5677466_pass_2.fastq will be downloaded.
Run FastViromeExplorer using the following commands:
```bash
mkdir /path-to-results/FVE-outputDirectory
cd FVE-project-directory
javac -d bin src/*.java
java -cp bin FastViromeExplorer -1 $read1File -2 $read2File -i /path-to-referencedb-folder/GOV-viral-populations/GOV_viral_populations.idx -l /path-to-referencedb-folder/GOV-viral-populations/gov_viral_populations-length.txt -o /path-to-results/FVE-outputDirectory
```

Then run FastViromeExplorer-novel using the following commands:
```bash
mkdir /path-to-results/FVE-novel-outputDirectory
cd FVE-novel-project-directory
javac -d bin src/*.java
java -cp bin FVENovel -1 $read1File -2 $read2File -o /path-to-results/FVE-novel-outputDirectory -fveres /path-to-FVE-res/FVE-outputDirectory -dbType gov -dbDir /path-to-referencedb-folder/GOV-viral-populations
```
To quickly run and test FastViromeExplorer-novel set the -topBins paprameter to 1. Then instead of running for top 100 bins, it will run for top 1 bin.
```bash
mkdir /path-to-results/FVE-novel-outputDirectory
cd FVE-novel-project-directory
javac -d bin src/*.java
java -cp bin FVENovel -1 $read1File -2 $read2File -o /path-to-results/FVE-novel-outputDirectory -fveres /path-to-results/FVE-outputDirectory -dbType gov -dbDir /path-to-referencedb-folder/GOV-viral-populations -topBins 1
```

# Output
The ouput folders and files will be generated in the `FVE-novel-outputDirectory` folder. The outputs will be stored in `final-results` folder inside `FVE-novel-outputDirectory` folder.

The output files are:
1. *final-scaffolds.fasta* : all scaffolds assembled from the reads, they are draft viral genomes
2. *bin-info.tsv* : For each scaffold, gives the bin information of that scaffold.
3. *ANI-info.tsv* : For each scaffold, gives the ANI(Average Nucleotide Identity) of this scaffold with the reference genome (and all other reference genomes present in the same bin) from where this scaffold got assembled.  