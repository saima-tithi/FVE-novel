#!/bin/bash
#bbqc in=<infile> out=<outfile>

usage(){
echo "
Written by Brian Bushnell
Last modified June 4, 2015

*** DEPRECATED! Please use RQCFilter instead. ***

Description:  Performs quality-trimming; artifact, human, and phiX removal; adapter-trimming; error-correction and normalization.
Designed for Illumina fragment libraries only.
NOTE!  This program uses hard-coded paths and will only run on Genepool.

Usage:        bbqc.sh in=<input file> out=<output file>

Input parameters:
in=<file>           Input reads.
in2=<file>          Use this if 2nd read of pairs are in a different file.
ref=<file,file>     Comma-delimited list of additional reference files for filtering.
artifactdb=<file>   Override default Illumina artifacts file.
rnadb=<file>        Override default rna spikein file.
dnadb=<file>        Override default dna spikein file.
phixref=<file>      Override default phiX reference file.

Output parameters:
path=null           Set to the directory to use for all output files.
out=null            Read output file name.  If this is left blank, the input filename will be used with '.filtered' inserted before the extension.
scafstats=scaffoldStats.txt     Scaffold stats file name (how many reads matched which reference scaffold) .
kmerstats=kmerStats.txt         Kmer stats file name (duk-like output).
log=status.log                  Progress log file name.
filelist=file-list.txt          Progress log file name.
stats=filterStats.txt           Overall stats file name.
reproduceName=reproduce.sh      Name of shellscript to reproduce these results.

Processing parameters:
rna=f               Set to 't' for rna libraries, 'f' for dna libraries.
phix=t              Set to 't' to remove reads containing phiX kmers.
pjet=t              Set to 'f' to skip removal of reads containing PJET kmers.
threads=auto        (t) Set number of threads to use; default is number of logical processors.
filterk=27          Kmer length for finding contaminants.  Contaminants shorter than k will not be found..
trimk=23            Kmer length for linker/adapter trimming.
mink=11             Minimum kmer length for short kmers when trimming.
mapk=13             Kmer length for mapping to human.
forcetrimmod=5      (ftm) If positive, right-trim length to be equal to zero, modulo this number.
normalizek=27       Kmer length for normalization/error correction.
rcomp=t             Look for reverse-complements of kmers in addition to forward kmers.
maskmiddle=t        (mm) Treat the middle base of a kmer as a wildcard, to increase sensitivity in the presence of errors.
maxbadkmers=0       (mbk) Reads with more than this many contaminant kmers will be discarded.
filterhdist=1       Hamming distance used for filtering.
trimhdist=1         Hamming distance used for trimming.
trimhdist2=         Hamming distance used for trimming with short kmers.  If unset, trimhdist will be used.
mapindex=           Remove contaminants by mapping to the index at this location; default is an HG19 index.
mapref=             Remove contaminants by mapping to this fasta file.  Overrides mapindex flag.
human=t             Perform human contaminant removal.

Quality trimming parameters:
qtrim=rl            Trim read ends to remove bases with quality below minq.  Performed AFTER looking for kmers.
                    Values: t (trim both ends), f (neither end), r (right end only), l (left end only).
trimq=10            Trim quality threshold.
minlength=40        (ml) Reads shorter than this after trimming will be discarded.  Pairs will be discarded only if both are shorter.
minlengthfraction=0.6   (mlf) Reads shorter than this fraction of original length after trimming will be discarded.
minavgquality=8     (maq) Reads with average quality (after trimming) below this will be discarded.
maxns=1             Reads with more Ns than this (after trimming) will be discarded.
forcetrimmod=5      (ftm) If positive, trim length to be equal to zero modulo this number.

Normalization/Error Correction Paramters:
ecc=f               Correct errors.
aecc=f              Aggressive error correction (fixes more errors; less conservative).
cecc=f              Conservative error correction.
normalize=f         Normalize.
target=50           Target kmer depth after normalization.
min=-1              Throw away reads below this depth (default is 1/8th of target).
max=-1              Retain all reads between above target but below this depth (default is 11/10ths of target).
passes=2            Override number of passes for normalization.  By default it is 2 for normalization and 1 for error-correction only.

Java Parameters:
-Xmx                This will be passed to Java to set memory usage, overriding the program's automatic memory detection.
                    -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  The max is typically 85% of physical memory.
-eoom               This flag will cause the process to exit if an out-of-memory exception occurs.  Requires Java 8u92+.
-da                 Disable assertions.
"
}

pushd . > /dev/null
DIR="${BASH_SOURCE[0]}"
while [ -h "$DIR" ]; do
  cd "$(dirname "$DIR")"
  DIR="$(readlink "$(basename "$DIR")")"
done
cd "$(dirname "$DIR")"
DIR="$(pwd)/"
popd > /dev/null

#DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"
NATIVELIBDIR="$DIR""jni/"

z="-Xmx1g"
z2="-Xms1g"
EA="-ea"
EOOM=""
set=0

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

calcXmx () {
	source "$DIR""/calcmem.sh"
	parseXmx "$@"
	if [[ $set == 1 ]]; then
	return
	fi
	freeRam 9500m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

bbqc() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_144_64bit
		module load pigz
	export TZ="America/Los_Angeles" 
	fi
	local CMD="java -Djava.library.path=$NATIVELIBDIR $EA $z -cp $CP jgi.BBQC usejni $@"
	echo $CMD >&2
	eval $CMD
}

bbqc "$@"
