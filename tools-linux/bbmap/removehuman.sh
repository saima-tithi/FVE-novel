#!/bin/bash

function usage(){
echo "
Written by Brian Bushnell
Last modified October 26, 2015
This script requires at least 16GB RAM.

Description:  Removes all reads that map to the human genome with at least 95% identity after quality trimming.
Removes approximately 98.6% of human 2x150bp reads, with zero false-positives to non-animals.
NOTE!  This program uses hard-coded paths and will only run on Nersc systems unless you change the path.

Usage:  removehuman.sh in=<input file> outu=<clean output file>

Input may be fasta or fastq, compressed or uncompressed.

Parameters:
threads=auto        (t) Set number of threads to use; default is number of logical processors.
overwrite=t         (ow) Set to false to force the program to abort rather than overwrite an existing file.
interleaved=auto    (int) If true, forces fastq input to be paired and interleaved.
trim=t              Trim read ends to remove bases with quality below minq.
                    Values: t (trim both ends), f (neither end), r (right end only), l (left end only).
untrim=t            Undo the trimming after mapping.
minq=4              Trim quality threshold.
ziplevel=2          (zl) Set to 1 (lowest) through 9 (max) to change compression level; lower compression is faster.
outm=<file>         File to output the reads that mapped to human.
path=               Set the path to an indexed human genome.

***** All BBMap parameters can be used; run bbmap.sh for more details. *****

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
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

z="-Xmx14000m"
EA="-ea"
EOOM=""
set=0

calcXmx () {
	source "$DIR""/calcmem.sh"
	parseXmx "$@"
}
calcXmx "$@"

function removehuman() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_144_64bit
		module load pigz
	fi
	local CMD="java -Djava.library.path=$NATIVELIBDIR $EA $z -cp $CP align2.BBMap minratio=0.9 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 path=/global/projectb/sandbox/gaag/bbtools/hg19 pigz unpigz zl=6 qtrim=r trimq=10 untrim idtag usemodulo printunmappedcount usejni ztd=2 kfilter=25 maxsites=1 k=14 $@"
	echo $CMD >&2
	eval $CMD
}

removehuman "$@"