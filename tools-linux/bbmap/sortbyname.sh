#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified June 14, 2017

Description:  Sorts reads by name.
Writes temp files if memory is exceeded.

Usage:   sortbyname.sh in=<file> out=<file>

Input may be fasta, fastq, or sam, compressed or uncompressed.

Parameters:

in=<file>           Input file.
out=<file>          Output file.
delete=t            Delete temporary files.
allowtemp=t         Allow writing temporary files.
name=t              Sort reads by name.
length=f            Sort reads by length.
quality=f           Sort reads by quality.
position=f          Sort reads by position (for mapped reads).
taxa=f              Sort reads by taxonomy (for NCBI naming convention).
sequence=f          Sort reads by sequence, alphabetically.
flowcell=f          Sort reads by flowcell coordinates.
shuffle=f           Shuffle reads randomly (untested).
list=<file>         Sort reads according to this list of names.
ascending=t         Sort ascending.

Taxonomy-sorting parameters:
tree=               Specify a taxtree file.  On Genepool, use 'auto'.
gi=                 Specify a gitable file.  On Genepool, use 'auto'.
accession=          Specify one or more comma-delimited NCBI accession to
                    taxid files.  On Genepool, use 'auto'.

Note: name, length, and quality are mutually exclusive.
Sorting by quality actually sorts by average expected error rate,
so ascending will place the highest-quality reads first.

Java Parameters:
-Xmx                This will be passed to Java to set memory usage, overriding
                    the program's automatic memory detection.  -Xmx20g will 
                    specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
                    The max is typically 85% of physical memory.
-eoom               This flag will cause the process to exit if an 
                    out-of-memory exception occurs.  Requires Java 8u92+.
-da                 Disable assertions.

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

z="-Xmx2g"
z2="-Xms2g"
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
	freeRam 2000m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

sortbyname() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_64bit
		module load pigz
	fi
	local CMD="java $EA $EOOM $z $z2 -cp $CP sort.SortByName $@"
	echo $CMD >&2
	eval $CMD
}

sortbyname "$@"
