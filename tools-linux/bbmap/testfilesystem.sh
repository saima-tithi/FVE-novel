#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified July 25, 2017

Logs filesystem performance by copying a file.

Usage:  testfilesystem.sh <in> <out> <log> <interval in seconds>

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

z="-Xmx50m"
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
}
calcXmx "$@"

function testfs() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_144_64bit
	fi
	local CMD="java $EA $EOOM $z -cp $CP jgi.TestFilesystem $@"
	echo $CMD >&2
	eval $CMD
}

testfs "$@"
