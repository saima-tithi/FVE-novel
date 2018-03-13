#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified October 15, 2015

Description:  Normalizes and error-corrects reads.
Designed for Jigsaw and uses hard-coded values.

Usage:
normandcorrectwrapper.sh in=reads.fq out=corrected.fq

Parameters:
in=<file>           Input reads filename.
out=<file>          Output reads filename.
correctfirst=f      Correct reads before normalization.
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

z="-Xmx14g"
z2="-Xms14g"
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
	freeRam 15000m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

normandcorrectwrapper() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_144_64bit
		module load pigz
	fi
	local CMD="java $EA $EOOM $z $z2 -cp $CP jgi.NormAndCorrectWrapper $@"
	echo $CMD >&2
	eval $CMD
}

normandcorrectwrapper "$@"
