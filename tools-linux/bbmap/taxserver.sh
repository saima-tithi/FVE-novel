#!/bin/bash

usage(){
echo "
Written by Brian Bushnell and Shijie Yao
Last modified September 28, 2017

Description:   Starts a server that translates NCBI taxonomy.

Usage:  taxserver.sh tree=<taxtree file> table=<gitable file> port=<number>

Usage examples:
taxserver.sh tree=tree.taxtree.gz table=gitable.int1d.gz port=1234

On Genepool:
taxserver.sh tree=auto table=auto port=1234

For accession number support, add accession=<something>  E.g.:

External:
taxserver.sh -Xmx45g tree=tree.taxtree.gz table=gitable.int1d.gz accession=prot.accession2taxid.gz,nucl_wgs.accession2taxid.gz port=1234

On Genepool:
taxserver.sh tree=auto table=auto accession=auto port=1234

If all expected files are in some specific location, you can also do this:
taxserver.sh -Xmx45g tree=auto table=auto accession=auto port=1234 taxpath=/path/to/files

To kill remotely, launch with the flag kill=<password>, then access /kill/password

Parameters:

tree=auto           taxtree path.  Always necessary.
table=auto          gitable path.  Necessary for gi number support.
accession=null      Comma-delimited paths of accession files.
                    Necessary for accession support.
img=null            IMG dump file.
port=3068           Port number.
domain=             Domain to be displayed in the help message.
                    Default is taxonomy.jgi-psf.org.
dbname=             Set the name of the database in the help message.
killcode=           Set a password for remote killing.
oldcode=            Set the password of a prior instance.
oldaddress=         Attempt to kill a prior instance after initialization,
                    by sending the old code to this address.  For example,
                    taxonomy.jgi-psf.org/kill/
sketchthreads=16    Max concurrent threads for a sketch comparison is the
                    minimum of this and the threads argument.
sketchonly=f        Don't hash taxa names.
k=31                Kmer length, 1-32.  To maximize sensitivity and 
                    specificity, dual kmer lengths may be used:  k=31,24

Unrecognized parameters with no = symbol will be treated as sketch files.
Other sketch parameters such as index and k are also allowed.
Please consult bbmap/docs/guides/TaxonomyGuide.txt and BBSketchGuide.txt for more information.

Java Parameters:
-Xmx                This will be passed to Java to set memory usage, overriding the program's automatic memory detection.
                    -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  The max is typically 85% of physical memory.
-eoom               This flag will cause the process to exit if an out-of-memory exception occurs.  Requires Java 8u92+.
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

z="-Xmx45g"
z2="-Xms45g"
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

taxserver() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_64bit
		module load pigz
	fi
	local CMD="java $EA $EOOM $z $z2 -cp $CP tax.TaxServer $@"
	echo $CMD >&2
	eval $CMD
}

taxserver "$@"
