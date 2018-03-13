#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified October 24, 2017

Description:  Compares query sketches to others, and prints their kmer identity.
The input can be sketches made by sketch.sh, or fasta/fastq files.
It's recommended to first sketch references with sketch.sh for large files,
or when taxonomic information is desired.

Please read bbmap/docs/guides/BBSketchGuide.txt for more information.

Usage:  comparesketch.sh in=<file,file,file...> ref=<file,file,file...>
Alternative:  comparesketch.sh in=<file,file,file...> file file file
Alternative:  comparesketch.sh in=<file,file,file...> *.sketch
Alternative:  comparesketch.sh alltoall *.sketch.gz

Standard parameters:
in=<file,file...>   Sketches or fasta files to compare.
out=stdout          Comparison output.  Can be set to a file instead.
ref=<file,file...>  List of sketches to compare against.  Files given without
                    a prefix (ref=) will be treated as references,
                    so you can use *.sketch without ref=.
                    You can also do ref=nt#.sketch to load all numbered files
                    fitting that pattern.
                    On NERSC, you can also specify these abbreviations:
                       nt:      Hard-coded path to nt sketches
                       refseq:  Hard-coded path to Refseq sketches
                       silva:   Hard-coded path to Silva sketches
                       img:     Hard-coded path to IMG sketches
                    Using an abbreviation automatically sets the blacklist, 
                    and k.
blacklist=<file>    Ignore keys in this sketch file.  Additionaly, there are
                    built-in blacklists that can be specified:
                       nt:      Blacklist for nt
                       refseq:  Blacklist for Refseq
                       silva:   Blacklist for Silva
                       img:     Blacklist for IMG
whitelist=f         Ignore keys that are not in the index.  Requires index=t.
threads=auto        Use this many threads for comparison.
index=auto          Index the sketches for much faster searching.
                    Requires more memory and adds startup time.
                    Recommended true for many query sketches, false for few.
amino=f             Use amino acid mode.
alltoall            (ata) Compare all refs to all.  Must be sketches.
compareself=f       In all-to-all mode, compare a sketch to itself.

Sketch-making parameters:
mode=single         Possible modes, for fasta input:
                       single: Generate one sketch per file.
                       sequence: Generate one sketch per sequence.
size=10000          Size of sketches to generate, if autosize=f.
                    For raw PacBio data, 100k is suggested.
maxfraction=0.01    (mgf) Max fraction of genomic kmers to use.
autosize=t          Use flexible sizing instead of fixed-length.
sizemult=1          Multiply the default size of sketches by this factor.
k=31                Kmer length, 1-32.  To maximize sensitivity and 
                    specificity, dual kmer lengths may be used:  k=31,24
                    Dual kmers are fastest if the shorter is a multiple 
                    of 4.  Query and reference k must match.
keyfraction=0.2     Only consider this upper fraction of keyspace.
samplerate=1        Set to a lower value to sample a fraction of input reads.
                    For raw reads (rather than an assembly), 1-3x coverage
                    gives best results, by reducing error kmers.  Somewhat
                    higher is better for high-error-rate data like PacBio.
minkeycount=1       Ignore kmers that occur fewer times than this.  Values
                    over 1 can be used with raw reads to avoid error kmers.
sketchheapfactor=4  If minkeycount>1, temporarily track this many kmers until
                    counts are known and low-count kmers are discarded.
minprob=0.0001      Ignore kmers below this probability of correctness.
minqual=0           Ignore kmers spanning bases below this quality.
entropy=0.66        Ignore sequence with entropy below this value.
merge=f             Merge paired reads prior to sketching.

Taxonomy-related parameters:
tree=<file>         Specify a TaxTree file.  On Genepool, use tree=auto.
                    Only necessary for use with printtaxa and level.
                    Assumes comparisons are done against reference sketches
                    with known taxonomy information.
level=2             Only report the best record per taxa at this level.
                    Either level names or numbers may be used.
                       -1: disabled
                        1: subspecies
                        2: species
                        3: genus
                       ...etc
taxfilter=          List of NCBI taxIDs to filter from output.
taxfilterlevel=0    Filter everything sharing an ancestor at this level.
taxfilterinclude=t  Include only filtered records. 'f' will discard them.

Output columns:
printall=f          Enable all output columns.
printani=t          (ani) Print average nucleotide identity estimate.
completeness=t      Genome completeness estimate.
score=f             Score (used for sorting the output).
printmatches=t      Number of kmer matches to reference.
printlength=f       Number of kmers compared.
printtaxid=t        NCBI taxID.
printimg=f          IMG identifier (only for IMG data).
printgbases=f       Number of genomic bases.
printgkmers=f       Number of genomic kmers.
printgsize=t        Estimated number of unique genomic kmers.
printgseqs=t        Number of sequences (scaffolds/reads).
printtaxname=t      Name associated with this taxID.
printname0=t        (pn0) Original seqeuence name.
printfname=t        Query filename.
printtaxa=f         Full taxonomy of each record.
printcontam=t       Print contamination estimate, and factor contaminant kmers
                    into calculations.  Kmers are considered contaminant if
                    present in some ref sketch but not the current one.
printunique=t       Number of matches unique to this reference.
printunique2=f      Number of matches unique to this reference's taxa.
printunique3=f      Number of query kmers unique to this reference's taxa,
                    regardless of whether they are in this reference sketch.
printnohit=t        Number of kmers that don't hit anything.
printucontam=f      Contam hits that hit exactly one reference sketch.
printcontam2=f      Print contamination estimate using only kmer hits
                    to unrelated taxa.
contamlevel=species Taxonomic level to use for contam2/unique2/unique3.
NOTE: unique2/unique3/contam2 require an index.

printdepth=f        (depth) Print average depth of sketch kmers; intended
                    for shotgun read input.
printdepth2=f       (depth2) Print depth compensating for genomic repeats.
                    Requires reference sketches to be generated with depth.
actualdepth=t       If this is false, the raw average count is printed.
                    If true, the raw average (observed depth) is converted 
                    to estimated actual depth (including uncovered areas).
printvolume=f       (volume) Product of average depth and matches.

Sorting:
sortbyscore=t       Default sort order is by score.
sortbydepth=f       Include depth as a factor in sort order.
sortbydepth2=f      Include depth2 as a factor in sort order.
sortbyvolume=f      Include volume as a factor in sort order.

Other output parameters:
minhits=3           (hits) Only report records with at least this many hits.
minani=0            (ani) Only report records with at least this ANI (0-1).
minwkid=0.0001      (wkid) Only report records with at least this WKID (0-1).
records=20          Report at most this many best-matching records.
color=family        Color records at the family level.  color=f will disable.
                    Colors work in most terminals but may cause odd characters
                    to appear in text editors.  So, color defaults to f if 
                    writing to a file.  Requires the taxtree to be loaded.
format=2            Available formats:
                      1: Original format; deprecated.
                      2: Customizable with different columns.
                      3: 3-column format for alltoall: query, ref, ANI.
intersect=f         Print sketch intersections.  delta=f is suggested.

Metadata flags (optional, for the query sketch header):
taxid=-1            Set the NCBI taxid.
imgid=-1            Set the IMG id.
spid=-1             Set the sequencing project id (JGI-specific).
name=               Set the name (taxname).
name0=              Set name0 (normally the first sequence header).
fname=              Set fname (normally the file name).
meta_=              Set an arbitrary metadata field.
                    For example, meta_Month=March.

Java Parameters:
-Xmx                This will be passed to Java to set memory usage, overriding the program's automatic memory detection.
                    -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  The max is typically 85% of physical memory.
-eoom               This flag will cause the process to exit if an out-of-memory exception occurs.  Requires Java 8u92+.
-da                 Disable assertions.

For more detailed information, please read /bbmap/docs/guides/BBSketchGuide.txt.
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

z="-Xmx4g"
z2="-Xms4g"
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
	freeRam 3200m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

comparesketch() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_64bit
		module load pigz
	fi
	local CMD="java $EA $EOOM $z $z2 -cp $CP sketch.CompareSketch $@"
#	echo $CMD >&2
	eval $CMD
}

comparesketch "$@"
