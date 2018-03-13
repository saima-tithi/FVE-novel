#!/bin/bash
/global/projectb/sandbox/gaag/bbtools/jgi-bbtools/taxserver.sh -da -Xmx12g port=3072 verbose tree=auto sketchonly index domain=https://refseq-sketch.jgi-psf.org killcode=xxx oldcode=xxx oldaddress=https://refseq-sketch.jgi-psf.org/kill/ /global/projectb/sandbox/gaag/bbtools/refseq/current/taxa*.sketch dbname=RefSeq blacklist=refseq k=31,24 1>>refseqlog02.o 2>&1 &

