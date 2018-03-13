#!/bin/bash
/global/projectb/sandbox/gaag/bbtools/jgi-bbtools/taxserver.sh -da -Xmx45g port=3068 verbose accession=auto tree=auto table=auto domain=http://taxonomy.jgi-psf.org killcode=xxx 1>>taxlog01.o 2>&1 &
