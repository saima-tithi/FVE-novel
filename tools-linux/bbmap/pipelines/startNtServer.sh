#!/bin/bash
/global/projectb/sandbox/gaag/bbtools/jgi-bbtools/taxserver.sh -da -Xmx8g port=3071 verbose tree=auto sketchonly index domain=https://nt-sketch.jgi-psf.org killcode=xxx oldcode=xxx oldaddress=https://nt-sketch.jgi-psf.org/kill/ blacklist=nt /global/projectb/sandbox/gaag/bbtools/nt/current/taxa*.sketch dbname=nt k=31,24 1>>ntlog02.o 2>&1 &
