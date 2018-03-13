#!/bin/bash
/global/projectb/sandbox/gaag/bbtools/jgi-bbtools/taxserver.sh -da -Xmx8g port=3073 verbose tree=auto sketchonly index whitelist domain=https://ribo-sketch.jgi-psf.org killcode=xxx oldcode=xxx oldaddress=https://ribo-sketch.jgi-psf.org/kill/ /global/projectb/sandbox/gaag/bbtools/silva/bl_seq*.sketch dbname=Silva blacklist=silva k=31,0 1>>silvalog02.o 2>&1 &


