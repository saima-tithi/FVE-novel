#!/bin/bash

# exit script if one command fails
set -o errexit
git clone https://github.com/saima-tithi/FVE-novel.git

CUR_DIR=$(pwd)

#Install bbmap
export PATH=$PATH:${CUR_DIR}/FVE-novel/tools-linux/bbmap
#Install Spades
export PATH=$PATH:${CUR_DIR}/FVE-novel/tools-linux/SPAdes-3.10.1-Linux/bin
#Install Salmon
export PATH=$PATH:${CUR_DIR}/FVE-novel/tools-linux/Salmon-latest_linux_x86_64/bin
#Install CD-HIT
export PATH=$PATH:${CUR_DIR}/FVE-novel/tools-linux/cd-hit-v4.6.8-2017-0621
#Install Bedtools
export PATH=$PATH:${CUR_DIR}/FVE-novel/tools-linux/bedtools2/bin
#Install Mummer
export PATH=$PATH:${CUR_DIR}/FVE-novel/tools-linux/MUMmer3.23