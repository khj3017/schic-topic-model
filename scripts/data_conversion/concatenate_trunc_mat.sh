#!/bin/bash

source $HOME/.bashrc
source $HOME/.bash_profile

MATDIR=$1
OUTDIR=$2

LIBIDs=("GM12878_IMR90") # "IMR90_Hap1_rep1" "IMR90_Hap1_rep2")

for LIBNAME in "${LIBIDs[@]}"; do

find $MATDIR/truncMat_${LIBNAME}_* -print0 | sort -z | xargs -0 cat > $OUTDIR/allMat_${LIBNAME}.txt
find $MATDIR/truncLPNames_${LIBNAME}_* -print0 | sort -z | xargs -0 cat > $OUTDIR/allLPNames_${LIBNAME}.txt

done

