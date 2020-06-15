#!/bin/bash

source $HOME/.bashrc
source $HOME/.bash_profile

MATDIR=$1
BIN=$2
BINNAME=`expr $BIN \* 2`

LIBIDs=("ESR10" "ESR7" "ESR8" "GM12878_IMR90" "H1Esc_1" "H1Esc_2" "H1Esc_HFF_1" "H1Esc_HFF_2" "HFF_GM_1" "HFF_GM_2" "IMR90_Hap1_rep1" "IMR90_Hap1_rep2")
COUNT=0

## create output file
OUTFILE=$MATDIR/all_sparse_matrices_${BIN}.txt
touch $OUTFILE

## add to output file
for LIBNAME in "${LIBIDs[@]}"; do

awk -v count="$COUNT" '{print ($1+count), '\t', $2, '\t', $3}' $MATDIR/sparse_matrices_${LIBNAME}_${BIN}.txt | cat >> $OUTFILE

NUMCELLS=`awk 'END{print $1}' $MATDIR/sparse_matrices_${LIBNAME}_${BIN}.txt`
COUNT=`expr $COUNT + $NUMCELLS`

done

