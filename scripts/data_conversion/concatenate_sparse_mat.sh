#!/bin/bash

source $HOME/.bashrc
source $HOME/.bash_profile

MATDIR=$1
RESOL=$2
BIN=$3
OUTDIR=$4


LIBIDs=("ESR8" "ESR10" "H1Esc_1" "H1Esc_2" "ESR7" "H1Esc_HFF_1" "H1Esc_HFF_2" "HFF_GM_1" "HFF_GM_2" "IMR90_Hap1_rep1" "IMR90_Hap1_rep2" "GM12878_IMR90")

for LIBNAME in "${LIBIDs[@]}"; do

cat $MATDIR/${LIBNAME}/${RESOL}_matrices/*.sparse.matrix_${BIN} > $OUTDIR/sparse_matrices_${LIBNAME}_${BIN}.txt

done

