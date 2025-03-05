#!/bin/bash

# To run with GNU parallel, use :
#parallel -a ./02_info/regions.txt -j 20 srun -p batch_48h -c 2 --mem=25G --time=2-00:00 -J doPlink_{} -o  doPlink_{}_%j.log /bin/sh 01_scripts/16_call_vcf.sh {} &
#

# VARIABLES
REGIONS="$1" # positional argument 1, refers to 1st arg provided when calling the script (03_saf_maf_gl_all_parallel_LL.sh HERE)
NB_CPU=2

# LOAD MODULES
module load angsd

ulimit -S -n 2048

angsd -P $NB_CPU \
 -GL 1 -doMaf 2 -doMajorMinor 1 -doPlink 2 -doGeno 4 -doPost 1 -postCutoff 0.8 -SNP_pval 1e-6 \
 -ref 02_info/genome.fasta \
 -r $REGIONS -sites 04_ngsParalog/singletons_${REGIONS}.txt \
 -b 02_info/valeria_bam/bam.filelist \
 -out 16_vcf/singleton_${REGIONS}


