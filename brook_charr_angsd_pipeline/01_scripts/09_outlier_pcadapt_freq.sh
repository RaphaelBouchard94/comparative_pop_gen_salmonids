#!/bin/bash

#To run
#parallel -a 02_info/pop.txt -j 10 srun -p ibis_medium -c 5 --mem=100G --time=3-00:00 -J maf_{}_outliers -o outliers_pcadapt_{}_%j.log /bin/sh 01_scripts/09_outlier_pcadapt_freq.sh {} &


POP="$1"
NCPU=5

module load angsd

source 01_scripts/01_config.sh

angsd -P $NCPU \
-GL 2 -doMaf 1 -doMajorMinor 1 \
-minInd $MIN_IND -minMaf $MIN_MAF \
-b 02_info/"$POP"bam.filelist \
-anc 02_info/genome.fasta \
-rf 02_info/regions.txt -sites 09_pcadapt/outliers_pcadapt.tsv \
-out 09_pcadapt/"$POP"outliers

