#!/bin/bash

#To run:
#parallel -a 02_info/sub_pop.txt -j 2 srun -p medium -c 2 --mem=20G --time=5-00:00 -J saf_{} -o 11_saf_maf_per_pop_{}_%j.log /bin/sh 01_scripts/11_calculate_saf_maf_per_site.sh {} &

#maybe edit
NB_CPU=5 #change accordingly in SLURM header
POP="$1"

# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR
module load angsd/0.931

#prepare variables - avoid to modify
source 01_scripts/01_config.sh


angsd -P $NB_CPU \
-dosaf 1 -doMajorMinor 4 -doMaf 1 -GL 2 \
-anc 02_info/genome.fasta -ref 02_info/genome.fasta \
-rf 02_info/regions.txt -sites 04_ngsParalog/singletons_ALL_CHR.sites \
-b 10_sites_bamfilelist/"$POP".cluster.sampled \
-out 11_saf_maf_by_pop/"$POP" 



