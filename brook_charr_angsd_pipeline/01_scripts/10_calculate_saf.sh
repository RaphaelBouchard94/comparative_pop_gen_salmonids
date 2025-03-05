#!/bin/bash

#To run:
#parallel -a 02_info/pop.txt -j 10 srun -p ibis_medium -c 10 --mem=100G --time=5-00:00 -J saf_{} -o saf_{}_%j.log /bin/sh 01_scripts/10_calculate_saf.sh {} &

#maybe edit
NB_CPU=10 #change accordingly in SLURM header
POP="$1"

# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR
module load angsd/0.931

#prepare variables - avoid to modify
source 01_scripts/01_config.sh


angsd -P $NB_CPU \
-dosaf 1 -GL 2 \
-anc 02_info/genome.fasta \
-remove_bads 1 -minMapQ 10 -minQ 20 \
-rf 02_info/regions.txt -sites 05_ngsParalog/sites_singletons_all_chr \
-b 02_info/"$POP"bam.filelist \
-out 10_fst/"$POP" 



