#!/bin/bash

#this requires NGSadmix to be installed and its path export in the bashrc
#To run:
#parallel -a ./08_ngsadmix/k.txt -j 2 srun -p medium -c 20 --mem=200G --time=7-00:00 -J k_{} -o ngsadmixk_%j.log 01_scripts/08_ngsadmix.sh {} &

#maybe edit
NB_CPU=20 #change accordingly in SLURM header
K="$1"

# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

#prepare variables - avoid to modify
source 01_scripts/01_config.sh
source ~/.bashrc

for i in {1..50}; do NGSadmix -P $NB_CPU -likes 06_pruned_singletons/all_maf0.05pctind0.8_maxdepth10_ALL_CHR.singletons.pruned.clean.beagle.gz -K $K -outfiles 08_ngsadmix/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_pruned_singletons_"$K"_"$i" ; done
