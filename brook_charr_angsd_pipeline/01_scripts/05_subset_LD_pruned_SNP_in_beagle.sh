#!/bin/bash
#SBATCH -J "subset_ld_pruned_beagle"
#SBATCH -o log_%j
#SBATCH -c 1
#SBATCH -p ibis_small
#SBATCH --mail-type=ALL
#SBATCH --mail-user=raphael.bouchard.3@ulaval.ca
#SBATCH --time=1:00:00
#SBATCH --mem=10G

#This script subsets the LD pruned file list in a beagle file containing info on all chromosomes. 
#The LD_PRUNED_LIST file contains two columns chromosome and position separated by a tabulation.
#The script uses a python script beagle_extract_wanted_snps.py developped by Eric Normandeau.


source 01_scripts/01_config.sh
N_IND=$(wc -l 02_info/valeria_bam/bam.filelist | cut -d " " -f 1)
MIN_IND_FLOAT=$(echo "($N_IND * $PERCENT_IND)"| bc -l)
MIN_IND=${MIN_IND_FLOAT%.*}
MAX_DEPTH=$(echo "($N_IND * $MAX_DEPTH_FACTOR)" |bc -l)


BEAGLE_FILE=all_maf"$MIN_MAF"pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_ALL_CHR.beagle.gz
LD_PRUNED_LIST=all_maf"$MIN_MAF"pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_ALL_CHR.singletons.pruned.ld
LD_PRUNED_BEAGLE=all_maf"$MIN_MAF"pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_ALL_CHR.singletons.pruned.beagle.gz 

python 01_scripts/Scripts/beagle_extract_wanted_snps.py 03_saf_maf_gl_all/$BEAGLE_FILE 05_ngsLD/$LD_PRUNED_LIST 06_pruned_singletons/$LD_PRUNED_BEAGLE

