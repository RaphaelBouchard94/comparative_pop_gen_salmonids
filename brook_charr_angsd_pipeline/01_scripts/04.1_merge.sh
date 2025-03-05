#!/bin/bash

# Combine each chromosome's output (.beagle and .mafs) from previous step (03_saf_maf_gl_all_parallel_LL.sh)

# VARIABLES
CHR_LIST="02_info/regions.txt"

#prepare variables - avoid to modify
source 01_scripts/01_config.sh
N_IND=$(wc -l 02_info/valeria_bam/bam.filelist | cut -d " " -f 1)
MIN_IND_FLOAT=$(echo "($N_IND * $PERCENT_IND)"| bc -l)
MIN_IND=${MIN_IND_FLOAT%.*}
MAX_DEPTH=$(echo "($N_IND * $MAX_DEPTH_FACTOR)" |bc -l)


# BEAGLE
# 1. Extract header for 1st chr : 
## what is the first chromosome ? 
FIRST_CHR=$(cat $CHR_LIST | head -n1) ## will be Chr01
## Extract header from beagle for first chr and initialize output file
zcat 04_ngsParalog/all_maf"$MIN_MAF"pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_"$FIRST_CHR"_singletons.beagle.gz | head -n1 > 04_ngsParalog/all_maf"$MIN_MAF"pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_ALL_CHR_singletons.beagle

# 2. Append beagle’s contents for all chromosomes
for i in $(cat 02_info/regions.txt)
do
  # extract the right beagle file for a given chr 
	# Extract all lines except first one and append to ALL_CHR.beagle
  zcat 04_ngsParalog/all_maf"$MIN_MAF"pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_${i}_singletons.beagle.gz  | tail -n+2  >> 04_ngsParalog/all_maf"$MIN_MAF"pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_ALL_CHR_singletons.beagle
done

gzip 04_ngsParalog/all_maf"$MIN_MAF"pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_ALL_CHR_singletons.beagle



