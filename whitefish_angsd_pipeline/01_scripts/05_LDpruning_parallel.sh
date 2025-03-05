#!/bin/bash

#To run on VALERIA:
#parallel -a ./02_info/regions.txt -j 10 srun -p ibis_small -c 1 --mem=30G --time=1-00:00 -J 05_LDpruning_{} -o 05_LDpruning_{}_%j.log /bin/sh 01_scripts/05_LDpruning_parallel.sh {} &

module load StdEnv/2020 gcc/9.3.0 python/3.8 graph-tool/2.45

REGIONS="$1"

01_scripts/ngsLD/scripts/prune_ngsLD.py --input 05_ngsLD/all_maf0.05pctind0.8_maxdepth10_"$REGIONS"_singletons.without_nan.ld --max_dist 250000 --min_weight 0.05 --field_dist 3 --field_weight 7 --weight_precision 4 --output 05_ngsLD/all_maf0.05pctind0.8_maxdepth10_"$REGIONS".singletons.pruned_maxdist_200kb_weight_0.05.ld

