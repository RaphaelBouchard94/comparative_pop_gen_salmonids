#!/bin/bash
#SBATCH -J "19_ngsF"
#SBATCH -o log_%j
#SBATCH -c 2
#SBATCH -p medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=raphael.bouchard.3@ulaval.ca
#SBATCH --time=1-00:00
#SBATCH --mem=24G

#conda activate ngsf

./01_scripts/ngsF/ngsF --n_threads 2 --init_values r --min_epsilon 1e-9 --n_ind 362 --n_sites 118860 --glf 06_pruned_singletons/all_maf0.05pctind0.8_maxdepth10_ALL_CHR.singletons.pruned.clean.glf --out 19_ngsf/ngsf_safo
