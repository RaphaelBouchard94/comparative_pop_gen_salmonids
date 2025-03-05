#!/bin/bash
#SBATCH -J "04_pca"
#SBATCH -o log_%j
#SBATCH -c 10
#SBATCH -p medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=raphael.bouchard.3@ulaval.ca
#SBATCH --time=1-00:00
#SBATCH --mem=100G

# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

#maybe edit
NB_CPU=10 #change accordingly in SLURM header

#prepare variables - avoid to modify
source 01_scripts/01_config.sh

#Module to load on Valeria
#module load StdEnv/2020 python/3.8 ; pip install --no-index numpy cython scipy

echo "analyse covariance matrix on all individuals"

pcangsd -b 06_pruned_singletons/all_maf0.05pctind0.8_maxdepth10_ALL_CHR.singletons.pruned.clean.beagle.gz \
	-e 4 \
	-t $NB_CPU \
	-o 07_pca/all_maf0.05pctind0.85_maxdepth10_ALL_CHR_singletons.pruned
