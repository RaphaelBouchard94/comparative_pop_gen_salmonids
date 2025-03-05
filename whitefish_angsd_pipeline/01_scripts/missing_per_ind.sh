#!/bin/bash
#SBATCH -J "missing_data_per_ind"
#SBATCH -o missing_per_ind.log
#SBATCH -c 2
#SBATCH -p small
#SBATCH --mail-type=ALL
#SBATCH --mail-user=raphael.bouchard.3@ulaval.ca
#SBATCH --time=1-00:00
#SBATCH --mem=10G

cd $SLURM_SUBMIT_DIR

zcat 03_saf_maf_gl_all/all_maf0.05pctind0.8_maxdepth10_ALL_CHR.beagle.gz | java -jar 01_scripts/gprobssamplemissing.jar 0.5 > missing_per_ind.txt
