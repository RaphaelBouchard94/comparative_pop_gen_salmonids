#!/bin/bash
#SBATCH -J "09_pcadapt"
#SBATCH -o log_%j
#SBATCH -c 5
#SBATCH -p medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=raphael.bouchard.3@ulaval.ca
#SBATCH --time=1-00:00
#SBATCH --mem=100G

NCPU=5

pcangsd -b 04_ngsParalog/all_maf0.05pctind0.85_maxdepth10_ALL_CHR_singletons.beagle.gz -e 4 -t 6 --selection -o 09_fastpca/all_maf0.05pctind0.85_maxdepth10_ALL_CHR_singletons.fastpca
