#!/bin/bash
#SBATCH -J "09_pcadapt"
#SBATCH -o log_%j
#SBATCH -c 20
#SBATCH -p medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=raphael.bouchard.3@ulaval.ca
#SBATCH --time=1-00:00
#SBATCH --mem=200G

NCPU=20

module load pcangsd

pcangsd -b 04_ngsParalog/all_maf0.05pctind0.8_maxdepth10_ALL_CHR_singletons.beagle.gz -e 4 -t 20 --pcadapt -o 09_pcadapt/all_maf0.05pctind0.8_maxdepth10_ALL_CHR_singletons
