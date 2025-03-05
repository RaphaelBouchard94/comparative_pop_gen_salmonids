#!/bin/bash
#SBATCH -J "05_ld_decay"
#SBATCH -o LD_decay
#SBATCH -c 1
#SBATCH -p medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=raphael.bouchard.3@ulaval.ca
#SBATCH --time=1-00:00
#SBATCH --mem=400G

#FILE variable contains a file which list the files to analyze
FILE="$1"

Rscript ../01_scripts/ngsLD/scripts/fit_LDdecay.R --ld_files $FILE --ld r2 --n_ind 378 --max_kb_dist 1000 --fit_level 100 --fit_boot 100  --plot_data -o all_chr_ld_decay.pdf

