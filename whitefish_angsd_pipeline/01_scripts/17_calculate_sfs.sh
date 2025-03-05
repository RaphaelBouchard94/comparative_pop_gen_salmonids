#!/bin/bash
#SBATCH -J "het_by_indv"
#SBATCH -o log_%j
#SBATCH -c 20
#SBATCH -p batch_72h
#SBATCH --mail-type=ALL
#SBATCH --mail-user=YOURMAIL
#SBATCH --time=3-00:00
#SBATCH --mem=100G

#cat 02_info/all_pop.txt | parallel -j 6 srun -c 20 --mem 100G -p large -o log_%j --time 3-00:00 ./01_scripts/17_calculate_sfs.sh {} &


POP=$1

module load angsd

realSFS 11_saf_maf_by_pop/${POP}.saf.idx -P 20 -fold 1 > 17_genetic_diversity/${POP}.sfs
