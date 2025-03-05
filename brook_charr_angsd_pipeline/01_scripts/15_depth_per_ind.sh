#!/bin/bash
#SBATCH -J "depth_per_ind"
#SBATCH -o log_%j
#SBATCH -c 1
#SBATCH -p small
#SBATCH --mail-type=ALL
#SBATCH --mail-user=YOURMAIL
#SBATCH --time=1-00:00
#SBATCH --mem=10G


#cat ../00_remove_missing/bam2keep.filelist | parallel -j20 srun -c 1 --mem 10G -p small -o ../log_%j --time 1-00:00 ../01_scripts/15_depth_per_ind.sh {}

BAM=$1
ID=$(grep $BAM ../00_remove_missing/bam2keep.filelist | grep -o [A-Z][A-Z][A-Z][m-s]_[0-9][0-9][0-9]*-[0-9][0-9]) 

mosdepth --by depth.bed $ID $BAM

