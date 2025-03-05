#!/bin/bash
#SBATCH -J "depth_per_ind"
#SBATCH -o log_%j
#SBATCH -c 1
#SBATCH -p small
#SBATCH --mail-type=ALL
#SBATCH --mail-user=YOURMAIL
#SBATCH --time=1-00:00
#SBATCH --mem=10G
#
#

for BAM in $(cat 00_remove_missing/bam2keep.filelist)
	do
		ID=$(grep $BAM 00_remove_missing/bam2keep.filelist | grep -o [A-Z][A-Z][A-Z][m-s]_[0-9][0-9][0-9]*-[0-9][0-9])
		zless 15_depth_per_ind/${ID}.per-base.bed.gz | cut -f 4 | sort -n | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }' >> 15_depth_per_ind/median_depth.per_ind.txt
	done


