#!/bin/bash
#SBATCH -J "het_by_indv"
#SBATCH -o log_%j
#SBATCH -c 1
#SBATCH -p batch_72h
#SBATCH --mail-type=ALL
#SBATCH --mail-user=YOURMAIL
#SBATCH --time=3-00:00
#SBATCH --mem=10G



for i in $(ls *.saf.idx); do realSFS $i > ${i%.*}.est.ml ; done
	
