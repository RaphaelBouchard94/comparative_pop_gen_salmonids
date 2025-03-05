#!/bin/bash
#SBATCH -J "het_by_indv"
#SBATCH -o log_%j
#SBATCH -c 1
#SBATCH -p bigmem_72h
#SBATCH --mail-type=ALL
#SBATCH --mail-user=YOURMAIL
#SBATCH --time=3-00:00
#SBATCH --mem=50G



for i in $(cat ../het_to_do.tsv); do realSFS ${i}.saf.idx > ${i}.saf.est.ml ; done
	
