#!/bin/bash
#BATCH -J "het_by_indv"
#SBATCH -o log_%j
#SBATCH -c 10
#SBATCH -p bigmem_96h
#SBATCH --mail-type=ALL
#SBATCH --mail-user=YOURMAIL
#SBATCH --time=4-00:00
#SBATCH --mem=100G


#maybe edit
NB_CPU=10 #change accordingly in SLURM header

module load angsd

angsd -P $NB_CPU \
 -dosaf 1 -anc 02_info/genome.fasta -gl 1\
 -rf 02_info/regions.txt -sites 04_ngsParalog/singletons_ALL_CHR.sites \
 -b 02_info/valeria_bam/bam.filelist \
 -out 16_vcf/all_singletons_cocl
