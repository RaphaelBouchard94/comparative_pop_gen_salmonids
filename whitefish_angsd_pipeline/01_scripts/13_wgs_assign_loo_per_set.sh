#!/bin/bash
#SBATCH -J "WGSassign_cocl"
#SBATCH -o log_%j
#SBATCH -c 10
#SBATCH -p medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=raphael.bouchard.3@ulaval.ca
#SBATCH --time=1-00:00
#SBATCH --mem=10G


for i in $(seq 1 10); 
do

breeding_beagle=13_wgs_assign/essbp_loo_split/set_${i}_essbp_samples.singletons.pruned.clean.beagle.gz 
breeding_IDs=13_wgs_assign/eesbp_samples_cluster.txt
outname=13_wgs_assign/essbp_loo_split/set_${i}_cocl_180_ind_eesbp_230697_pruned_singleton_snp

WGSassign --beagle ${breeding_beagle} --pop_af_IDs ${breeding_IDs} --get_reference_af --loo --out ${outname} --threads 10

done


