#!/bin/bash
#SBATCH -J "WGSassign_safo"
#SBATCH -o log_%j
#SBATCH -c 10
#SBATCH -p medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=raphael.bouchard.3@ulaval.ca
#SBATCH --time=1-00:00
#SBATCH --mem=10G

for i in $(seq 1 10); 
do
nonbreeding_beagle=13_wgs_assign/non_essbp_loo_split/set_${i}_non_essbp_samples.singletons.pruned.clean.beagle.gz
breeding_af=13_wgs_assign/essbp_loo_split/set_${i}_cocl_180_ind_eesbp_230697_pruned_singleton_snp.pop_af.npy
outname=13_wgs_assign/non_essbp_loo_split/cocl_276_remaining_samples_assignment_set_${i}

WGSassign --beagle ${nonbreeding_beagle} --pop_af_file ${breeding_af} --get_pop_like --out ${outname} --threads 10
done
