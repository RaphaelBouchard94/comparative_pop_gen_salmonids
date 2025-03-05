#!/bin/bash
#SBATCH -J "WGSassign_safo"
#SBATCH -o log_%j
#SBATCH -c 10
#SBATCH -p medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=raphael.bouchard.3@ulaval.ca
#SBATCH --time=1-00:00
#SBATCH --mem=10G


nonbreeding_beagle=13_wgs_assign/beagle/non_essbp_samples.singletons.pruned.clean.beagle.gz
breeding_af=13_wgs_assign/safo_254_essbp_ind_118861_pruned_singleton_snp.pop_af.npy
outname=13_wgs_assign/safo_120_remaining_samples_assignment


WGSassign --beagle ${nonbreeding_beagle} --pop_af_file ${breeding_af} --get_pop_like --out ${outname} --threads 10
