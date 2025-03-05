#!/bin/bash
#SBATCH -J "WGSassign_safo"
#SBATCH -o log_%j
#SBATCH -c 10
#SBATCH -p medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=raphael.bouchard.3@ulaval.ca
#SBATCH --time=1-00:00
#SBATCH --mem=10G

breeding_beagle=13_wgs_assign/beagle/essbp_samples.singletons.pruned.clean.beagle.gz
breeding_IDs=13_wgs_assign/eesbp_samples_cluster.txt
outname=13_wgs_assign/cocl_98_ind_eesbp_230697_pruned_singleton_snp


WGSassign --beagle ${breeding_beagle} --pop_af_IDs ${breeding_IDs} --get_reference_af --loo --out ${outname} --threads 10



