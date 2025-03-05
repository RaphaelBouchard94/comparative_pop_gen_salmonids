#!/bin/bash
#SBATCH -J "WGSassign_safo"
#SBATCH -o log_%j
#SBATCH -c 20
#SBATCH -p large
#SBATCH --mail-type=ALL
#SBATCH --mail-user=raphael.bouchard.3@ulaval.ca
#SBATCH --time=2-00:00
#SBATCH --mem=10G

breeding_beagle=06_pruned_singletons/all_maf0.05pctind0.8_maxdepth10_ALL_CHR.singletons.pruned.clean.beagle.gz
breeding_IDs=13_wgs_assign/get_ref_zscore/trout_wgs_assign_site_id.txt
pop_names=13_wgs_assign/get_ref_zscore/pop_names.txt
breeding_ad=18_allele_count/allele_count_ALL_CHR.majmin.counts.txt.gz
outname=13_wgs_assign/get_ref_zscore/safo_362_ind_118861_pruned_singleton

WGSassign --beagle ${breeding_beagle} --pop_af_IDs ${breeding_IDs} --pop_names ${pop_names} --ind_ad_file ${breeding_ad} --get_reference_z_score --out ${outname} --threads 20
