#!/bin/bash
cd $SLURM_SUBMIT_DIR

#To run
#parallel -a ./02_info/regions.txt -j 4 srun -p ibis_medium -c 5 --mem=20G --time=1-00:00 -J 05_ngsLD_{} -o 05_ngsLD_per_chr_{}_%j.log /bin/sh 01_scripts/05_run_ngsLD_parallel.sh {} &

#module to load on Valeria
module load gsl StdEnv/2018.3 gcc/7.3.0 rstudio-server/1.2.1335


#Important ngsLD parameters:

#--probs: specification of whether the input is genotype probabilities (likelihoods or posteriors)?
#--n_ind INT: sample size (number of individuals).
#--n_sites INT: total number of sites.
#--max_kb_dist DOUBLE: maximum distance between SNPs (in Kb) to calculate LD. Set to 0(zero) to disable filter. [100]
#--max_snp_dist INT: maximum distance between SNPs (in number of SNPs) to calculate LD. Set to 0 (zero) to disable filter. [0]
#--rnd_sample= sample x% of comparison between snp
#--n_threads INT: number of threads to use. [1]
#--out FILE: output file name. [stdout]
#
REGION="$1"
NGSLD=01_scripts/ngsLD
NIND=$(wc -l 02_info/valeria_bam/bam.filelist)
NCPU=5

#Prepare pos file

awk '{print $1"\t"$2}' 04_ngsParalog/singletons_"$REGION".txt > 05_ngsLD/singletons_"$REGION".pos

#prepare variables - avoid to modify
source 01_scripts/01_config.sh
N_IND=$(wc -l 02_info/valeria_bam/bam.filelist | cut -d " " -f 1)
MIN_IND_FLOAT=$(echo "($N_IND * $PERCENT_IND)"| bc -l)
MIN_IND=${MIN_IND_FLOAT%.*}
MAX_DEPTH=$(echo "($N_IND * $MAX_DEPTH_FACTOR)" |bc -l)

NSITES=$(wc -l 05_ngsLD/singletons_"$REGION".pos)

$NGSLD/ngsLD --geno 04_ngsParalog/all_maf"$MIN_MAF"pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_"$REGION"_singletons.beagle.gz \
	--probs \
	--pos 05_ngsLD/singletons_"$REGION".pos \
	--n_ind "$NIND" \
	--n_sites "$NSITES" \
	--min_maf 0.05 \
	--max_kb_dist 1000 \
	--rnd_sample 0.5 \
	--n_threads "$NCPU" \
	--out 05_ngsLD/all_maf"$MIN_MAF"pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_"$REGION"_singletons.ld 
