#!/bin/bash 

# This script will parallelize genotype likelihood computation by chromosome or region. 
# It requires a chromosome list : less 02_info/genome.fasta.fai | cut -f1 | grep ^Chr > 02_info/chr_list.txt (adapt grep pattern as required)
# For other unplaced contigs : less 02_info/genome.fasta.fai | cut -f1 | grep -v ^Chr > 02_info/contigs_list.txt (adapt grep pattern as required)

# To run with GNU parallel, use :
#parallel -a ./02_info/regions.txt -j 10 srun -p small -c 2 --mem=25G --time=1-00:00 -J 03_saf_maf_gl_all_parallel_{} -o  03_saf_maf_gl_all_parallel_{}_%j.log /bin/sh 01_scripts/03_saf_maf_gl_all_parallel_manitou.sh {} &
#

# VARIABLES
REGIONS="$1" # positional argument 1, refers to 1st arg provided when calling the script (03_saf_maf_gl_all_parallel_LL.sh HERE)
NB_CPU=2

# LOAD MODULES
module load angsd

ulimit -S -n 2048


###############
###FEATURES####
###############
#-P: nb of threads 
#-nQueueSize maximum waiting in memory (necesary to optimize CPU usage)
#-doMaf 1: allele frequencies 
#-dosaf: (prior for SFS) 
#-GL: Genotype likelihood 2 GATK method - export GL in beagle format -doGLF2)
#-doMajorMinor: 1 use the most frequent allele as major
#-doHWE 1:Estimate the divination from HWE for each site
#-doCounts: Count the number A,C,G,T on all sites and all samples
#-dumpCounts: 
#### 1) Print overall depth in the .pos file 
#### 2) Prints the depth of each individual 
#### 3) Prints the depth for each of the four bases across all individuals 
#### 4) Prints the depth for each of the four bases for each indivial for each site

###############
####FILTERS####
###############
#-remove_bads : remove files with flag above 255
#-minMapQ: minimum map quality
#-minQ: minimum base quality score
#-skipTriallelic: calculates the pvalue of the site being triallelic by calculating the frequencies using the genotype likelihoods.
#The skipTriallelic is only used when you force a major and minor from external fastafiles.
#-minMAF: filter on minor allele frequency 
#-minInd: filter on minimum number of individuals with at least one read at this locus

#Prepare filter based on 01_config.sh file
source 01_scripts/01_config.sh
N_IND=$(wc -l 02_info/manitou_bam/bam.filelist | cut -d " " -f 1)
MIN_IND_FLOAT=$(echo "($N_IND * $PERCENT_IND)"| bc -l)
MIN_IND=${MIN_IND_FLOAT%.*}
MAX_DEPTH=$(echo "($N_IND * $MAX_DEPTH_FACTOR)" |bc -l)

echo " Calculate the MAF, GL, HWE, Counts for all individuals listed in 02_info/bam.filelist"
echo "keep loci with at leat one read for n individuals = $MIN_IND, which is $PERCENT_IND % of total $N_IND individuals"
echo "filter on allele frequency = $MIN_MAF"

###############
#####INPUTS####
###############
#-anc: provide an ancestral sequence = reference in our case
#-r: specify region or chromosome to work on
#-b: bam file list 
#-out: name of the output file


angsd -P $NB_CPU -nQueueSize 50 \
	-doMaf 1  -GL 2 -doGlf 2 -doMajorMinor 5 -doCounts 1 \
	-anc 02_info/genome.fasta \
	-remove_bads 1 -minMapQ 10 -minQ 20 -skipTriallelic 1 \
	-minInd $MIN_IND -minMaf $MIN_MAF -setMaxDepth $MAX_DEPTH -setMinDepthInd $MIN_DEPTH \
	-b 02_info/manitou_bam/bam.filelist \
	-r $REGIONS -out 03_saf_maf_gl_all/all_maf"$MIN_MAF"pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_"$REGIONS" #### ICI

#gunzip -c 03_saf_maf_gl_all/all_maf"$MIN_MAF"pctind"$PERCENT_IND".mafs.gz | cut -f -4 > 02_info/sites_all_maf"$MIN_MAF"pctind"$PERCENT_IND"

