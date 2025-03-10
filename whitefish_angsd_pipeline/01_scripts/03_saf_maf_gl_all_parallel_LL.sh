#!/bin/bash 

# This script will parallelize instructions in 03_saf_maf_gl_all.sh by CHROMOSOME or region. 
# It requires a chromosome list : less 02_info/genome.fasta.fai | cut -f1 | grep ^Chr > 02_info/chr_list.txt (adapt grep pattern as required)
# For other unplaced contigs : less 02_info/genome.fasta.fai | cut -f1 | grep -v ^Chr > 02_info/contigs_list.txt (adapt grep pattern as required)

# To run with GNU parallel, use :
#parallel -a ./02_info/regions.txt -j 10 srun -p medium -c 2 --mem=50G --time=7-00:00 -J 03_saf_maf_gl_all_parallel_{} -o  03_saf_maf_gl_all_parallel_{}_%j.log /bin/sh 01_scripts/03_saf_maf_gl_all_parallel_LL.sh {} &

# VARIABLES
REGIONS="$1" # positional argument 1, refers to 1st arg provided when calling the script (03_saf_maf_gl_all_parallel_LL.sh HERE)
NB_CPU=2

# LOAD MODULES
module load angsd/0.931

ulimit -S -n 2048

#prepare variables - avoid to modify
source 01_scripts/01_config.sh
N_IND=$(wc -l 02_info/bam.filelist | cut -d " " -f 1)
MIN_IND_FLOAT=$(echo "($N_IND * $PERCENT_IND)"| bc -l)
MIN_IND=${MIN_IND_FLOAT%.*}
MAX_DEPTH=$(echo "($N_IND * $MAX_DEPTH_FACTOR)" |bc -l)

echo " Calculate the SAF, MAF and GL for all individuals listed in 02_info/bam_corrected.filelist"
echo "keep loci with at leat one read for n individuals = $MIN_IND, which is $PERCENT_IND % of total $N_IND individuals"
echo "filter on allele frequency = $MIN_MAF"

####Calculate the SAF, MAF and GL
angsd -P $NB_CPU -nQueueSize 50 \
-doMaf 1  -GL 2 -doGlf 2 -doMajorMinor 5 -doHWE 1 \
-anc 02_info/genome.fasta -remove_bads 1 -minMapQ 10 -minQ 20 -skipTriallelic 1 \
-doCounts 1 -dumpCounts 4 -minInd $MIN_IND -minMaf $MIN_MAF -setMaxDepth $MAX_DEPTH -setMinDepthInd $MIN_DEPTH \
-b 02_info/bam.filelist \
-r $REGIONS -out 03_saf_maf_gl_all/all_maf"$MIN_MAF"pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_"$REGIONS" #### ICI

#main features

#-P nb of threads -nQueueSize maximum waiting in memory (necesary to optimize CPU usage

#-doMaf 1 (allele frequencies) -dosaf (prior for SFS) -GL (Genotype likelihood 2 GATK method - export GL in beagle format -doGLF2)
#-doMajorMinor 1 use the most frequent allele as major
#-anc provide a ancestral sequence = reference in our case -fold 1 (car on utilise la ref comme ancestral
#-rf (file with the region written) work on a defined region : OPTIONAL
#-b (bamlist) input file
#-out output file
#main filters

#filter on bam files -remove_bads (remove files with flag above 255) -minMapQ minimum mapquality -minQ (minimum quality of reads?)

#filter on frequency -minInd (minimum number of individuals with at least one read at this locus) we set it to 50%

#filter on allele frequency -minMaf, set to 0.05

#extract SNP which passed the MIN_MAF and PERCENT_IND filters & their Major-minor alleles

#order the sites file by chromosome names

#makes a region file matching the sites files and with same order

#index sites file
echo "from the maf file, extract a list of SNP chr, positoin, major all, minor all"
gunzip 03_saf_maf_gl_all/all_maf"$MIN_MAF"pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_"$REGIONS".mafs.gz ## added _$REGIONS

INFILE=03_saf_maf_gl_all/all_maf"$MIN_MAF"pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_"$REGIONS".mafs ## added _$REGIONS
OUTFILE_sites=02_info/sites_all_maf"$MIN_MAF"pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_"$REGIONS" ## added _$REGIONS
OUTFILE_regions=02_info/regions_all_maf"$MIN_MAF"pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_"$REGIONS" ## added _$REGIONS

Rscript 01_scripts/Rscripts/make_sites_list_maxdepth_simple.R "$INFILE" "$OUTFILE_sites" "$OUTFILE_regions"

angsd sites index $OUTFILE_sites

#Eric propose a much shorter version using bash to cut the 4 columns of the mafs.gz. but then angsd is unable to index it

#I can't find the problem, so I came back to older solution with R which works

#gunzip -c 03_saf_maf_gl_all/all_maf"$MIN_MAF"pctind"$PERCENT_IND".mafs.gz | cut -f -4 > 02_info/sites_all_maf"$MIN_MAF"pctind"$PERCENT_IND"



# Serial command (do not use, testing purpose only)
# srun -p small -J 03_saf_maf_gl_all_parallel_ssa01q -o 03_saf_maf_gl_all_parallel_ssa01q_%j.log /bin/sh 01_scripts/03_saf_maf_gl_all_parallel.sh ssa01q &
