#!/bin/bash
#SBATCH -J "het_by_indv"
#SBATCH -o log_%j
#SBATCH -c 1
#SBATCH -p small
#SBATCH --mail-type=ALL
#SBATCH --mail-user=YOURMAIL
#SBATCH --time=1-00:00
#SBATCH --mem=15G

#cat 00_remove_missing/bam2keep.filelist | parallel -j40 srun -c 1 --mem 10G -p small -o log_%j --time 1-00:00 ./01_scripts/14_het_by_ind_winsfs.sh {} &

###this script will work on each bamfiles and calculate individual saf and SFS for heterozygosity
#maybe edit
NB_CPU=1 #change accordingly in SLURM header

# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR
module load angsd
ulimit -S -n 2048

#prepare variables - avoid to modify
source 01_scripts/01_config.sh

BAM=$1
ID=$(grep $BAM 00_remove_missing/bam2keep.filelist | grep -o [A-Z][A-Z][A-Z][m-s]_[0-9][0-9][0-9]*-[0-9][0-9])

# Do SAF for one sample

echo $ID

 winsfs shuffle -v --output 14_heterozygosity/"$ID".saf.shuf \
    14_heterozygosity/"$ID".saf.idx

 winsfs -v 14_heterozygosity/"$ID".saf.shuf > 14_heterozygosity/"$ID".dsfs

 rm 14_heterozygosity/"$ID".saf.shuf

 #Remove header
 
 tail -n 1 14_heterozygosity/"$ID".dsfs > 14_heterozygosity/"$ID".dsfs.tmp
 mv 14_heterozygosity/"$ID".dsfs.tmp 14_heterozygosity/"$ID".dsfs

 ### Extract heterozygosity
 HET=$(awk '{ sum += $2; total += $0 } END { print sum / total }' 14_heterozygosity/"$ID".dsfs)
 echo -e $ID"\t"$HET >> 14_heterozygosity/summary.het
