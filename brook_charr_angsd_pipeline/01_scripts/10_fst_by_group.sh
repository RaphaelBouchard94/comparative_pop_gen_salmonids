#!/bin/bash
#SBATCH -J "10_FST_by_group" 
#SBATCH -o fst_by_group_%j
#SBATCH -c 5
#SBATCH -p medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=raphael.bouchard.3@ulaval.ca
#SBATCH --time=1-00:00
#SBATCH --mem=100G


#Input variables:
source 01_scripts/01_config.sh
NSITES=1000000
NB_CPU=10
POP_FILE=02_info/new_pop.txt 
##Calculate number of pop
num_pops=$(wc -l "$POP_FILE" | cut -d " " -f 1)


# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR
module load angsd
ulimit -S -n 2048


# Estimate pairwise FST for all populations listed

for i in $(seq $num_pops)
do
	pop1=$(cat "$POP_FILE" | head -"$i" | tail -1)
	for j in $(seq $[ $i + 1 ] $num_pops)
	do
		pop2=$(cat "$POP_FILE" | head -"$j" | tail -1)
		echo "FST between $pop1 and $pop2"
		echo "$pop1"
		echo "$pop2"
		
		echo "calcualte the 2dsfs priors"
		realSFS  10_fst/"$pop1".saf.idx 10_fst/"$pop2".saf.idx \
		-P $NB_CPU -maxIter 30 -nSites $NSITES > 10_fst/"$pop1"_"$pop2"."$NSITES"

		file=10_fst/"$pop1"_"$pop2"."$NSITES"

		Rscript 01_scripts/Rscripts/sum_sites_2dsfs.r "$file"
		
		echo " prepare the fst for easy window analysis etc"
		realSFS fst index 10_fst/"$pop1".saf.idx 10_fst/"$pop2".saf.idx \
		-sfs 10_fst/"$pop1"_"$pop2"."$NSITES".2dsfs \
		-P $NB_CPU -fstout 10_fst/"$pop1"_"$pop2"

		echo "print SFS priori for each position"
		realSFS fst print 10_fst/"$pop1"_"$pop2".fst.idx \ 
		-P $NB_CPU > 10_fst/"$pop1"_"$pop2".bypos.sfs
		
		echo "get the global estimate of FST throughout the genome"
		realSFS fst stats 10_fst/"$pop1"_"$pop2".fst.idx \
		-P $NB_CPU > 10_fst/$GROUP/"$pop1"_"$pop2".fst
		
		echo "calculate FST by slidingwindow, window size=$WINDOW and step=$WINDOW_STEP, as given in 01_config.sh"
		realSFS  fst stats2 10_fst/"$pop1"_"$pop2".fst.idx \
		-win $WINDOW -step $WINDOW_STEP -P $NB_CPU > 10_fst/"$pop1"_"$pop2".slidingwindow
	done
done
