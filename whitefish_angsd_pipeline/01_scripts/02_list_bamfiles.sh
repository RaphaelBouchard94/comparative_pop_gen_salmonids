###this script will make a list of all bamfiles and several list by population

#variables
source 01_scripts/01_config.sh

# Make a list of all bamfiles
ls /project/lbernatchez/01_projects/Projects/Raphael_Bouchard/2022-05-11_wgs_sample_preparation_coar/09_no_overlap/*.bam >> 02_info/bam.filelist

# Make a list of bamfiles by pop, as given by pop.txt
POP_FILE1=02_info/pop.txt #pop list n° 1 

cat $POP_FILE1 | while read i
do
echo $i
ls /project/lbernatchez/01_projects/Projects/Raphael_Bouchard/2022-05-12_wgs_sample_preparation_cocl/09_no_overlap/*$i*.bam >> 02_info/"$i"bam.filelist
done

##possible to re-do this for another list of pop, beware it is not the same names, or it will overwrite
#POP_FILE2=02_info/group.txt #pop list n° 2

#cat $POP_FILE2 | while read i
#do
#echo $i
#ls $BAM_PATH/*$i*.sorted.bam > 02_info/"$i"bam.filelist
#done
