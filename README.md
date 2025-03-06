# Shared Dispersal Patterns but Contrasting Levels of Gene Flow in Two Anadromous Salmonids Along a Broad Subarctic Coastal Gradient

This is the page that compiles the scripts needed to replicate the analysis in the paper.

The lcWGS raw data for lake whitefish (PRJNA1037535, PRJNA1051576), and brook charr (PRJNA1160917, PRJNA1168880) are deposited online on NCBI Short Read Archive (SRA). 

The raw reads were processed using Eric Normandeau's pipeline available through his github: https://github.com/enormandeau/wgs_sample_preparation.git

`brook_charr_angsd_pipeline` and `whitefish_angsd_pipeline` are the scripts that process the bam files using angsd to get genotype likelihood and run population structure analysis on the servers such as pca, admixture, WGSassign.

`R_analysis` is the collection of Rscripts that I ran on my local computer to produce figures that appear on the paper. 

Please contact me if you have any problems: raphael.bouchard.3@ulaval.ca
