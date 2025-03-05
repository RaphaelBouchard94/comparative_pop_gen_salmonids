rm(list=ls())

setwd(dir = "~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/")

#---------------#
#--- Library ---#
#---------------#

library(tidyverse)
library(magrittr)
library(vcfR)
library(dartR)
library(adegenet)
library(radiator)
library(hierfstat)

##Write gene pop

writeGenPop <- function(gi, file.name, comment) {
  
  if (is.list(gi)) {
    # do all genind objects have the same number of loci?
    if (length(unique(sapply(gi, nLoc))) != 1) stop("Number of loci per individual genind object in a list is not equal for all.")
    gi.char <- gi[[1]]
    loc.names <- locNames(gi[[1]])
  } else {
    gi.char <- gi
    loc.names <- locNames(gi)
  }
  
  # Calculate the length of two alleles.
  lng <- as.character(na.omit(genind2df(gi.char)[, locNames(gi.char)[1]]))
  lng <- unique(nchar(lng))
  
  stopifnot(length(lng) == 1)
  
  cat(paste(comment, "\n"), file = file.name)
  cat(paste(paste(loc.names, collapse = ", "), "\n"), file = file.name, append = TRUE)
  
  if (is.list(gi)) {
    pop.names <- seq_len(length(gi))
  } else {
    pop.names <- popNames(gi)
  }
  
  for (i in pop.names) {
    cat("pop\n", file = file.name, append = TRUE)
    if (is.list(gi)) {
      intm <- gi[[i]]
      loc.names <- locNames(gi[[i]])
    } else {
      intm <- gi[pop(gi) == i, drop = FALSE]
    }
    ind.names <- indNames(intm)
    intm <- genind2df(intm, sep = "")
    intm[is.na(intm)] <- paste(rep("0", lng), collapse = "")
    out <- cbind(names = paste(ind.names, ",", sep = ""), intm[, loc.names])
    write.table(out, file = file.name, row.names = FALSE, col.names = FALSE, append = TRUE, quote = FALSE)
  }
  
  return(NULL)
}

alleleConvert<-function(combinedData){
  genotypes<-combinedData[,1:length(combinedData)-1]
  numericAlleles<-str_replace_all(genotypes,"NA/NA",paste(rep("0",as.numeric(combinedData)*2),collapse=""))
  numericAlleles<-str_replace_all(numericAlleles,c("A"="01","C"="02","G"="03","T"="04","-"="05"))
  numericAlleles<-str_replace_all(numericAlleles,",","")
}

#---------------#
#----- Data ----#
#---------------#


cocl_vcf <- read.vcfR(file = "../data/wgs_ne/filtered_cocl.vcf")
cocl_genind <- vcfR2genind(cocl_vcf)

pop <- row.names(cocl_genind@tab) %>% substr( 5, 7)

cocl_genind@pop <- as.factor(pop)

strata_file <- data.frame(INDIVIDUALS = row.names(cocl_genind@tab),
                          STRATA = pop)


# Get unique population levels
pop_levels <- levels(cocl_genind@pop)

# Loop through each population
for (pop in pop_levels) {
  # Extract individuals belonging to the population
  indexes <- grep(pop, cocl_gl@ind.names)
  subset_names <- cocl_gl@ind.names[indexes]
  
  # Define new population in gl object
  cocl_gl <- gl.define.pop(cocl_gl, ind.list = subset_names, new = pop)
  
  # Print progress
  cat("Subset defined for:", pop, "\n")
}

cocl_gl@pop

bea_gl <- cocl_gl[1:25,]

gl.LDNe(bea_gl, 
        outfile = "cocl_ldne.txt",
        outpath = "../data/wgs_ne/", 
        neest.path = "~/Documents/10_Programmes/ne_estimator/")





gl2genepop(
  cocl_gl,
  outfile = "cocl_genepop.gen",
  outpath = "data/wgs_ne/",
  pop_order = "alphabetic",
  output_format = "2_digits",
  verbose = NULL
)



