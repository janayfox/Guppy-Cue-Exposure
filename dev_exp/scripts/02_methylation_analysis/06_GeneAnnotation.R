##############################################################
### Goal: Annotate CpGs, DMSs, and DMRs with genes parts
### Author: Janay Fox
### R script
#############################################################

## Set up ##
#install packages 
#install.packages("S4Vectors")
#install.packages("IRanges")
#install.packages("GenomicRanges")
#install.packages("methylKit")
#install.packages("genomation")

#load packages
library("S4Vectors", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")
library("IRanges", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")
library("GenomicRanges", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")
library("methylKit", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")
library("genomation", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")

setwd("/scratch/janayfox/guppyWGBS/")
load(file=".RData")

## Annotate DMRs ##
#load in annotation
gene.obj <- readTranscriptFeatures("/scratch/janayfox/guppyWGBS/annotation/guppy.sorted.bed", unique.prom = FALSE, up.flank = 1500, down.flank = 500)

#convert diffMeth to GRanges
DMR.diffmeth.gr <- as(DMR.diffMeth, "GRanges") #DMRs
DMS.diffmeth.gr <- as(DMS.diffMeth, "GRanges") #DMSs
CpG.gr <- as(DMS.myDiff, "GRanges") #CpGs

#change chromosome names to match 
#create function that renames chromosomes
renameChr <- function(gr.obj) {
  gr.obj.rename <- renameSeqlevels(gr.obj, c(NC_024331.1="LG1", NC_024332.1="LG2",
                                     NC_024333.1="LG3", NC_024334.1="LG4",
                                     NC_024335.1="LG5",NC_024336.1="LG6",
                                     NC_024337.1="LG7", NC_024338.1="LG8",
                                     NC_024339.1="LG9", NC_024340.1="LG10",
                                     NC_024341.1="LG11", NC_024343.1="LG13",
                                     NC_024344.1="LG14", NC_024345.1="LG15",
                                     NC_024346.1="LG16", NC_024347.1="LG17",
                                     NC_024348.1="LG18", NC_024349.1="LG19",
                                     NC_024350.1="LG20", NC_024351.1="LG21",
                                     NC_024352.1="LG22", NC_024353.1="LG23"))
  return(gr.obj.rename)
}

#rename chromosomes in each GRanges
DMR.diffmeth.gr.rename <- renameChr(DMR.diffmeth.gr)
DMS.diffmeth.gr.rename <- renameChr(DMS.diffmeth.gr)
CpG.gr.rename <- renameChr(CpG.gr)

#annotate with gene parts 
#DMR.diffMeth.an <- annotateWithGeneParts(DMR.diffmeth.gr.rename, gene.obj)
DMS.diffMeth.an <- annotateWithGeneParts(DMS.diffmeth.gr.rename, gene.obj)
CpG.an <- annotateWithGeneParts(CpG.gr.rename, gene.obj)

#get percentage of DMRs that overlap with different features 
DMR.ann.stats <- getTargetAnnotationStats(DMR.diffMeth.an, percentage = TRUE, precedence = TRUE)
DMS.ann.stats <- getTargetAnnotationStats(DMS.diffMeth.an, percentage = TRUE, precedence = TRUE)
CpG.ann.stats <- getTargetAnnotationStats(CpG.an, percentage = TRUE, precedence = TRUE)

#get nearest TSS
DMR.tss <- getAssociationWithTSS(DMR.diffMeth.an)
DMS.tss <- getAssociationWithTSS(DMS.diffMeth.an)
#CpG.tss <- getAssociationWithTSS(CpG.an)
#not sure that I need to do this for CpGs

## Save workspace image for later loading ##
save.image(file = ".RData")
save.image(file = "./backupRData/06_GeneAnnotation-backup.RData")
#maybe add in save rds??? 

q(save="yes")



