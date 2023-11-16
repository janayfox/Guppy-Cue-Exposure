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

setwd("/scratch/janayfox/guppyWGBS/methylKit/dev/perc20")

## Load in data and reformat ##
ref.anno <- readTranscriptFeatures("/scratch/janayfox/guppyWGBS/annotation/guppy.sorted.bed", unique.prom = FALSE, up.flank = 1500, down.flank = 500)
DMR.diffmeth <- readRDS("./DMR_res/DMRdiffmeth5X.RDS")
DMS.diffmeth <- readRDS("./DMS_res/DMSdiffmeth5X.RDS")
CpG <- readRDS("./DMS_res/DMSmydiff5X.RDS")

#convert to GRanges
DMR.diffmeth.gr <- as(DMR.diffmeth, "GRanges") 
DMS.diffmeth.gr <- as(DMS.diffmeth, "GRanges") 
CpG.gr <- as(CpG, "GRanges") 

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

## Annotate ##
#annotate with gene parts 
DMR.anno <- annotateWithGeneParts(DMR.diffmeth.gr.rename, ref.anno)
DMS.anno <- annotateWithGeneParts(DMS.diffmeth.gr.rename, ref.anno)
CpG.anno <- annotateWithGeneParts(CpG.gr.rename, ref.anno)

#get percentage of DMRs that overlap with different features 
DMR.ann.stats.perc <- getTargetAnnotationStats(DMR.anno, percentage = TRUE, precedence = TRUE)
DMR.ann.stats.perc 
DMS.ann.stats.perc <- getTargetAnnotationStats(DMS.anno, percentage = TRUE, precedence = TRUE)
DMS.ann.stats.perc
CpG.ann.stats.perc <- getTargetAnnotationStats(CpG.anno, percentage = TRUE, precedence = TRUE)
CpG.ann.stats.perc

DMR.ann.stats.num <- getTargetAnnotationStats(DMR.anno, percentage = FALSE, precedence = TRUE)
DMR.ann.stats.num
DMS.ann.stats.num <- getTargetAnnotationStats(DMS.anno, percentage = FALSE, precedence = TRUE)
DMS.ann.stats.num
CpG.ann.stats.num <- getTargetAnnotationStats(CpG.anno, percentage = FALSE, precedence = TRUE)
CpG.ann.stats.num

#get nearest TSS for DMS and DMRs
DMR.tss <- getAssociationWithTSS(DMR.anno)
DMS.tss <- getAssociationWithTSS(DMS.anno)

## Save data ## 
saveRDS(DMR.anno, file = "./DMR_anno.RDS")
saveRDS(DMS.anno, file = "./DMS_anno.RDS")
saveRDS(CpG.anno, file = "./CpG_anno.RDS")

saveRDS(DMR.ann.stats.perc, file = "./DMR_annStats_perc.RDS")
saveRDS(DMS.ann.stats.perc, file = "./DMS_annStats_perc.RDS")
saveRDS(CpG.ann.stats.perc, file = "./CpG_annStats_perc.RDS")

saveRDS(DMR.ann.stats.num, file = "./DMR_annStats_num.RDS")
saveRDS(DMS.ann.stats.num, file = "./DMS_annStats_num.RDS")
saveRDS(CpG.ann.stats.num, file = "./CpG_annStats_num.RDS")

saveRDS(DMR.tss, file = "./DMR_TSS.RDS")
saveRDS(DMS.tss, file = "./DMS_TSS.RDS")
