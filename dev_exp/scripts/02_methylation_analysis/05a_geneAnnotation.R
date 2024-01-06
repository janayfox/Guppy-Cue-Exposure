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

DMR.diffmeth.all <- readRDS("./DMR_res/DMRdiffmeth_all_5X.RDS")
DMR.diffmeth.fem <- readRDS("./DMR_res/DMRdiffmeth_fem_5X.RDS")
DMR.diffmeth.mal <- readRDS("./DMR_res/DMRdiffmeth_mal_5X.RDS")

DMS.diffmeth.all <- readRDS("./DMS_res/DMSdiffmeth_all_5X.RDS")
DMS.diffmeth.fem <- readRDS("./DMS_res/DMSdiffmeth_fem_5X.RDS")
DMS.diffmeth.mal <- readRDS("./DMS_res/DMSdiffmeth_mal_5X.RDS")

CpG.all <- readRDS("./DMS_res/DMSmydiff_all_5X.RDS")
CpG.fem <- readRDS("./DMS_res/DMSmydiff_fem_5X.RDS")
CpG.mal <- readRDS("./DMS_res/DMSmydiff_mal_5X.RDS")

regions.all <- readRDS("./DMR_res/DMRmydiff_all_5X.RDS")
regions.fem <- readRDS("./DMR_res/DMRmydiff_fem_5X.RDS")
regions.mal <- readRDS("./DMR_res/DMRmydiff_mal_5X.RDS")

#convert to GRanges
DMR.diffmeth.all.gr <- as(DMR.diffmeth.all, "GRanges") 
DMR.diffmeth.fem.gr <- as(DMR.diffmeth.fem, "GRanges") 
DMR.diffmeth.mal.gr <- as(DMR.diffmeth.mal, "GRanges") 

DMS.diffmeth.all.gr <- as(DMS.diffmeth.all, "GRanges") 
DMS.diffmeth.fem.gr <- as(DMS.diffmeth.fem, "GRanges") 
DMS.diffmeth.mal.gr <- as(DMS.diffmeth.mal, "GRanges") 

CpG.all.gr <- as(CpG.all, "GRanges") 
CpG.fem.gr <- as(CpG.fem, "GRanges") 
CpG.mal.gr <- as(CpG.mal, "GRanges") 

regions.all.gr <- as(regions.all, "GRanges") 
regions.fem.gr <- as(regions.fem, "GRanges") 
regions.mal.gr <- as(regions.mal, "GRanges") 

#change chromosome names to match 
#create function that renames chromosomes
renameChr <- function(gr.obj) {
  gr.obj.rename <- renameSeqlevels(gr.obj, c(NC_024331.1="LG1", NC_024332.1="LG2",
                                     NC_024333.1="LG3", NC_024334.1="LG4",
                                     NC_024335.1="LG5",NC_024336.1="LG6",
                                     NC_024337.1="LG7", NC_024338.1="LG8",
                                     NC_024339.1="LG9", NC_024340.1="LG10",
                                     NC_024341.1="LG11", NC_024342.1="LG12",
                                     NC_024343.1="LG13",
                                     NC_024344.1="LG14", NC_024345.1="LG15",
                                     NC_024346.1="LG16", NC_024347.1="LG17",
                                     NC_024348.1="LG18", NC_024349.1="LG19",
                                     NC_024350.1="LG20", NC_024351.1="LG21",
                                     NC_024352.1="LG22", NC_024353.1="LG23"))
  return(gr.obj.rename)
}

#rename chromosomes in each GRanges
DMR.diffmeth.all.gr.rename <- renameChr(DMR.diffmeth.all.gr)
DMR.diffmeth.fem.gr.rename <- renameChr(DMR.diffmeth.fem.gr)
DMR.diffmeth.mal.gr.rename <- renameChr(DMR.diffmeth.mal.gr)

DMS.diffmeth.all.gr.rename <- renameChr(DMS.diffmeth.all.gr)
DMS.diffmeth.fem.gr.rename <- renameChr(DMS.diffmeth.fem.gr)
DMS.diffmeth.mal.gr.rename <- renameChr(DMS.diffmeth.mal.gr)

CpG.all.gr.rename <- renameChr(CpG.all.gr)
CpG.fem.gr.rename <- renameChr(CpG.fem.gr)
CpG.mal.gr.rename <- renameChr(CpG.mal.gr)

regions.all.gr.rename <- renameChr(CpG.all.gr)
regions.fem.gr.rename <- renameChr(CpG.fem.gr)
regions.mal.gr.rename <- renameChr(CpG.mal.gr)

#get hyper and hypo methylated 
DMR.diffmeth.all.gr.rename.hyper <- subset(DMR.diffmeth.all.gr, meth.diff > 0)
DMR.diffmeth.fem.gr.rename.hyper <- subset(DMR.diffmeth.fem.gr, meth.diff > 0)
DMR.diffmeth.mal.gr.rename.hyper <- subset(DMR.diffmeth.mal.gr, meth.diff > 0)

DMS.diffmeth.all.gr.rename.hyper <- subset(DMS.diffmeth.all.gr, meth.diff > 0)
DMS.diffmeth.fem.gr.rename.hyper <- subset(DMS.diffmeth.fem.gr, meth.diff > 0)
DMS.diffmeth.mal.gr.rename.hyper <- subset(DMS.diffmeth.mal.gr, meth.diff > 0)

DMR.diffmeth.all.gr.rename.hypo <- subset(DMR.diffmeth.all.gr, meth.diff < 0)
DMR.diffmeth.fem.gr.rename.hypo <- subset(DMR.diffmeth.fem.gr, meth.diff < 0)
DMR.diffmeth.mal.gr.rename.hypo <- subset(DMR.diffmeth.mal.gr, meth.diff < 0)

DMS.diffmeth.all.gr.rename.hypo <- subset(DMS.diffmeth.all.gr, meth.diff < 0)
DMS.diffmeth.fem.gr.rename.hypo <- subset(DMS.diffmeth.fem.gr, meth.diff < 0)
DMS.diffmeth.mal.gr.rename.hypo<- subset(DMS.diffmeth.mal.gr, meth.diff < 0)

## Annotate ##
#annotate with gene parts 
DMR.all.anno <- annotateWithGeneParts(DMR.diffmeth.all.gr.rename, ref.anno)
DMR.all.anno.hyper <- annotateWithGeneParts(DMR.diffmeth.all.gr.rename.hyper, ref.anno)
DMR.all.anno.hypo <- annotateWithGeneParts(DMR.diffmeth.all.gr.rename.hypo, ref.anno)
DMS.all.anno <- annotateWithGeneParts(DMS.diffmeth.all.gr.rename, ref.anno)
DMS.all.anno.hyper <- annotateWithGeneParts(DMS.diffmeth.all.gr.rename.hyper, ref.anno)
DMS.all.anno.hypo <- annotateWithGeneParts(DMS.diffmeth.all.gr.rename.hypo, ref.anno)
CpG.all.anno <- annotateWithGeneParts(CpG.all.gr.rename, ref.anno)
regions.all.anno <- annotateWithGeneParts(regions.all.gr.rename, ref.anno)

DMR.fem.anno <- annotateWithGeneParts(DMR.diffmeth.fem.gr.rename, ref.anno)
DMR.fem.anno.hyper <- annotateWithGeneParts(DMR.diffmeth.fem.gr.rename.hyper, ref.anno)
DMR.fem.anno.hypo <- annotateWithGeneParts(DMR.diffmeth.fem.gr.rename.hyper, ref.anno)
DMS.fem.anno <- annotateWithGeneParts(DMS.diffmeth.fem.gr.rename, ref.anno)
DMS.fem.anno.hyper <- annotateWithGeneParts(DMS.diffmeth.fem.gr.rename.hypo, ref.anno)
DMS.fem.anno.hypo <- annotateWithGeneParts(DMS.diffmeth.fem.gr.rename.hypo, ref.anno)
CpG.fem.anno <- annotateWithGeneParts(CpG.fem.gr.rename, ref.anno)
regions.fem.anno <- annotateWithGeneParts(regions.fem.gr.rename, ref.anno)

DMR.mal.anno <- annotateWithGeneParts(DMR.diffmeth.mal.gr.rename, ref.anno)
DMR.mal.anno.hyper <- annotateWithGeneParts(DMR.diffmeth.mal.gr.rename.hyper, ref.anno)
DMR.mal.anno.hypo <- annotateWithGeneParts(DMR.diffmeth.mal.gr.rename.hypo, ref.anno)
DMS.mal.anno <- annotateWithGeneParts(DMS.diffmeth.mal.gr.rename, ref.anno)
DMS.mal.anno.hyper <- annotateWithGeneParts(DMS.diffmeth.mal.gr.rename.hyper, ref.anno)
DMS.mal.anno.hypo <- annotateWithGeneParts(DMS.diffmeth.mal.gr.rename.hypo, ref.anno)
CpG.mal.anno <- annotateWithGeneParts(CpG.mal.gr.rename, ref.anno)
regions.mal.anno <- annotateWithGeneParts(regions.mal.gr.rename, ref.anno)

#get percentage of DMRs that overlap with different features 
DMR.all.ann.stats.perc <- getTargetAnnotationStats(DMR.all.anno, percentage = TRUE, precedence = TRUE)
DMR.all.ann.stats.perc 
DMS.all.ann.stats.perc <- getTargetAnnotationStats(DMS.all.anno, percentage = TRUE, precedence = TRUE)
DMS.all.ann.stats.perc
CpG.all.ann.stats.perc <- getTargetAnnotationStats(CpG.all.anno, percentage = TRUE, precedence = TRUE)
CpG.all.ann.stats.perc
regions.all.ann.stats.perc <- getTargetAnnotationStats(regions.all.anno, percentage = TRUE, precedence = TRUE)
regions.all.ann.stats.perc

DMR.all.ann.stats.num <- getTargetAnnotationStats(DMR.all.anno, percentage = FALSE, precedence = TRUE)
DMR.all.ann.stats.num
DMS.all.ann.stats.num <- getTargetAnnotationStats(DMS.all.anno, percentage = FALSE, precedence = TRUE)
DMS.all.ann.stats.num
CpG.all.ann.stats.num <- getTargetAnnotationStats(CpG.all.anno, percentage = FALSE, precedence = TRUE)
CpG.all.ann.stats.num
regions.all.ann.stats.num <- getTargetAnnotationStats(regions.all.anno, percentage = FALSE, precedence = TRUE)
regions.all.ann.stats.num

DMR.fem.ann.stats.perc <- getTargetAnnotationStats(DMR.fem.anno, percentage = TRUE, precedence = TRUE)
DMR.fem.ann.stats.perc 
DMS.fem.ann.stats.perc <- getTargetAnnotationStats(DMS.fem.anno, percentage = TRUE, precedence = TRUE)
DMS.fem.ann.stats.perc
CpG.fem.ann.stats.perc <- getTargetAnnotationStats(CpG.fem.anno, percentage = TRUE, precedence = TRUE)
CpG.fem.ann.stats.perc
regions.fem.ann.stats.perc <- getTargetAnnotationStats(regions.fem.anno, percentage = TRUE, precedence = TRUE)
regions.fem.ann.stats.perc

DMR.fem.ann.stats.num <- getTargetAnnotationStats(DMR.fem.anno, percentage = FALSE, precedence = TRUE)
DMR.fem.ann.stats.num
DMS.fem.ann.stats.num <- getTargetAnnotationStats(DMS.fem.anno, percentage = FALSE, precedence = TRUE)
DMS.fem.ann.stats.num
CpG.fem.ann.stats.num <- getTargetAnnotationStats(CpG.fem.anno, percentage = FALSE, precedence = TRUE)
CpG.fem.ann.stats.num
regions.fem.ann.stats.num <- getTargetAnnotationStats(regions.fem.anno, percentage = FALSE, precedence = TRUE)
regions.fem.ann.stats.num

DMR.mal.ann.stats.perc <- getTargetAnnotationStats(DMR.mal.anno, percentage = TRUE, precedence = TRUE)
DMR.mal.ann.stats.perc 
DMS.mal.ann.stats.perc <- getTargetAnnotationStats(DMS.mal.anno, percentage = TRUE, precedence = TRUE)
DMS.mal.ann.stats.perc
CpG.mal.ann.stats.perc <- getTargetAnnotationStats(CpG.mal.anno, percentage = TRUE, precedence = TRUE)
CpG.mal.ann.stats.perc
regions.mal.ann.stats.perc <- getTargetAnnotationStats(regions.mal.anno, percentage = TRUE, precedence = TRUE)
regions.mal.ann.stats.perc

DMR.mal.ann.stats.num <- getTargetAnnotationStats(DMR.mal.anno, percentage = FALSE, precedence = TRUE)
DMR.mal.ann.stats.num
DMS.mal.ann.stats.num <- getTargetAnnotationStats(DMS.mal.anno, percentage = FALSE, precedence = TRUE)
DMS.mal.ann.stats.num
CpG.mal.ann.stats.num <- getTargetAnnotationStats(CpG.mal.anno, percentage = FALSE, precedence = TRUE)
CpG.mal.ann.stats.num
regions.mal.ann.stats.num <- getTargetAnnotationStats(regions.mal.anno, percentage = FALSE, precedence = TRUE)
regions.mal.ann.stats.num

#get nearest TSS for DMS and DMRs
DMR.all.tss <- getAssociationWithTSS(DMR.all.anno)
DMR.fem.tss <- getAssociationWithTSS(DMR.fem.anno)
DMR.mal.tss <- getAssociationWithTSS(DMR.mal.anno)

DMR.all.tss.hyper <- getAssociationWithTSS(DMR.all.anno.hyper)
DMR.fem.tss.hyper <- getAssociationWithTSS(DMR.fem.anno.hyper)
DMR.mal.tss.hyper <- getAssociationWithTSS(DMR.mal.anno.hyper)

DMR.all.tss.hypo <- getAssociationWithTSS(DMR.all.anno.hypo)
DMR.fem.tss.hypo <- getAssociationWithTSS(DMR.fem.anno.hypo)
DMR.mal.tss.hypo <- getAssociationWithTSS(DMR.mal.anno.hypo)

DMS.all.tss <- getAssociationWithTSS(DMS.all.anno)
DMS.fem.tss <- getAssociationWithTSS(DMS.fem.anno)
DMS.mal.tss <- getAssociationWithTSS(DMS.mal.anno)

DMS.all.tss.hyper <- getAssociationWithTSS(DMS.all.anno.hyper)
DMS.fem.tss.hyper <- getAssociationWithTSS(DMS.fem.anno.hyper)
DMS.mal.tss.hyper <- getAssociationWithTSS(DMS.mal.anno.hyper)

DMS.all.tss.hypo <- getAssociationWithTSS(DMS.all.anno.hypo)
DMS.fem.tss.hypo <- getAssociationWithTSS(DMS.fem.anno.hypo)
DMS.mal.tss.hypo <- getAssociationWithTSS(DMS.mal.anno.hypo)

## Save data ## 
saveRDS(DMR.all.anno, file = "./anno_res/DMR_anno_all.RDS")
saveRDS(DMR.all.anno.hyper, file = "./anno_res/DMR_anno_all_hyper.RDS")
saveRDS(DMR.all.anno.hypo, file = "./anno_res/DMR_anno_all_hypo.RDS")
saveRDS(DMS.all.anno, file = "./anno_res/DMS_anno_all.RDS")
saveRDS(DMS.all.anno.hyper, file = "./anno_res/DMS_anno_all_hyper.RDS")
saveRDS(DMS.all.anno.hypo, file = "./anno_res/DMS_anno_all_hypo.RDS")
saveRDS(CpG.all.anno, file = "./anno_res/CpG_anno_all.RDS")
saveRDS(regions.all.anno, file = "./anno_res/regions_anno_all.RDS")

saveRDS(DMR.fem.anno, file = "./anno_res/DMR_anno_fem.RDS")
saveRDS(DMR.fem.anno.hyper, file = "./anno_res/DMR_anno_fem_hyper.RDS")
saveRDS(DMR.fem.anno.hypo, file = "./anno_res/DMR_anno_fem_hypo.RDS")
saveRDS(DMS.fem.anno, file = "./anno_res/DMS_anno_fem.RDS")
saveRDS(DMS.fem.anno.hyper, file = "./anno_res/DMS_anno_fem_hyper.RDS")
saveRDS(DMS.fem.anno.hypo, file = "./anno_res/DMS_anno_fem_hypo.RDS")
saveRDS(CpG.fem.anno, file = "./anno_res/CpG_anno_fem.RDS")
saveRDS(regions.fem.anno, file = "./anno_res/regions_anno_fem.RDS")

saveRDS(DMR.mal.anno, file = "./anno_res/DMR_anno_mal.RDS")
saveRDS(DMR.mal.anno.hyper, file = "./anno_res/DMR_anno_mal_hyper.RDS")
saveRDS(DMR.mal.anno.hypo, file = "./anno_res/DMR_anno_mal_hypo.RDS")
saveRDS(DMS.mal.anno, file = "./anno_res/DMS_anno_mal.RDS")
saveRDS(DMS.mal.anno.hyper, file = "./anno_res/DMS_anno_mal_hyper.RDS")
saveRDS(DMS.mal.anno.hypo, file = "./anno_res/DMS_anno_mal_hypo.RDS")
saveRDS(CpG.mal.anno, file = "./anno_res/CpG_anno_mal.RDS")
saveRDS(regions.mal.anno, file = "./anno_res/regions_anno_mal.RDS")

saveRDS(DMR.all.ann.stats.perc, file = "./anno_res/DMR_annStats_perc_all.RDS")
saveRDS(DMS.all.ann.stats.perc, file = "./anno_res/DMS_annStats_perc_all.RDS")
saveRDS(CpG.all.ann.stats.perc, file = "./anno_res/CpG_annStats_perc_all.RDS")
saveRDS(regions.all.ann.stats.perc, file = "./anno_res/regions_annStats_perc_all.RDS")

saveRDS(DMR.fem.ann.stats.perc, file = "./anno_res/DMR_annStats_perc_fem.RDS")
saveRDS(DMS.fem.ann.stats.perc, file = "./anno_res/DMS_annStats_perc_fem.RDS")
saveRDS(CpG.fem.ann.stats.perc, file = "./anno_res/CpG_annStats_perc_fem.RDS")
saveRDS(regions.fem.ann.stats.perc, file = "./anno_res/regions_annStats_perc_fem.RDS")

saveRDS(DMR.mal.ann.stats.perc, file = "./anno_res/DMR_annStats_perc_mal.RDS")
saveRDS(DMS.mal.ann.stats.perc, file = "./anno_res/DMS_annStats_perc_mal.RDS")
saveRDS(CpG.mal.ann.stats.perc, file = "./anno_res/CpG_annStats_perc_mal.RDS")
saveRDS(regions.mal.ann.stats.perc, file = "./anno_res/regions_annStats_perc_mal.RDS")

saveRDS(DMR.all.ann.stats.num, file = "./anno_res/DMR_annStats_num_all.RDS")
saveRDS(DMS.all.ann.stats.num, file = "./anno_res/DMS_annStats_num_all.RDS")
saveRDS(CpG.all.ann.stats.num, file = "./CpG_annStats_num_all.RDS")
saveRDS(regions.all.ann.stats.num, file = "./anno_res/regions_annStats_num_all.RDS")

saveRDS(DMR.fem.ann.stats.num, file = "./anno_res/DMR_annStats_num_fem.RDS")
saveRDS(DMS.fem.ann.stats.num, file = "./anno_res/DMS_annStats_num_fem.RDS")
saveRDS(CpG.fem.ann.stats.num, file = "./anno_res/CpG_annStats_num_fem.RDS")
saveRDS(regions.fem.ann.stats.num, file = "./anno_res/regions_annStats_num_fem.RDS")

saveRDS(DMR.mal.ann.stats.num, file = "./anno_res/DMR_annStats_num_mal.RDS")
saveRDS(DMS.mal.ann.stats.num, file = "./anno_res/DMS_annStats_num_mal.RDS")
saveRDS(CpG.mal.ann.stats.num, file = "./anno_res/CpG_annStats_num_mal.RDS")
saveRDS(regions.mal.ann.stats.num, file = "./anno_res/regions_annStats_num_mal.RDS")

saveRDS(DMR.all.tss, file = "./anno_res/DMR_TSS_all.RDS")
saveRDS(DMR.all.tss.hyper, file = "./anno_res/DMR_TSS_all_hyper.RDS")
saveRDS(DMR.all.tss.hypo, file = "./anno_res/DMR_TSS_all_hypo.RDS")
saveRDS(DMS.all.tss, file = "./anno_res/DMS_TSS_all.RDS")
saveRDS(DMS.all.tss.hyper, file = "./anno_res/DMS_TSS_all_hyper.RDS")
saveRDS(DMS.all.tss.hypo, file = "./anno_res/DMS_TSS_all_hypo.RDS")

saveRDS(DMR.fem.tss, file = "./anno_res/DMR_TSS_fem.RDS")
saveRDS(DMR.fem.tss.hyper, file = "./anno_res/DMR_TSS_fem_hyper.RDS")
saveRDS(DMR.fem.tss.hypo, file = "./anno_res/DMR_TSS_fem_hypo.RDS")
saveRDS(DMS.fem.tss, file = "./anno_res/DMS_TSS_fem.RDS")
saveRDS(DMS.fem.tss.hyper, file = "./anno_res/DMS_TSS_fem_hyper.RDS")
saveRDS(DMS.fem.tss.hypo, file = "./anno_res/DMS_TSS_fem_hypo.RDS")

saveRDS(DMR.mal.tss, file = "./anno_res/DMR_TSS_mal.RDS")
saveRDS(DMR.mal.tss.hyper, file = "./anno_res/DMR_TSS_mal_hyper.RDS")
saveRDS(DMR.mal.tss.hypo, file = "./anno_res/DMR_TSS_mal_hypo.RDS")
saveRDS(DMS.mal.tss, file = "./anno_res/DMS_TSS_mal.RDS")
saveRDS(DMS.mal.tss.hyper, file = "./anno_res/DMS_TSS_mal_hyper.RDS")
saveRDS(DMS.mal.tss.hypo, file = "./anno_res/DMS_TSS_mal_hypo.RDS")