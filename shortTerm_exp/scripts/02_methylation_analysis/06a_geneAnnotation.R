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

setwd("/scratch/janayfox/guppyWGBS/methylKit/st")

## Load in data and reformat ##
ref.anno <- readTranscriptFeatures("/scratch/janayfox/guppyWGBS/annotation/guppy.sorted.bed", unique.prom = FALSE, up.flank = 1500, down.flank = 500)

DMR.diffmeth.05h <- readRDS("./05h/DMR_res/DMRdiffMeth_05h_5X.RDS")
DMR.diffmeth.1h <- readRDS("./1h/DMR_res/DMRdiffMeth_1h_5X.RDS")
DMR.diffmeth.4h <- readRDS("./4h/DMR_res/DMRdiffMeth_4h_5X.RDS")
DMR.diffmeth.24h <- readRDS("./24h/DMR_res/DMRdiffMeth_24h_5X.RDS")
DMR.diffmeth.72h <- readRDS("./72h/DMR_res/DMRdiffMeth_72h_5X.RDS")
DMR.diffmeth.all <- readRDS("./all/DMR_res/DMRdiffMeth_all_5X.RDS")

DMR.diffmeth.05h.fem <- readRDS("./05h/DMR_res/DMRdiffMeth_05h_fem_5X.RDS")
DMR.diffmeth.1h.fem <- readRDS("./1h/DMR_res/DMRdiffMeth_1h_fem_5X.RDS")
DMR.diffmeth.4h.fem <- readRDS("./4h/DMR_res/DMRdiffMeth_4h_fem_5X.RDS")
DMR.diffmeth.24h.fem <- readRDS("./24h/DMR_res/DMRdiffMeth_24h_fem_5X.RDS")
DMR.diffmeth.72h.fem <- readRDS("./72h/DMR_res/DMRdiffMeth_72h_fem_5X.RDS")
DMR.diffmeth.all.fem <- readRDS("./all/DMR_res/DMRdiffMeth_all_fem_5X.RDS")

DMR.diffmeth.05h.mal <- readRDS("./05h/DMR_res/DMRdiffMeth_05h_mal_5X.RDS")
DMR.diffmeth.1h.mal <- readRDS("./1h/DMR_res/DMRdiffMeth_1h_mal_5X.RDS")
DMR.diffmeth.4h.mal <- readRDS("./4h/DMR_res/DMRdiffMeth_4h_mal_5X.RDS")
DMR.diffmeth.24h.mal <- readRDS("./24h/DMR_res/DMRdiffMeth_24h_mal_5X.RDS")
DMR.diffmeth.72h.mal <- readRDS("./72h/DMR_res/DMRdiffMeth_72h_mal_5X.RDS")
DMR.diffmeth.all.mal <- readRDS("./all/DMR_res/DMRdiffMeth_all_mal_5X.RDS")

DMS.diffmeth.05h <- readRDS("./05h/DMS_res/shortterm_DMSDiffMeth_05h_5X.RDS")
DMS.diffmeth.1h <- readRDS("./1h/DMS_res/shortterm_DMSDiffMeth_1h_5X.RDS")
DMS.diffmeth.4h <- readRDS("./4h/DMS_res/shortterm_DMSDiffMeth_4h_5X.RDS")
DMS.diffmeth.24h <- readRDS("./24h/DMS_res/shortterm_DMSDiffMeth_24h_5X.RDS")
DMS.diffmeth.72h <- readRDS("./72h/DMS_res/shortterm_DMSDiffMeth_72h_5X.RDS")
DMS.diffmeth.all <- readRDS("./all/DMS_res/shortterm_DMSDiffMeth_all_5X.RDS")

DMS.diffmeth.05h.fem <- readRDS("./05h/DMS_res/shortterm_DMSDiffMeth_05h_fem_5X.RDS")
DMS.diffmeth.1h.fem <- readRDS("./1h/DMS_res/shortterm_DMSDiffMeth_1h_fem_5X.RDS")
DMS.diffmeth.4h.fem <- readRDS("./4h/DMS_res/shortterm_DMSDiffMeth_4h_fem_5X.RDS")
DMS.diffmeth.24h.fem <- readRDS("./24h/DMS_res/shortterm_DMSDiffMeth_24h_fem_5X.RDS")
DMS.diffmeth.72h.fem <- readRDS("./72h/DMS_res/shortterm_DMSDiffMeth_72h_fem_5X.RDS")
DMS.diffmeth.all.fem <- readRDS("./all/DMS_res/shortterm_DMSDiffMeth_all_fem_5X.RDS")

DMS.diffmeth.05h.mal <- readRDS("./05h/DMS_res/shortterm_DMSDiffMeth_05h_mal_5X.RDS")
DMS.diffmeth.1h.mal <- readRDS("./1h/DMS_res/shortterm_DMSDiffMeth_1h_mal_5X.RDS")
DMS.diffmeth.4h.mal <- readRDS("./4h/DMS_res/shortterm_DMSDiffMeth_4h_mal_5X.RDS")
DMS.diffmeth.24h.mal <- readRDS("./24h/DMS_res/shortterm_DMSDiffMeth_24h_mal_5X.RDS")
DMS.diffmeth.72h.mal <- readRDS("./72h/DMS_res/shortterm_DMSDiffMeth_72h_mal_5X.RDS")
DMS.diffmeth.all.mal <- readRDS("./all/DMS_res/shortterm_DMSDiffMeth_all_mal_5X.RDS")

CpG.05h <- readRDS("./05h/DMS_res/shortterm_DMSmyDiff_05h_5X.RDS")
CpG.1h <- readRDS("./1h/DMS_res/shortterm_DMSmyDiff_1h_5X.RDS")
CpG.4h <- readRDS("./4h/DMS_res/shortterm_DMSmyDiff_4h_5X.RDS")
CpG.24h <- readRDS("./24h/DMS_res/shortterm_DMSmyDiff_24h_5X.RDS")
CpG.72h <- readRDS("./72h/DMS_res/shortterm_DMSmyDiff_72h_5X.RDS")
CpG.all <- readRDS("./all/DMS_res/shortterm_DMSmyDiff_all_5X.RDS")

CpG.05h.fem <- readRDS("./05h/DMS_res/shortterm_DMSmyDiff_05h_fem_5X.RDS")
CpG.1h.fem <- readRDS("./1h/DMS_res/shortterm_DMSmyDiff_1h_fem_5X.RDS")
CpG.4h.fem <- readRDS("./4h/DMS_res/shortterm_DMSmyDiff_4h_fem_5X.RDS")
CpG.24h.fem <- readRDS("./24h/DMS_res/shortterm_DMSmyDiff_24h_fem_5X.RDS")
CpG.72h.fem <- readRDS("./72h/DMS_res/shortterm_DMSmyDiff_72h_fem_5X.RDS")
CpG.all.fem <- readRDS("./all/DMS_res/shortterm_DMSmyDiff_all_fem_5X.RDS")

CpG.05h.mal <- readRDS("./05h/DMS_res/shortterm_DMSmyDiff_05h_mal_5X.RDS")
CpG.1h.mal <- readRDS("./1h/DMS_res/shortterm_DMSmyDiff_1h_mal_5X.RDS")
CpG.4h.mal <- readRDS("./4h/DMS_res/shortterm_DMSmyDiff_4h_mal_5X.RDS")
CpG.24h.mal <- readRDS("./24h/DMS_res/shortterm_DMSmyDiff_24h_mal_5X.RDS")
CpG.72h.mal <- readRDS("./72h/DMS_res/shortterm_DMSmyDiff_72h_mal_5X.RDS")
CpG.all.mal <- readRDS("./all/DMS_res/shortterm_DMSmyDiff_all_mal_5X.RDS")

#change chromosome names to match 
#create function that renames chromosomes and converts to G ranges 
renameChr <- function(obj) {
  obj.gr <- as(obj, "GRanges")
  gr.obj.rename <- renameSeqlevels(obj.gr, c(NC_024331.1="LG1", NC_024332.1="LG2",
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

renameChr.noXY <- function(obj) {
  obj.gr <- as(obj, "GRanges")
  gr.obj.rename <- renameSeqlevels(obj.gr, c(NC_024331.1="LG1", NC_024332.1="LG2",
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
DMR.diffmeth.05h.gr.rename <- renameChr.noXY(DMR.diffmeth.05h)
DMR.diffmeth.1h.gr.rename <- renameChr.noXY(DMR.diffmeth.1h)
DMR.diffmeth.4h.gr.rename <- renameChr.noXY(DMR.diffmeth.4h)
DMR.diffmeth.24h.gr.rename <- renameChr.noXY(DMR.diffmeth.24h)
DMR.diffmeth.72h.gr.rename <- renameChr.noXY(DMR.diffmeth.72h)
DMR.diffmeth.all.gr.rename <- renameChr.noXY(DMR.diffmeth.all)

DMR.diffmeth.05h.fem.gr.rename <- renameChr(DMR.diffmeth.05h.fem)
DMR.diffmeth.1h.fem.gr.rename <- renameChr(DMR.diffmeth.1h.fem)
DMR.diffmeth.4h.fem.gr.rename <- renameChr(DMR.diffmeth.4h.fem)
DMR.diffmeth.24h.fem.gr.rename <- renameChr(DMR.diffmeth.24h.fem)
DMR.diffmeth.72h.fem.gr.rename <- renameChr(DMR.diffmeth.72h.fem)
DMR.diffmeth.all.fem.gr.rename <- renameChr(DMR.diffmeth.all.fem)

DMR.diffmeth.05h.mal.gr.rename <- renameChr(DMR.diffmeth.05h.mal)
DMR.diffmeth.1h.mal.gr.rename <- renameChr(DMR.diffmeth.1h.mal)
DMR.diffmeth.4h.mal.gr.rename <- renameChr(DMR.diffmeth.4h.mal)
DMR.diffmeth.24h.mal.gr.rename <- renameChr(DMR.diffmeth.24h.mal)
DMR.diffmeth.72h.mal.gr.rename <- renameChr(DMR.diffmeth.72h.mal)
DMR.diffmeth.all.mal.gr.rename <- renameChr(DMR.diffmeth.all.mal)

DMS.diffmeth.05h.gr.rename <- renameChr.noXY(DMS.diffmeth.05h)
DMS.diffmeth.1h.gr.rename <- renameChr.noXY(DMS.diffmeth.1h)
DMS.diffmeth.4h.gr.rename <- renameChr.noXY(DMS.diffmeth.4h)
DMS.diffmeth.24h.gr.rename <- renameChr.noXY(DMS.diffmeth.24h)
DMS.diffmeth.72h.gr.rename <- renameChr.noXY(DMS.diffmeth.72h)
DMS.diffmeth.all.gr.rename <- renameChr.noXY(DMS.diffmeth.all)

DMS.diffmeth.05h.fem.gr.rename <- renameChr(DMS.diffmeth.05h.fem)
DMS.diffmeth.1h.fem.gr.rename <- renameChr(DMS.diffmeth.1h.fem)
DMS.diffmeth.4h.fem.gr.rename <- renameChr(DMS.diffmeth.4h.fem)
DMS.diffmeth.24h.fem.gr.rename <- renameChr(DMS.diffmeth.24h.fem)
DMS.diffmeth.72h.fem.gr.rename <- renameChr(DMS.diffmeth.72h.fem)
DMS.diffmeth.all.fem.gr.rename <- renameChr(DMS.diffmeth.all.fem)

DMS.diffmeth.05h.mal.gr.rename <- renameChr(DMS.diffmeth.05h.mal)
DMS.diffmeth.1h.mal.gr.rename <- renameChr(DMS.diffmeth.1h.mal)
DMS.diffmeth.4h.mal.gr.rename <- renameChr(DMS.diffmeth.4h.mal)
DMS.diffmeth.24h.mal.gr.rename <- renameChr(DMS.diffmeth.24h.mal)
DMS.diffmeth.72h.mal.gr.rename <- renameChr(DMS.diffmeth.72h.mal)
DMS.diffmeth.all.mal.gr.rename <- renameChr(DMS.diffmeth.all.mal)

CpG.05h.gr.rename <- renameChr.noXY(CpG.05h)
CpG.1h.gr.rename <- renameChr.noXY(CpG.1h)
CpG.4h.gr.rename <- renameChr.noXY(CpG.4h)
CpG.24h.gr.rename <- renameChr.noXY(CpG.24h)
CpG.72h.gr.rename <- renameChr.noXY(CpG.72h)
CpG.all.gr.rename <- renameChr.noXY(CpG.all)

CpG.05h.fem.gr.rename <- renameChr(CpG.05h.fem)
CpG.1h.fem.gr.rename <- renameChr(CpG.1h.fem)
CpG.4h.fem.gr.rename <- renameChr(CpG.4h.fem)
CpG.24h.fem.gr.rename <- renameChr(CpG.24h.fem)
CpG.72h.fem.gr.rename <- renameChr(CpG.72h.fem)
CpG.all.fem.gr.rename <- renameChr(CpG.all.fem)

CpG.05h.mal.gr.rename <- renameChr(CpG.05h.mal)
CpG.1h.mal.gr.rename <- renameChr(CpG.1h.mal)
CpG.4h.mal.gr.rename <- renameChr(CpG.4h.mal)
CpG.24h.mal.gr.rename <- renameChr(CpG.24h.mal)
CpG.72h.mal.gr.rename <- renameChr(CpG.72h.mal)
CpG.all.mal.gr.rename <- renameChr(CpG.all.mal)

## Annotate ##
anno.func <- function(DMR, DMS, CpG, 
                      DMR_anno_name, DMS_anno_name, CpG_anno_name,
                      DMR_perc_name, DMS_perc_name, CpG_perc_name,
                      DMR_num_name, DMS_num_name, CpG_num_name,
                      DMR_tss_name, DMS_tss_name) {
  #annotate with gene parts 
  DMR.anno <- annotateWithGeneParts(DMR, ref.anno)
  DMS.anno <- annotateWithGeneParts(DMS, ref.anno)
  CpG.anno <- annotateWithGeneParts(CpG, ref.anno)

  #get percentage of DMRs that overlap with different features 
  DMR.ann.stats.perc <- getTargetAnnotationStats(DMR.anno, percentage = TRUE, precedence = TRUE)
  DMS.ann.stats.perc <- getTargetAnnotationStats(DMS.anno, percentage = TRUE, precedence = TRUE)
  CpG.ann.stats.perc <- getTargetAnnotationStats(CpG.anno, percentage = TRUE, precedence = TRUE)

  DMR.ann.stats.num <- getTargetAnnotationStats(DMR.anno, percentage = FALSE, precedence = TRUE)
  DMS.ann.stats.num <- getTargetAnnotationStats(DMS.anno, percentage = FALSE, precedence = TRUE)
  CpG.ann.stats.num <- getTargetAnnotationStats(CpG.anno, percentage = FALSE, precedence = TRUE)
  
  #get nearest TSS for DMS and DMRs
  DMR.tss <- getAssociationWithTSS(DMR.anno)
  DMS.tss <- getAssociationWithTSS(DMS.anno)

  ## Save data ## 
  saveRDS(DMR.anno, file = DMR_anno_name)
  saveRDS(DMS.anno, file = DMS_anno_name)
  saveRDS(CpG.anno, file = CpG_anno_name)

  saveRDS(DMR.ann.stats.perc, file = DMR_perc_name)
  saveRDS(DMS.ann.stats.perc, file = DMS_perc_name)
  saveRDS(CpG.ann.stats.perc, file = CpG_perc_name)

  saveRDS(DMR.ann.stats.num, file = DMR_num_name)
  saveRDS(DMS.ann.stats.num, file = DMS_num_name)
  saveRDS(CpG.ann.stats.num, file = CpG_num_name)

  saveRDS(DMR.tss, file = DMR_tss_name)
  saveRDS(DMS.tss, file = DMS_tss_name)
}

anno.func(DMR.diffmeth.05h.gr.rename, DMS.diffmeth.05h.gr.rename, CpG.05h.gr.rename,
          "./05h/anno_res/DMR_anno_05h.RDS", "./05h/anno_res/DMS_anno_05h.RDS", "./05h/anno_res/CpG_anno_05h.RDS",
          "./05h/anno_res/DMR_annStats_perc_05h.RDS", "./05h/anno_res/DMS_annStats_perc_05h.RDS", "./05h/anno_res/CpG_annStats_perc_05h.RDS",
          "./05h/anno_res/DMR_annStats_num_05h.RDS", "./05h/anno_res/DMS_annStats_num_05h.RDS", "./05h/anno_res/CpG_annStats_num_05h.RDS",
          "./05h/anno_res/DMR_TSS_05h.RDS", "./05h/anno_res/DMS_TSS_05h.RDS")

anno.func(DMR.diffmeth.05h.fem.gr.rename, DMS.diffmeth.05h.fem.gr.rename, CpG.05h.fem.gr.rename,
          "./05h/anno_res/DMR_anno_05h_fem.RDS", "./05h/anno_res/DMS_anno_05h_fem.RDS", "./05h/anno_res/CpG_anno_05h_fem.RDS",
          "./05h/anno_res/DMR_annStats_perc_05h_fem.RDS", "./05h/anno_res/DMS_annStats_perc_05h_fem.RDS", "./05h/anno_res/CpG_annStats_perc_05h_fem.RDS",
          "./05h/anno_res/DMR_annStats_num_05h_fem.RDS", "./05h/anno_res/DMS_annStats_num_05h_fem.RDS", "./05h/anno_res/CpG_annStats_num_05h_fem.RDS",
          "./05h/anno_res/DMR_TSS_05h_fem.RDS", "./05h/anno_res/DMS_TSS_05h_fem.RDS")

anno.func(DMR.diffmeth.05h.mal.gr.rename, DMS.diffmeth.05h.mal.gr.rename, CpG.05h.mal.gr.rename,
          "./05h/anno_res/DMR_anno_05h_mal.RDS", "./05h/anno_res/DMS_anno_05h_mal.RDS", "./05h/anno_res/CpG_anno_05h_mal.RDS",
          "./05h/anno_res/DMR_annStats_perc_05h_mal.RDS", "./05h/anno_res/DMS_annStats_perc_05h_mal.RDS", "./05h/anno_res/CpG_annStats_perc_05h_mal.RDS",
          "./05h/anno_res/DMR_annStats_num_05h_mal.RDS", "./05h/anno_res/DMS_annStats_num_05h_mal.RDS", "./05h/anno_res/CpG_annStats_num_05h_mal.RDS",
          "./05h/anno_res/DMR_TSS_05h_mal.RDS", "./05h/anno_res/DMS_TSS_05h_mal.RDS")

anno.func(DMR.diffmeth.1h.gr.rename, DMS.diffmeth.1h.gr.rename, CpG.1h.gr.rename,
          "./1h/anno_res/DMR_anno_1h.RDS", "./1h/anno_res/DMS_anno_1h.RDS", "./1h/anno_res/CpG_anno_1h.RDS",
          "./1h/anno_res/DMR_annStats_perc_1h.RDS", "./1h/anno_res/DMS_annStats_perc_1h.RDS", "./1h/anno_res/CpG_annStats_perc_1h.RDS",
          "./1h/anno_res/DMR_annStats_num_1h.RDS", "./1h/anno_res/DMS_annStats_num_1h.RDS", "./1h/anno_res/CpG_annStats_num_1h.RDS",
          "./1h/anno_res/DMR_TSS_1h.RDS", "./1h/anno_res/DMS_TSS_1h.RDS")

anno.func(DMR.diffmeth.1h.fem.gr.rename, DMS.diffmeth.1h.fem.gr.rename, CpG.1h.fem.gr.rename,
          "./1h/anno_res/DMR_anno_1h_fem.RDS", "./1h/anno_res/DMS_anno_1h_fem.RDS", "./1h/anno_res/CpG_anno_1h_fem.RDS",
          "./1h/anno_res/DMR_annStats_perc_1h_fem.RDS", "./1h/anno_res/DMS_annStats_perc_1h_fem.RDS", "./1h/anno_res/CpG_annStats_perc_1h_fem.RDS",
          "./1h/anno_res/DMR_annStats_num_1h_fem.RDS", "./1h/anno_res/DMS_annStats_num_1h_fem.RDS", "./1h/anno_res/CpG_annStats_num_1h_fem.RDS",
          "./1h/anno_res/DMR_TSS_1h_fem.RDS", "./1h/anno_res/DMS_TSS_1h_fem.RDS")

anno.func(DMR.diffmeth.1h.mal.gr.rename, DMS.diffmeth.1h.mal.gr.rename, CpG.1h.mal.gr.rename,
          "./1h/anno_res/DMR_anno_1h_mal.RDS", "./1h/anno_res/DMS_anno_1h_mal.RDS", "./1h/anno_res/CpG_anno_1h_mal.RDS",
          "./1h/anno_res/DMR_annStats_perc_1h_mal.RDS", "./1h/anno_res/DMS_annStats_perc_1h_mal.RDS", "./1h/anno_res/CpG_annStats_perc_1h_mal.RDS",
          "./1h/anno_res/DMR_annStats_num_1h_mal.RDS", "./1h/anno_res/DMS_annStats_num_1h_mal.RDS", "./1h/anno_res/CpG_annStats_num_1h_mal.RDS",
          "./1h/anno_res/DMR_TSS_1h_mal.RDS", "./1h/anno_res/DMS_TSS_1h_mal.RDS")

anno.func(DMR.diffmeth.4h.gr.rename, DMS.diffmeth.4h.gr.rename, CpG.4h.gr.rename,
          "./4h/anno_res/DMR_anno_4h.RDS", "./4h/anno_res/DMS_anno_4h.RDS", "./4h/anno_res/CpG_anno_4h.RDS",
          "./4h/anno_res/DMR_annStats_perc_4h.RDS", "./4h/anno_res/DMS_annStats_perc_4h.RDS", "./4h/anno_res/CpG_annStats_perc_4h.RDS",
          "./4h/anno_res/DMR_annStats_num_4h.RDS", "./4h/anno_res/DMS_annStats_num_4h.RDS", "./4h/anno_res/CpG_annStats_num_4h.RDS",
          "./4h/anno_res/DMR_TSS_4h.RDS", "./4h/anno_res/DMS_TSS_4h.RDS")

anno.func(DMR.diffmeth.4h.fem.gr.rename, DMS.diffmeth.4h.fem.gr.rename, CpG.4h.fem.gr.rename,
          "./4h/anno_res/DMR_anno_4h_fem.RDS", "./4h/anno_res/DMS_anno_4h_fem.RDS", "./4h/anno_res/CpG_anno_4h_fem.RDS",
          "./4h/anno_res/DMR_annStats_perc_4h_fem.RDS", "./4h/anno_res/DMS_annStats_perc_4h_fem.RDS", "./4h/anno_res/CpG_annStats_perc_4h_fem.RDS",
          "./4h/anno_res/DMR_annStats_num_4h_fem.RDS", "./4h/anno_res/DMS_annStats_num_4h_fem.RDS", "./4h/anno_res/CpG_annStats_num_4h_fem.RDS",
          "./4h/anno_res/DMR_TSS_4h_fem.RDS", "./4h/anno_res/DMS_TSS_4h_fem.RDS")

anno.func(DMR.diffmeth.4h.mal.gr.rename, DMS.diffmeth.4h.mal.gr.rename, CpG.4h.mal.gr.rename,
          "./4h/anno_res/DMR_anno_4h_mal.RDS", "./4h/anno_res/DMS_anno_4h_mal.RDS", "./4h/anno_res/CpG_anno_4h_mal.RDS",
          "./4h/anno_res/DMR_annStats_perc_4h_mal.RDS", "./4h/anno_res/DMS_annStats_perc_4h_mal.RDS", "./4h/anno_res/CpG_annStats_perc_4h_mal.RDS",
          "./4h/anno_res/DMR_annStats_num_4h_mal.RDS", "./4h/anno_res/DMS_annStats_num_4h_mal.RDS", "./4h/anno_res/CpG_annStats_num_4h_mal.RDS",
          "./4h/anno_res/DMR_TSS_4h_mal.RDS", "./4h/anno_res/DMS_TSS_4h_mal.RDS")

anno.func(DMR.diffmeth.24h.gr.rename, DMS.diffmeth.24h.gr.rename, CpG.24h.gr.rename,
          "./24h/anno_res/DMR_anno_24h.RDS", "./24h/anno_res/DMS_anno_24h.RDS", "./24h/anno_res/CpG_anno_24h.RDS",
          "./24h/anno_res/DMR_annStats_perc_24h.RDS", "./24h/anno_res/DMS_annStats_perc_24h.RDS", "./24h/anno_res/CpG_annStats_perc_24h.RDS",
          "./24h/anno_res/DMR_annStats_num_24h.RDS", "./24h/anno_res/DMS_annStats_num_24h.RDS", "./24h/anno_res/CpG_annStats_num_24h.RDS",
          "./24h/anno_res/DMR_TSS_24h.RDS", "./24h/anno_res/DMS_TSS_24h.RDS")

anno.func(DMR.diffmeth.24h.fem.gr.rename, DMS.diffmeth.24h.fem.gr.rename, CpG.24h.fem.gr.rename,
          "./24h/anno_res/DMR_anno_24h_fem.RDS", "./24h/anno_res/DMS_anno_24h_fem.RDS", "./24h/anno_res/CpG_anno_24h_fem.RDS",
          "./24h/anno_res/DMR_annStats_perc_24h_fem.RDS", "./24h/anno_res/DMS_annStats_perc_24h_fem.RDS", "./24h/anno_res/CpG_annStats_perc_24h_fem.RDS",
          "./24h/anno_res/DMR_annStats_num_24h_fem.RDS", "./24h/anno_res/DMS_annStats_num_24h_fem.RDS", "./24h/anno_res/CpG_annStats_num_24h_fem.RDS",
          "./24h/anno_res/DMR_TSS_24h_fem.RDS", "./24h/anno_res/DMS_TSS_24h_fem.RDS")

anno.func(DMR.diffmeth.24h.mal.gr.rename, DMS.diffmeth.24h.mal.gr.rename, CpG.24h.mal.gr.rename,
          "./24h/anno_res/DMR_anno_24h_mal.RDS", "./24h/anno_res/DMS_anno_24h_mal.RDS", "./24h/anno_res/CpG_anno_24h_mal.RDS",
          "./24h/anno_res/DMR_annStats_perc_24h_mal.RDS", "./24h/anno_res/DMS_annStats_perc_24h_mal.RDS", "./24h/anno_res/CpG_annStats_perc_24h_mal.RDS",
          "./24h/anno_res/DMR_annStats_num_24h_mal.RDS", "./24h/anno_res/DMS_annStats_num_24h_mal.RDS", "./24h/anno_res/CpG_annStats_num_24h_mal.RDS",
          "./24h/anno_res/DMR_TSS_24h_mal.RDS", "./24h/anno_res/DMS_TSS_24h_mal.RDS")

anno.func(DMR.diffmeth.72h.gr.rename, DMS.diffmeth.72h.gr.rename, CpG.72h.gr.rename,
          "./72h/anno_res/DMR_anno_72h.RDS", "./72h/anno_res/DMS_anno_72h.RDS", "./72h/anno_res/CpG_anno_72h.RDS",
          "./72h/anno_res/DMR_annStats_perc_72h.RDS", "./72h/anno_res/DMS_annStats_perc_72h.RDS", "./72h/anno_res/CpG_annStats_perc_72h.RDS",
          "./72h/anno_res/DMR_annStats_num_72h.RDS", "./72h/anno_res/DMS_annStats_num_72h.RDS", "./72h/anno_res/CpG_annStats_num_72h.RDS",
          "./72h/anno_res/DMR_TSS_72h.RDS", "./72h/anno_res/DMS_TSS_72h.RDS")

anno.func(DMR.diffmeth.72h.fem.gr.rename, DMS.diffmeth.72h.fem.gr.rename, CpG.72h.fem.gr.rename,
          "./72h/anno_res/DMR_anno_72h_fem.RDS", "./72h/anno_res/DMS_anno_72h_fem.RDS", "./72h/anno_res/CpG_anno_72h_fem.RDS",
          "./72h/anno_res/DMR_annStats_perc_72h_fem.RDS", "./72h/anno_res/DMS_annStats_perc_72h_fem.RDS", "./72h/anno_res/CpG_annStats_perc_72h_fem.RDS",
          "./72h/anno_res/DMR_annStats_num_72h_fem.RDS", "./72h/anno_res/DMS_annStats_num_72h_fem.RDS", "./72h/anno_res/CpG_annStats_num_72h_fem.RDS",
          "./72h/anno_res/DMR_TSS_72h_fem.RDS", "./72h/anno_res/DMS_TSS_72h_fem.RDS")

anno.func(DMR.diffmeth.72h.mal.gr.rename, DMS.diffmeth.72h.mal.gr.rename, CpG.72h.mal.gr.rename,
          "./72h/anno_res/DMR_anno_72h_mal.RDS", "./72h/anno_res/DMS_anno_72h_mal.RDS", "./72h/anno_res/CpG_anno_72h_mal.RDS",
          "./72h/anno_res/DMR_annStats_perc_72h_mal.RDS", "./72h/anno_res/DMS_annStats_perc_72h_mal.RDS", "./72h/anno_res/CpG_annStats_perc_72h_mal.RDS",
          "./72h/anno_res/DMR_annStats_num_72h_mal.RDS", "./72h/anno_res/DMS_annStats_num_72h_mal.RDS", "./72h/anno_res/CpG_annStats_num_72h_mal.RDS",
          "./72h/anno_res/DMR_TSS_72h_mal.RDS", "./72h/anno_res/DMS_TSS_72h_mal.RDS")

anno.func(DMR.diffmeth.all.gr.rename, DMS.diffmeth.all.gr.rename, CpG.all.gr.rename,
          "./all/anno_res/DMR_anno_all.RDS", "./all/anno_res/DMS_anno_all.RDS", "./all/anno_res/CpG_anno_all.RDS",
          "./all/anno_res/DMR_annStats_perc_all.RDS", "./all/anno_res/DMS_annStats_perc_all.RDS", "./all/anno_res/CpG_annStats_perc_all.RDS",
          "./all/anno_res/DMR_annStats_num_all.RDS", "./all/anno_res/DMS_annStats_num_all.RDS", "./all/anno_res/CpG_annStats_num_all.RDS",
          "./all/anno_res/DMR_TSS_all.RDS", "./all/anno_res/DMS_TSS_all.RDS")

anno.func(DMR.diffmeth.all.fem.gr.rename, DMS.diffmeth.all.fem.gr.rename, CpG.all.fem.gr.rename,
          "./all/anno_res/DMR_anno_all_fem.RDS", "./all/anno_res/DMS_anno_all_fem.RDS", "./all/anno_res/CpG_anno_all_fem.RDS",
          "./all/anno_res/DMR_annStats_perc_all_fem.RDS", "./all/anno_res/DMS_annStats_perc_all_fem.RDS", "./all/anno_res/CpG_annStats_perc_all_fem.RDS",
          "./all/anno_res/DMR_annStats_num_all_fem.RDS", "./all/anno_res/DMS_annStats_num_all_fem.RDS", "./all/anno_res/CpG_annStats_num_all_fem.RDS",
          "./all/anno_res/DMR_TSS_all_fem.RDS", "./all/anno_res/DMS_TSS_all_fem.RDS")

anno.func(DMR.diffmeth.all.mal.gr.rename, DMS.diffmeth.all.mal.gr.rename, CpG.all.mal.gr.rename,
          "./all/anno_res/DMR_anno_all_mal.RDS", "./all/anno_res/DMS_anno_all_mal.RDS", "./all/anno_res/CpG_anno_all_mal.RDS",
          "./all/anno_res/DMR_annStats_perc_all_mal.RDS", "./all/anno_res/DMS_annStats_perc_all_mal.RDS", "./all/anno_res/CpG_annStats_perc_all_mal.RDS",
          "./all/anno_res/DMR_annStats_num_all_mal.RDS", "./all/anno_res/DMS_annStats_num_all_mal.RDS", "./all/anno_res/CpG_annStats_num_all_mal.RDS",
          "./all/anno_res/DMR_TSS_all_mal.RDS", "./all/anno_res/DMS_TSS_all_mal.RDS")