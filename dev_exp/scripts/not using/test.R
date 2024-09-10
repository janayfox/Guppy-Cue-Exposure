
#load packages
library("S4Vectors")
library("IRanges")
library("GenomicRanges")
library("methylKit")
library("genomation")

## Load in data and reformat ##
ref.anno <- readTranscriptFeatures("./data/annotation/guppy.sorted.bed", unique.prom = FALSE, up.flank = 1500, down.flank = 500)

DMR.diffmeth.all <- readRDS("./data/methylkit_res/DMR_res/DMRdiffmeth_all_5X.RDS")
DMR.diffmeth.fem <- readRDS("./data/methylkit_res/DMR_res/DMRdiffmeth_fem_5X.RDS")
DMR.diffmeth.mal <- readRDS("./data/methylkit_res/DMR_res/DMRdiffmeth_mal_5X.RDS")

DMS.diffmeth.all <- readRDS("./data/methylkit_res/DMS_res/DMSdiffmeth_all_5X.RDS")
DMS.diffmeth.fem <- readRDS("./data/methylkit_res/DMS_res/DMSdiffmeth_fem_5X.RDS")
DMS.diffmeth.mal <- readRDS("./data/methylkit_res/DMS_res/DMSdiffmeth_mal_5X.RDS")

CpG.all <- readRDS("./data/methylkit_res/DMS_res/DMSmydiff_all_5X.RDS")
CpG.fem <- readRDS("./data/methylkit_res/DMS_res/DMSmydiff_fem_5X.RDS")
CpG.mal <- readRDS("./data/methylkit_res/DMS_res/DMSmydiff_mal_5X.RDS")

regions.all <- readRDS("./data/methylkit_res/DMR_res/DMRmydiff_all_5X.RDS")
regions.fem <- readRDS("./data/methylkit_res/DMR_res/DMRmydiff_fem_5X.RDS")
regions.mal <- readRDS("./data/methylkit_res/DMR_res/DMRmydiff_mal_5X.RDS")

#change chromosome names to match 
#create function that renames chromosomes
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

#rename chromosomes in each GRanges
DMR.diffmeth.all.gr.rename <- renameChr(DMR.diffmeth.all)
DMR.diffmeth.fem.gr.rename <- renameChr(DMR.diffmeth.fem)
DMR.diffmeth.mal.gr.rename <- renameChr(DMR.diffmeth.mal)

DMS.diffmeth.all.gr.rename <- renameChr(DMS.diffmeth.all)
DMS.diffmeth.fem.gr.rename <- renameChr(DMS.diffmeth.fem)
DMS.diffmeth.mal.gr.rename <- renameChr(DMS.diffmeth.mal)

CpG.all.gr.rename <- renameChr(CpG.all)
CpG.fem.gr.rename <- renameChr(CpG.fem)
CpG.mal.gr.rename <- renameChr(CpG.mal)

regions.all.gr.rename <- renameChr(regions.all)
regions.fem.gr.rename <- renameChr(regions.fem)
regions.mal.gr.rename <- renameChr(regions.mal)

#get hyper and hypo methylated 
DMR.diffmeth.all.gr.rename.hyper <- subset(DMR.diffmeth.all.gr.rename, meth.diff > 0)
DMR.diffmeth.fem.gr.rename.hyper <- subset(DMR.diffmeth.fem.gr.rename, meth.diff > 0)
DMR.diffmeth.mal.gr.rename.hyper <- subset(DMR.diffmeth.mal.gr.rename, meth.diff > 0)

DMS.diffmeth.all.gr.rename.hyper <- subset(DMS.diffmeth.all.gr.rename, meth.diff > 0)
DMS.diffmeth.fem.gr.rename.hyper <- subset(DMS.diffmeth.fem.gr.rename, meth.diff > 0)
DMS.diffmeth.mal.gr.rename.hyper <- subset(DMS.diffmeth.mal.gr.rename, meth.diff > 0)

DMR.diffmeth.all.gr.rename.hypo <- subset(DMR.diffmeth.all.gr.rename, meth.diff < 0)
DMR.diffmeth.fem.gr.rename.hypo <- subset(DMR.diffmeth.fem.gr.rename, meth.diff < 0)
DMR.diffmeth.mal.gr.rename.hypo <- subset(DMR.diffmeth.mal.gr.rename, meth.diff < 0)

DMS.diffmeth.all.gr.rename.hypo <- subset(DMS.diffmeth.all.gr.rename, meth.diff < 0)
DMS.diffmeth.fem.gr.rename.hypo <- subset(DMS.diffmeth.fem.gr.rename, meth.diff < 0)
DMS.diffmeth.mal.gr.rename.hypo <- subset(DMS.diffmeth.mal.gr.rename, meth.diff < 0)

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
DMR.fem.anno.hypo <- annotateWithGeneParts(DMR.diffmeth.fem.gr.rename.hypo, ref.anno)
DMS.fem.anno <- annotateWithGeneParts(DMS.diffmeth.fem.gr.rename, ref.anno)
DMS.fem.anno.hyper <- annotateWithGeneParts(DMS.diffmeth.fem.gr.rename.hyper, ref.anno)
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
DMR.all.ann.stats.perc <- genomation::getTargetAnnotationStats(DMR.all.anno, percentage = TRUE, precedence = TRUE)
DMR.all.ann.stats.perc 
DMS.all.ann.stats.perc <- genomation::getTargetAnnotationStats(DMS.all.anno, percentage = TRUE, precedence = TRUE)
DMS.all.ann.stats.perc
CpG.all.ann.stats.perc <- genomation::getTargetAnnotationStats(CpG.all.anno, percentage = TRUE, precedence = TRUE)
CpG.all.ann.stats.perc
regions.all.ann.stats.perc <- genomation::getTargetAnnotationStats(regions.all.anno, percentage = TRUE, precedence = TRUE)
regions.all.ann.stats.perc

DMR.all.ann.stats.num <- genomation::getTargetAnnotationStats(DMR.all.anno, percentage = FALSE, precedence = TRUE)
DMR.all.ann.stats.num
DMS.all.ann.stats.num <- genomation::getTargetAnnotationStats(DMS.all.anno, percentage = FALSE, precedence = TRUE)
DMS.all.ann.stats.num
CpG.all.ann.stats.num <- genomation::getTargetAnnotationStats(CpG.all.anno, percentage = FALSE, precedence = TRUE)
CpG.all.ann.stats.num
regions.all.ann.stats.num <- genomation::getTargetAnnotationStats(regions.all.anno, percentage = FALSE, precedence = TRUE)
regions.all.ann.stats.num

DMR.fem.ann.stats.perc <- genomation::getTargetAnnotationStats(DMR.fem.anno, percentage = TRUE, precedence = TRUE)
DMR.fem.ann.stats.perc 
DMS.fem.ann.stats.perc <- genomation::getTargetAnnotationStats(DMS.fem.anno, percentage = TRUE, precedence = TRUE)
DMS.fem.ann.stats.perc
DMR.fem.ann.stats.perc.hyper <- genomation::getTargetAnnotationStats(DMR.fem.anno.hyper, percentage = TRUE, precedence = TRUE)
DMR.fem.ann.stats.perc.hyper 
DMS.fem.ann.stats.perc.hyper <- genomation::getTargetAnnotationStats(DMS.fem.anno.hyper, percentage = TRUE, precedence = TRUE)
DMS.fem.ann.stats.perc.hyper
DMR.fem.ann.stats.perc.hypo <- genomation::getTargetAnnotationStats(DMR.fem.anno.hypo, percentage = TRUE, precedence = TRUE)
DMR.fem.ann.stats.perc.hypo 
DMS.fem.ann.stats.perc.hypo <- genomation::getTargetAnnotationStats(DMS.fem.anno.hypo, percentage = TRUE, precedence = TRUE)
DMS.fem.ann.stats.perc.hypo
CpG.fem.ann.stats.perc <- genomation::getTargetAnnotationStats(CpG.fem.anno, percentage = TRUE, precedence = TRUE)
CpG.fem.ann.stats.perc
regions.fem.ann.stats.perc <- genomation::getTargetAnnotationStats(regions.fem.anno, percentage = TRUE, precedence = TRUE)
regions.fem.ann.stats.perc

DMR.fem.ann.stats.num <- genomation::getTargetAnnotationStats(DMR.fem.anno, percentage = FALSE, precedence = TRUE)
DMR.fem.ann.stats.num
DMS.fem.ann.stats.num <- genomation::getTargetAnnotationStats(DMS.fem.anno, percentage = FALSE, precedence = TRUE)
DMS.fem.ann.stats.num
DMR.fem.ann.stats.num.hyper <- genomation::getTargetAnnotationStats(DMR.fem.anno.hyper, percentage = FALSE, precedence = TRUE)
DMR.fem.ann.stats.num.hyper
DMS.fem.ann.stats.num.hyper <- genomation::getTargetAnnotationStats(DMS.fem.anno.hyper, percentage = FALSE, precedence = TRUE)
DMS.fem.ann.stats.num.hyper
DMR.fem.ann.stats.num.hypo <- genomation::getTargetAnnotationStats(DMR.fem.anno.hypo, percentage = FALSE, precedence = TRUE)
DMR.fem.ann.stats.num.hypo
DMS.fem.ann.stats.num.hypo <- genomation::getTargetAnnotationStats(DMS.fem.anno.hypo, percentage = FALSE, precedence = TRUE)
DMS.fem.ann.stats.num.hypo
CpG.fem.ann.stats.num <- genomation::getTargetAnnotationStats(CpG.fem.anno, percentage = FALSE, precedence = TRUE)
CpG.fem.ann.stats.num
regions.fem.ann.stats.num <- genomation::getTargetAnnotationStats(regions.fem.anno, percentage = FALSE, precedence = TRUE)
regions.fem.ann.stats.num

DMR.mal.ann.stats.perc <- genomation::getTargetAnnotationStats(DMR.mal.anno, percentage = TRUE, precedence = TRUE)
DMR.mal.ann.stats.perc 
DMS.mal.ann.stats.perc <- genomation::getTargetAnnotationStats(DMS.mal.anno, percentage = TRUE, precedence = TRUE)
DMS.mal.ann.stats.perc
DMR.mal.ann.stats.perc.hyper <- genomation::getTargetAnnotationStats(DMR.mal.anno.hyper, percentage = TRUE, precedence = TRUE)
DMR.mal.ann.stats.perc.hyper 
DMS.mal.ann.stats.perc.hyper <- genomation::getTargetAnnotationStats(DMS.mal.anno.hyper, percentage = TRUE, precedence = TRUE)
DMS.mal.ann.stats.perc.hyper
DMR.mal.ann.stats.perc.hypo <- genomation::getTargetAnnotationStats(DMR.mal.anno.hypo, percentage = TRUE, precedence = TRUE)
DMR.mal.ann.stats.perc.hypo 
DMS.mal.ann.stats.perc.hypo <- genomation::getTargetAnnotationStats(DMS.mal.anno.hypo, percentage = TRUE, precedence = TRUE)
DMS.mal.ann.stats.perc.hypo
CpG.mal.ann.stats.perc <- genomation::getTargetAnnotationStats(CpG.mal.anno, percentage = TRUE, precedence = TRUE)
CpG.mal.ann.stats.perc
regions.mal.ann.stats.perc <- genomation::getTargetAnnotationStats(regions.mal.anno, percentage = TRUE, precedence = TRUE)
regions.mal.ann.stats.perc

DMR.mal.ann.stats.num <- genomation::getTargetAnnotationStats(DMR.mal.anno, percentage = FALSE, precedence = TRUE)
DMR.mal.ann.stats.num
DMS.mal.ann.stats.num <- genomation::getTargetAnnotationStats(DMS.mal.anno, percentage = FALSE, precedence = TRUE)
DMS.mal.ann.stats.num
DMR.mal.ann.stats.num.hyper <- genomation::getTargetAnnotationStats(DMR.mal.anno.hyper, percentage = FALSE, precedence = TRUE)
DMR.mal.ann.stats.num.hyper
DMS.mal.ann.stats.num.hyper <- genomation::getTargetAnnotationStats(DMS.mal.anno.hyper, percentage = FALSE, precedence = TRUE)
DMS.mal.ann.stats.num.hyper
DMR.mal.ann.stats.num.hypo <- genomation::getTargetAnnotationStats(DMR.mal.anno.hypo, percentage = FALSE, precedence = TRUE)
DMR.mal.ann.stats.num.hypo
DMS.mal.ann.stats.num.hypo <- genomation::getTargetAnnotationStats(DMS.mal.anno.hypo, percentage = FALSE, precedence = TRUE)
DMS.mal.ann.stats.num.hypo
CpG.mal.ann.stats.num <- genomation::getTargetAnnotationStats(CpG.mal.anno, percentage = FALSE, precedence = TRUE)
CpG.mal.ann.stats.num
regions.mal.ann.stats.num <- genomation::getTargetAnnotationStats(regions.mal.anno, percentage = FALSE, precedence = TRUE)
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
saveRDS(DMR.fem.ann.stats.perc.hyper, file = "./anno_res/DMR_annStats_perc_fem_hyper.RDS")
saveRDS(DMS.fem.ann.stats.perc.hyper, file = "./anno_res/DMS_annStats_perc_fem_hyper.RDS")
saveRDS(DMR.fem.ann.stats.perc.hypo, file = "./anno_res/DMR_annStats_perc_fem_hypo.RDS")
saveRDS(DMS.fem.ann.stats.perc.hypo, file = "./anno_res/DMS_annStats_perc_fem_hypo.RDS")
saveRDS(CpG.fem.ann.stats.perc, file = "./anno_res/CpG_annStats_perc_fem.RDS")
saveRDS(regions.fem.ann.stats.perc, file = "./anno_res/regions_annStats_perc_fem.RDS")

saveRDS(DMR.mal.ann.stats.perc, file = "./anno_res/DMR_annStats_perc_mal.RDS")
saveRDS(DMS.mal.ann.stats.perc, file = "./anno_res/DMS_annStats_perc_mal.RDS")
saveRDS(DMR.mal.ann.stats.perc.hyper, file = "./anno_res/DMR_annStats_perc_mal_hyper.RDS")
saveRDS(DMS.mal.ann.stats.perc.hyper, file = "./anno_res/DMS_annStats_perc_mal_hyper.RDS")
saveRDS(DMR.mal.ann.stats.perc.hypo, file = "./anno_res/DMR_annStats_perc_mal_hypo.RDS")
saveRDS(DMS.mal.ann.stats.perc.hypo, file = "./anno_res/DMS_annStats_perc_mal_hypo.RDS")
saveRDS(CpG.mal.ann.stats.perc, file = "./anno_res/CpG_annStats_perc_mal.RDS")
saveRDS(regions.mal.ann.stats.perc, file = "./anno_res/regions_annStats_perc_mal.RDS")

saveRDS(DMR.all.ann.stats.num, file = "./anno_res/DMR_annStats_num_all.RDS")
saveRDS(DMS.all.ann.stats.num, file = "./anno_res/DMS_annStats_num_all.RDS")
saveRDS(CpG.all.ann.stats.num, file = "./CpG_annStats_num_all.RDS")
saveRDS(regions.all.ann.stats.num, file = "./anno_res/regions_annStats_num_all.RDS")

saveRDS(DMR.fem.ann.stats.num, file = "./anno_res/DMR_annStats_num_fem.RDS")
saveRDS(DMS.fem.ann.stats.num, file = "./anno_res/DMS_annStats_num_fem.RDS")
saveRDS(DMR.fem.ann.stats.num.hyper, file = "./anno_res/DMR_annStats_num_fem_hyper.RDS")
saveRDS(DMS.fem.ann.stats.num.hyper, file = "./anno_res/DMS_annStats_num_fem_hyper.RDS")
saveRDS(DMR.fem.ann.stats.num.hypo, file = "./anno_res/DMR_annStats_num_fem_hypo.RDS")
saveRDS(DMS.fem.ann.stats.num.hypo, file = "./anno_res/DMS_annStats_num_fem_hypo.RDS")
saveRDS(CpG.fem.ann.stats.num, file = "./anno_res/CpG_annStats_num_fem.RDS")
saveRDS(regions.fem.ann.stats.num, file = "./anno_res/regions_annStats_num_fem.RDS")

saveRDS(DMR.mal.ann.stats.num, file = "./anno_res/DMR_annStats_num_mal.RDS")
saveRDS(DMS.mal.ann.stats.num, file = "./anno_res/DMS_annStats_num_mal.RDS")
saveRDS(DMR.mal.ann.stats.num.hyper, file = "./anno_res/DMR_annStats_num_mal_hyper.RDS")
saveRDS(DMS.mal.ann.stats.num.hyper, file = "./anno_res/DMS_annStats_num_mal_hyper.RDS")
saveRDS(DMR.mal.ann.stats.num.hypo, file = "./anno_res/DMR_annStats_num_mal_hypo.RDS")
saveRDS(DMS.mal.ann.stats.num.hypo, file = "./anno_res/DMS_annStats_num_mal_hypo.RDS")
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




DMS_anno_stat_all <- readRDS("./dev_exp/data/methylkit_res/anno_res/DMS_annStats_perc_all.RDS")
DMS_anno_stat_fem <- readRDS("./dev_exp/data/methylkit_res/anno_res/DMS_annStats_perc_fem.RDS")
DMS_anno_stat_mal <- readRDS("./dev_exp/data/methylkit_res/anno_res/DMS_annStats_perc_mal.RDS")

DMR_anno_stat_all <- readRDS("./dev_exp/data/methylkit_res/anno_res/DMR_annStats_perc_all.RDS")
DMR_anno_stat_fem <- readRDS("./dev_exp/data/methylkit_res/anno_res/DMR_annStats_perc_fem.RDS")
DMR_anno_stat_fem <- readRDS("./dev_exp/data/methylkit_res/anno_res/DMR_annStats_perc_fem_hyper.RDS")
DMR_anno_stat_fem <- readRDS("./dev_exp/data/methylkit_res/anno_res/DMR_annStats_perc_fem_hypo.RDS")

DMR_anno_stat_mal <- readRDS("./dev_exp/data/methylkit_res/anno_res/DMR_annStats_perc_mal.RDS")

CpG_anno_stat_all <- readRDS("./dev_exp/data/methylkit_res/anno_res/CpG_annStats_perc_all.RDS")
CpG_anno_stat_fem <- readRDS("./dev_exp/data/methylkit_res/anno_res/CpG_annStats_perc_fem.RDS")
CpG_anno_stat_mal <- readRDS("./dev_exp/data/methylkit_res/anno_res/CpG_annStats_perc_mal.RDS")

region_anno_stat_all <- readRDS("./dev_exp/data/methylkit_res/anno_res/regions_annStats_perc_all.RDS")
regin_anno_stat_fem <- readRDS("./dev_exp/data/methylkit_res/anno_res/regions_annStats_perc_fem.RDS")
region_anno_stat_mal <- readRDS("./dev_exp/data/methylkit_res/anno_res/regions_annStats_perc_mal.RDS")
