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

setwd("/scratch/janayfox/guppyWGBS_final/methylKit/st/perc20_od")

## Load in data and reformat ##
ref.anno <- readTranscriptFeatures("/scratch/janayfox/guppyWGBS_final/annotation/guppy.sorted.bed", unique.prom = FALSE, up.flank = 1500, down.flank = 500)

DMR.diffmeth.05h <- readRDS("./05h/DMR_res_05h/DMRdiffmeth_05h_5X.RDS")
DMR.diffmeth.1h <- readRDS("./1h/DMR_res_1h/DMRdiffmeth_1h_5X.RDS")
DMR.diffmeth.4h <- readRDS("./4h/DMR_res_4h/DMRdiffmeth_4h_5X.RDS")
DMR.diffmeth.24h <- readRDS("./24h/DMR_res_24h/DMRdiffmeth_24h_5X.RDS")
DMR.diffmeth.72h <- readRDS("./72h/DMR_res_72h/DMRdiffmeth_72h_5X.RDS")
DMR.diffmeth.all <- readRDS("./all/DMR_res_all/DMRdiffmeth_all_5X.RDS")

DMR.diffmeth.05h.fem <- readRDS("./05h/DMR_res_05h/DMRdiffmeth_05h_fem_5X.RDS")
DMR.diffmeth.1h.fem <- readRDS("./1h/DMR_res_1h/DMRdiffmeth_1h_fem_5X.RDS")
DMR.diffmeth.4h.fem <- readRDS("./4h/DMR_res_4h/DMRdiffmeth_4h_fem_5X.RDS")
DMR.diffmeth.24h.fem <- readRDS("./24h/DMR_res_24h/DMRdiffmeth_24h_fem_5X.RDS")
DMR.diffmeth.72h.fem <- readRDS("./72h/DMR_res_72h/DMRdiffmeth_72h_fem_5X.RDS")
DMR.diffmeth.all.fem <- readRDS("./all/DMR_res_all/DMRdiffmeth_all_fem_5X.RDS")

DMR.diffmeth.05h.mal <- readRDS("./05h/DMR_res_05h/DMRdiffmeth_05h_mal_5X.RDS")
DMR.diffmeth.1h.mal <- readRDS("./1h/DMR_res_1h/DMRdiffmeth_1h_mal_5X.RDS")
DMR.diffmeth.4h.mal <- readRDS("./4h/DMR_res_4h/DMRdiffmeth_4h_mal_5X.RDS")
DMR.diffmeth.24h.mal <- readRDS("./24h/DMR_res_24h/DMRdiffmeth_24h_mal_5X.RDS")
DMR.diffmeth.72h.mal <- readRDS("./72h/DMR_res_72h/DMRdiffmeth_72h_mal_5X.RDS")
DMR.diffmeth.all.mal <- readRDS("./all/DMR_res_all/DMRdiffmeth_all_mal_5X.RDS")

DMS.diffmeth.05h <- readRDS("./05h/DMS_res_05h/DMSdiffmeth_05h_5X.RDS")
DMS.diffmeth.1h <- readRDS("./1h/DMS_res_1h/DMSdiffmeth_1h_5X.RDS")
DMS.diffmeth.4h <- readRDS("./4h/DMS_res_4h/DMSdiffmeth_4h_5X.RDS")
DMS.diffmeth.24h <- readRDS("./24h/DMS_res_24h/DMSdiffmeth_24h_5X.RDS")
DMS.diffmeth.72h <- readRDS("./72h/DMS_res_72h/DMSdiffmeth_72h_5X.RDS")
DMS.diffmeth.all <- readRDS("./all/DMS_res_all/DMSdiffmeth_all_5X.RDS")

DMS.diffmeth.05h.fem <- readRDS("./05h/DMS_res_05h/DMSdiffmeth_05h_fem_5X.RDS")
DMS.diffmeth.1h.fem <- readRDS("./1h/DMS_res_1h/DMSdiffmeth_1h_fem_5X.RDS")
DMS.diffmeth.4h.fem <- readRDS("./4h/DMS_res_4h/DMSdiffmeth_4h_fem_5X.RDS")
DMS.diffmeth.24h.fem <- readRDS("./24h/DMS_res_24h/DMSdiffmeth_24h_fem_5X.RDS")
DMS.diffmeth.72h.fem <- readRDS("./72h/DMS_res_72h/DMSdiffmeth_72h_fem_5X.RDS")
DMS.diffmeth.all.fem <- readRDS("./all/DMS_res_all/DMSdiffmeth_all_fem_5X.RDS")

DMS.diffmeth.05h.mal <- readRDS("./05h/DMS_res_05h/DMSdiffmeth_05h_mal_5X.RDS")
DMS.diffmeth.1h.mal <- readRDS("./1h/DMS_res_1h/DMSdiffmeth_1h_mal_5X.RDS")
DMS.diffmeth.4h.mal <- readRDS("./4h/DMS_res_4h/DMSdiffmeth_4h_mal_5X.RDS")
DMS.diffmeth.24h.mal <- readRDS("./24h/DMS_res_24h/DMSdiffmeth_24h_mal_5X.RDS")
DMS.diffmeth.72h.mal <- readRDS("./72h/DMS_res_72h/DMSdiffmeth_72h_mal_5X.RDS")
DMS.diffmeth.all.mal <- readRDS("./all/DMS_res_all/DMSdiffmeth_all_mal_5X.RDS")

CpG.05h <- readRDS("./05h/DMS_res_05h/DMSmydiff_05h_5X.RDS")
CpG.1h <- readRDS("./1h/DMS_res_1h/DMSmydiff_1h_5X.RDS")
CpG.4h <- readRDS("./4h/DMS_res_4h/DMSmydiff_4h_5X.RDS")
CpG.24h <- readRDS("./24h/DMS_res_24h/DMSmydiff_24h_5X.RDS")
CpG.72h <- readRDS("./72h/DMS_res_72h/DMSmydiff_72h_5X.RDS")
CpG.all <- readRDS("./all/DMS_res_all/DMSmydiff_all_5X.RDS")

CpG.05h.fem <- readRDS("./05h/DMS_res_05h/DMSmydiff_05h_fem_5X.RDS")
CpG.1h.fem <- readRDS("./1h/DMS_res_1h/DMSmydiff_1h_fem_5X.RDS")
CpG.4h.fem <- readRDS("./4h/DMS_res_4h/DMSmydiff_4h_fem_5X.RDS")
CpG.24h.fem <- readRDS("./24h/DMS_res_24h/DMSmydiff_24h_fem_5X.RDS")
CpG.72h.fem <- readRDS("./72h/DMS_res_72h/DMSmydiff_72h_fem_5X.RDS")
CpG.all.fem <- readRDS("./all/DMS_res_all/DMSmydiff_all_fem_5X.RDS")

CpG.05h.mal <- readRDS("./05h/DMS_res_05h/DMSmydiff_05h_mal_5X.RDS")
CpG.1h.mal <- readRDS("./1h/DMS_res_1h/DMSmydiff_1h_mal_5X.RDS")
CpG.4h.mal <- readRDS("./4h/DMS_res_4h/DMSmydiff_4h_mal_5X.RDS")
CpG.24h.mal <- readRDS("./24h/DMS_res_24h/DMSmydiff_24h_mal_5X.RDS")
CpG.72h.mal <- readRDS("./72h/DMS_res_72h/DMSmydiff_72h_mal_5X.RDS")
CpG.all.mal <- readRDS("./all/DMS_res_all/DMSmydiff_all_mal_5X.RDS")

regions.05h <- readRDS("./05h/DMR_res_05h/DMRmydiff_05h_5X.RDS")
regions.1h <- readRDS("./1h/DMR_res_1h/DMRmydiff_1h_5X.RDS")
regions.4h <- readRDS("./4h/DMR_res_4h/DMRmydiff_4h_5X.RDS")
regions.24h <- readRDS("./24h/DMR_res_24h/DMRmydiff_24h_5X.RDS")
regions.72h <- readRDS("./72h/DMR_res_72h/DMRmydiff_72h_5X.RDS")
regions.all <- readRDS("./all/DMR_res_all/DMRmydiff_all_5X.RDS")

regions.05h.fem <- readRDS("./05h/DMR_res_05h/DMRmydiff_05h_fem_5X.RDS")
regions.1h.fem <- readRDS("./1h/DMR_res_1h/DMRmydiff_1h_fem_5X.RDS")
regions.4h.fem <- readRDS("./4h/DMR_res_4h/DMRmydiff_4h_fem_5X.RDS")
regions.24h.fem <- readRDS("./24h/DMR_res_24h/DMRmydiff_24h_fem_5X.RDS")
regions.72h.fem <- readRDS("./72h/DMR_res_72h/DMRmydiff_72h_fem_5X.RDS")
regions.all.fem <- readRDS("./all/DMR_res_all/DMRmydiff_all_fem_5X.RDS")

regions.05h.mal <- readRDS("./05h/DMR_res_05h/DMRmydiff_05h_mal_5X.RDS")
regions.1h.mal <- readRDS("./1h/DMR_res_1h/DMRmydiff_1h_mal_5X.RDS")
regions.4h.mal <- readRDS("./4h/DMR_res_4h/DMRmydiff_4h_mal_5X.RDS")
regions.24h.mal <- readRDS("./24h/DMR_res_24h/DMRmydiff_24h_mal_5X.RDS")
regions.72h.mal <- readRDS("./72h/DMR_res_72h/DMRmydiff_72h_mal_5X.RDS")
regions.all.mal <- readRDS("./all/DMR_res_all/DMRmydiff_all_mal_5X.RDS")

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

regions.05h.gr.rename <- renameChr.noXY(regions.05h)
regions.1h.gr.rename <- renameChr.noXY(regions.1h)
regions.4h.gr.rename <- renameChr.noXY(regions.4h)
regions.24h.gr.rename <- renameChr.noXY(regions.24h)
regions.72h.gr.rename <- renameChr.noXY(regions.72h)
regions.all.gr.rename <- renameChr.noXY(regions.all)

regions.05h.fem.gr.rename <- renameChr(regions.05h.fem)
regions.1h.fem.gr.rename <- renameChr(regions.1h.fem)
regions.4h.fem.gr.rename <- renameChr(regions.4h.fem)
regions.24h.fem.gr.rename <- renameChr(regions.24h.fem)
regions.72h.fem.gr.rename <- renameChr(regions.72h.fem)
regions.all.fem.gr.rename <- renameChr(regions.all.fem)

regions.05h.mal.gr.rename <- renameChr(regions.05h.mal)
regions.1h.mal.gr.rename <- renameChr(regions.1h.mal)
regions.4h.mal.gr.rename <- renameChr(regions.4h.mal)
regions.24h.mal.gr.rename <- renameChr(regions.24h.mal)
regions.72h.mal.gr.rename <- renameChr(regions.72h.mal)
regions.all.mal.gr.rename <- renameChr(regions.all.mal)

# get lists of hypo and hyper methylation
DMR.diffmeth.05h.gr.rename.hyper <- subset(DMR.diffmeth.05h.gr.rename, meth.diff > 0)
DMR.diffmeth.1h.gr.rename.hyper <- subset(DMR.diffmeth.1h.gr.rename, meth.diff > 0)
DMR.diffmeth.4h.gr.rename.hyper <- subset(DMR.diffmeth.4h.gr.rename, meth.diff > 0)
DMR.diffmeth.24h.gr.rename.hyper <- subset(DMR.diffmeth.24h.gr.rename, meth.diff > 0)
DMR.diffmeth.72h.gr.rename.hyper <- subset(DMR.diffmeth.72h.gr.rename, meth.diff > 0)
DMR.diffmeth.all.gr.rename.hyper <- subset(DMR.diffmeth.all.gr.rename, meth.diff > 0)

DMR.diffmeth.05h.fem.gr.rename.hyper <- subset(DMR.diffmeth.05h.fem.gr.rename, meth.diff > 0)
DMR.diffmeth.1h.fem.gr.rename.hyper <- subset(DMR.diffmeth.1h.fem.gr.rename, meth.diff > 0)
DMR.diffmeth.4h.fem.gr.rename.hyper <- subset(DMR.diffmeth.4h.fem.gr.rename, meth.diff > 0)
DMR.diffmeth.24h.fem.gr.rename.hyper <- subset(DMR.diffmeth.24h.fem.gr.rename, meth.diff > 0)
DMR.diffmeth.72h.fem.gr.rename.hyper <- subset(DMR.diffmeth.72h.fem.gr.rename, meth.diff > 0)
DMR.diffmeth.all.fem.gr.rename.hyper <- subset(DMR.diffmeth.all.fem.gr.rename, meth.diff > 0)

DMR.diffmeth.05h.mal.gr.rename.hyper <- subset(DMR.diffmeth.05h.mal.gr.rename, meth.diff > 0)
DMR.diffmeth.1h.mal.gr.rename.hyper <- subset(DMR.diffmeth.1h.mal.gr.rename, meth.diff > 0)
DMR.diffmeth.4h.mal.gr.rename.hyper <- subset(DMR.diffmeth.4h.mal.gr.rename, meth.diff > 0)
DMR.diffmeth.24h.mal.gr.rename.hyper <- subset(DMR.diffmeth.24h.mal.gr.rename, meth.diff > 0)
DMR.diffmeth.72h.mal.gr.rename.hyper <- subset(DMR.diffmeth.72h.mal.gr.rename, meth.diff > 0)
DMR.diffmeth.all.mal.gr.rename.hyper <- subset(DMR.diffmeth.all.mal.gr.rename, meth.diff > 0)

DMS.diffmeth.05h.gr.rename.hyper <- subset(DMS.diffmeth.05h.gr.rename, meth.diff > 0)
DMS.diffmeth.1h.gr.rename.hyper <- subset(DMS.diffmeth.1h.gr.rename, meth.diff > 0)
DMS.diffmeth.4h.gr.rename.hyper <- subset(DMS.diffmeth.4h.gr.rename, meth.diff > 0)
DMS.diffmeth.24h.gr.rename.hyper <- subset(DMS.diffmeth.24h.gr.rename, meth.diff > 0)
DMS.diffmeth.72h.gr.rename.hyper <- subset(DMS.diffmeth.72h.gr.rename, meth.diff > 0)
DMS.diffmeth.all.gr.rename.hyper <- subset(DMS.diffmeth.all.gr.rename, meth.diff > 0)

DMS.diffmeth.05h.fem.gr.rename.hyper <- subset(DMS.diffmeth.05h.fem.gr.rename, meth.diff > 0)
DMS.diffmeth.1h.fem.gr.rename.hyper <- subset(DMS.diffmeth.1h.fem.gr.rename, meth.diff > 0)
DMS.diffmeth.4h.fem.gr.rename.hyper <- subset(DMS.diffmeth.4h.fem.gr.rename, meth.diff > 0)
DMS.diffmeth.24h.fem.gr.rename.hyper <- subset(DMS.diffmeth.24h.fem.gr.rename, meth.diff > 0)
DMS.diffmeth.72h.fem.gr.rename.hyper <- subset(DMS.diffmeth.72h.fem.gr.rename, meth.diff > 0)
DMS.diffmeth.all.fem.gr.rename.hyper <- subset(DMS.diffmeth.all.fem.gr.rename, meth.diff > 0)

DMS.diffmeth.05h.mal.gr.rename.hyper <- subset(DMS.diffmeth.05h.mal.gr.rename, meth.diff > 0)
DMS.diffmeth.1h.mal.gr.rename.hyper <- subset(DMS.diffmeth.1h.mal.gr.rename, meth.diff > 0)
DMS.diffmeth.4h.mal.gr.rename.hyper <- subset(DMS.diffmeth.4h.mal.gr.rename, meth.diff > 0)
DMS.diffmeth.24h.mal.gr.rename.hyper <- subset(DMS.diffmeth.24h.mal.gr.rename, meth.diff > 0)
DMS.diffmeth.72h.mal.gr.rename.hyper <- subset(DMS.diffmeth.72h.mal.gr.rename, meth.diff > 0)
DMS.diffmeth.all.mal.gr.rename.hyper <- subset(DMS.diffmeth.all.mal.gr.rename, meth.diff > 0)

# get hypo methyaltion 
DMR.diffmeth.05h.gr.rename.hypo <- subset(DMR.diffmeth.05h.gr.rename, meth.diff < 0)
DMR.diffmeth.1h.gr.rename.hypo <- subset(DMR.diffmeth.1h.gr.rename, meth.diff < 0)
DMR.diffmeth.4h.gr.rename.hypo <- subset(DMR.diffmeth.4h.gr.rename, meth.diff < 0)
DMR.diffmeth.24h.gr.rename.hypo <- subset(DMR.diffmeth.24h.gr.rename, meth.diff < 0)
DMR.diffmeth.72h.gr.rename.hypo <- subset(DMR.diffmeth.72h.gr.rename, meth.diff < 0)
DMR.diffmeth.all.gr.rename.hypo <- subset(DMR.diffmeth.all.gr.rename, meth.diff < 0)

DMR.diffmeth.05h.fem.gr.rename.hypo <- subset(DMR.diffmeth.05h.fem.gr.rename, meth.diff < 0)
DMR.diffmeth.1h.fem.gr.rename.hypo <- subset(DMR.diffmeth.1h.fem.gr.rename, meth.diff < 0)
DMR.diffmeth.4h.fem.gr.rename.hypo <- subset(DMR.diffmeth.4h.fem.gr.rename, meth.diff < 0)
DMR.diffmeth.24h.fem.gr.rename.hypo <- subset(DMR.diffmeth.24h.fem.gr.rename, meth.diff < 0)
DMR.diffmeth.72h.fem.gr.rename.hypo <- subset(DMR.diffmeth.72h.fem.gr.rename, meth.diff < 0)
DMR.diffmeth.all.fem.gr.rename.hypo <- subset(DMR.diffmeth.all.fem.gr.rename, meth.diff < 0)

DMR.diffmeth.05h.mal.gr.rename.hypo <- subset(DMR.diffmeth.05h.mal.gr.rename, meth.diff < 0)
DMR.diffmeth.1h.mal.gr.rename.hypo <- subset(DMR.diffmeth.1h.mal.gr.rename, meth.diff < 0)
DMR.diffmeth.4h.mal.gr.rename.hypo <- subset(DMR.diffmeth.4h.mal.gr.rename, meth.diff < 0)
DMR.diffmeth.24h.mal.gr.rename.hypo <- subset(DMR.diffmeth.24h.mal.gr.rename, meth.diff < 0)
DMR.diffmeth.72h.mal.gr.rename.hypo <- subset(DMR.diffmeth.72h.mal.gr.rename, meth.diff < 0)
DMR.diffmeth.all.mal.gr.rename.hypo <- subset(DMR.diffmeth.all.mal.gr.rename, meth.diff < 0)

DMS.diffmeth.05h.gr.rename.hypo <- subset(DMS.diffmeth.05h.gr.rename, meth.diff < 0)
DMS.diffmeth.1h.gr.rename.hypo <- subset(DMS.diffmeth.1h.gr.rename, meth.diff < 0)
DMS.diffmeth.4h.gr.rename.hypo <- subset(DMS.diffmeth.4h.gr.rename, meth.diff < 0)
DMS.diffmeth.24h.gr.rename.hypo <- subset(DMS.diffmeth.24h.gr.rename, meth.diff < 0)
DMS.diffmeth.72h.gr.rename.hypo <- subset(DMS.diffmeth.72h.gr.rename, meth.diff < 0)
DMS.diffmeth.all.gr.rename.hypo <- subset(DMS.diffmeth.all.gr.rename, meth.diff < 0)

DMS.diffmeth.05h.fem.gr.rename.hypo <- subset(DMS.diffmeth.05h.fem.gr.rename, meth.diff < 0)
DMS.diffmeth.1h.fem.gr.rename.hypo <- subset(DMS.diffmeth.1h.fem.gr.rename, meth.diff < 0)
DMS.diffmeth.4h.fem.gr.rename.hypo <- subset(DMS.diffmeth.4h.fem.gr.rename, meth.diff < 0)
DMS.diffmeth.24h.fem.gr.rename.hypo <- subset(DMS.diffmeth.24h.fem.gr.rename, meth.diff < 0)
DMS.diffmeth.72h.fem.gr.rename.hypo <- subset(DMS.diffmeth.72h.fem.gr.rename, meth.diff < 0)
DMS.diffmeth.all.fem.gr.rename.hypo <- subset(DMS.diffmeth.all.fem.gr.rename, meth.diff < 0)

DMS.diffmeth.05h.mal.gr.rename.hypo <- subset(DMS.diffmeth.05h.mal.gr.rename, meth.diff < 0)
DMS.diffmeth.1h.mal.gr.rename.hypo <- subset(DMS.diffmeth.1h.mal.gr.rename, meth.diff < 0)
DMS.diffmeth.4h.mal.gr.rename.hypo <- subset(DMS.diffmeth.4h.mal.gr.rename, meth.diff < 0)
DMS.diffmeth.24h.mal.gr.rename.hypo <- subset(DMS.diffmeth.24h.mal.gr.rename, meth.diff < 0)
DMS.diffmeth.72h.mal.gr.rename.hypo <- subset(DMS.diffmeth.72h.mal.gr.rename, meth.diff < 0)
DMS.diffmeth.all.mal.gr.rename.hypo <- subset(DMS.diffmeth.all.mal.gr.rename, meth.diff < 0)

## Annotate ##
anno.func <- function(DMR, DMR.hyper, DMR.hypo, DMS, DMS.hyper, DMS.hypo, CpG, regions,
                      DMR_anno_name, DMR_hyper_anno_name, DMR_hypo_anno_name, 
                      DMS_anno_name, DMS_hyper_anno_name, DMS_hypo_anno_name,
                      CpG_anno_name, regions_anno_name,
                      DMR_perc_name, DMS_perc_name, CpG_perc_name, regions_perc_name,
                      DMR_num_name, DMS_num_name, CpG_num_name, regions_num_name,
                      DMR_tss_name, DMR_hyper_tss_name, DMR_hypo_tss_name,
                      DMS_tss_name, DMS_hyper_tss_name, DMS_hypo_tss_name) {

  #annotate with gene parts 
  DMR.anno <- annotateWithGeneParts(DMR, ref.anno)
  DMR.hyper.anno <- annotateWithGeneParts(DMR.hyper, ref.anno)
  DMR.hypo.anno <- annotateWithGeneParts(DMR.hypo, ref.anno)

  DMS.anno <- annotateWithGeneParts(DMS, ref.anno)
  DMS.hyper.anno <- annotateWithGeneParts(DMS.hyper, ref.anno)
  DMS.hypo.anno <- annotateWithGeneParts(DMS.hypo, ref.anno)  
  
  CpG.anno <- annotateWithGeneParts(CpG, ref.anno)
  regions.anno <- annotateWithGeneParts(regions, ref.anno)
  
  #get percentage of DMRs that overlap with different features 
  DMR.ann.stats.perc <- getTargetAnnotationStats(DMR.anno, percentage = TRUE, precedence = TRUE)
  DMS.ann.stats.perc <- getTargetAnnotationStats(DMS.anno, percentage = TRUE, precedence = TRUE)
  CpG.ann.stats.perc <- getTargetAnnotationStats(CpG.anno, percentage = TRUE, precedence = TRUE)
  regions.ann.stats.perc <- getTargetAnnotationStats(regions.anno, percentage = TRUE, precedence = TRUE)

  DMR.ann.stats.num <- getTargetAnnotationStats(DMR.anno, percentage = FALSE, precedence = TRUE)
  DMS.ann.stats.num <- getTargetAnnotationStats(DMS.anno, percentage = FALSE, precedence = TRUE)
  CpG.ann.stats.num <- getTargetAnnotationStats(CpG.anno, percentage = FALSE, precedence = TRUE)
  regions.ann.stats.num <- getTargetAnnotationStats(regions.anno, percentage = FALSE, precedence = TRUE)
 
  #get nearest TSS for DMS and DMRs
  DMR.tss <- getAssociationWithTSS(DMR.anno)
  DMR.hyper.tss <- getAssociationWithTSS(DMR.hyper.anno)
  DMR.hypo.tss <- getAssociationWithTSS(DMR.hypo.anno)

  DMS.tss <- getAssociationWithTSS(DMS.anno)
  DMS.hyper.tss <- getAssociationWithTSS(DMS.hyper.anno)
  DMS.hypo.tss <- getAssociationWithTSS(DMS.hypo.anno)

  ## Save data ## 
  saveRDS(DMR.anno, file = DMR_anno_name)
  saveRDS(DMR.hyper.anno, file = DMR_hyper_anno_name)
  saveRDS(DMR.hypo.anno, file = DMR_hypo_anno_name)
  
  saveRDS(DMS.anno, file = DMS_anno_name)
  saveRDS(DMS.hyper.anno, file = DMS_hyper_anno_name)
  saveRDS(DMS.hypo.anno, file = DMS_hypo_anno_name)

  saveRDS(CpG.anno, file = CpG_anno_name)
  saveRDS(regions.anno, file = regions_anno_name)
  
  saveRDS(DMR.ann.stats.perc, file = DMR_perc_name)
  saveRDS(DMS.ann.stats.perc, file = DMS_perc_name)
  saveRDS(CpG.ann.stats.perc, file = CpG_perc_name)
  saveRDS(regions.ann.stats.perc, file = regions_perc_name)

  saveRDS(DMR.ann.stats.num, file = DMR_num_name)
  saveRDS(DMS.ann.stats.num, file = DMS_num_name)
  saveRDS(CpG.ann.stats.num, file = CpG_num_name)
  saveRDS(regions.ann.stats.num, file = regions_num_name)

  saveRDS(DMR.tss, file = DMR_tss_name)
  saveRDS(DMR.hyper.tss, file = DMR_hyper_tss_name)
  saveRDS(DMR.hypo.tss, file = DMR_hypo_tss_name)

  saveRDS(DMS.tss, file = DMS_tss_name)
  saveRDS(DMS.hyper.tss, file = DMS_hyper_tss_name)
  saveRDS(DMS.hypo.tss, file = DMS_hypo_tss_name)
}

#run on all sites 
anno.func(DMR.diffmeth.05h.gr.rename, DMR.diffmeth.05h.gr.rename.hyper, DMR.diffmeth.05h.gr.rename.hypo, 
          DMS.diffmeth.05h.gr.rename, DMS.diffmeth.05h.gr.rename.hyper, DMS.diffmeth.05h.gr.rename.hypo,
          CpG.05h.gr.rename, regions.05h.gr.rename,
          "./05h/anno_res_05h/DMR_anno_05h.RDS", "./05h/anno_res_05h/DMR_anno_05h_hyper.RDS", "./05h/anno_res_05h/DMR_anno_05h_hypo.RDS",  
          "./05h/anno_res_05h/DMS_anno_05h.RDS", "./05h/anno_res_05h/DMS_anno_05h_hyper.RDS", "./05h/anno_res_05h/DMS_anno_05h_hypo.RDS",  
          "./05h/anno_res_05h/CpG_anno_05h.RDS", "./05h/anno_res_05h/regions_anno_05h.RDS",
          "./05h/anno_res_05h/DMR_annStats_perc_05h.RDS", "./05h/anno_res_05h/DMS_annStats_perc_05h.RDS", 
          "./05h/anno_res_05h/CpG_annStats_perc_05h.RDS", "./05h/anno_res_05h/regions_annStats_perc_05h.RDS",
          "./05h/anno_res_05h/DMR_annStats_num_05h.RDS", "./05h/anno_res_05h/DMS_annStats_num_05h.RDS", 
          "./05h/anno_res_05h/CpG_annStats_num_05h.RDS", "./05h/anno_res_05h/regions_annStats_num_05h.RDS",
          "./05h/anno_res_05h/DMR_TSS_05h.RDS", "./05h/anno_res_05h/DMR_TSS_05h_hyper.RDS", "./05h/anno_res_05h/DMR_TSS_05h_hypo.RDS",
          "./05h/anno_res_05h/DMS_TSS_05h.RDS", "./05h/anno_res_05h/DMS_TSS_05h_hyper.RDS", "./05h/anno_res_05h/DMS_TSS_05h_hypo.RDS")

anno.func(DMR.diffmeth.05h.fem.gr.rename, DMR.diffmeth.05h.fem.gr.rename.hyper, DMR.diffmeth.05h.fem.gr.rename.hypo, 
          DMS.diffmeth.05h.fem.gr.rename, DMS.diffmeth.05h.fem.gr.rename.hyper, DMS.diffmeth.05h.fem.gr.rename.hypo,
          CpG.05h.fem.gr.rename, regions.05h.fem.gr.rename,
          "./05h/anno_res_05h/DMR_anno_05h_fem.RDS", "./05h/anno_res_05h/DMR_anno_05h_fem_hyper.RDS", "./05h/anno_res_05h/DMR_anno_05h_fem_hypo.RDS", 
          "./05h/anno_res_05h/DMS_anno_05h_fem.RDS", "./05h/anno_res_05h/DMS_anno_05h_fem_hyper.RDS", "./05h/anno_res_05h/DMS_anno_05h_fem_hypo.RDS",
          "./05h/anno_res_05h/CpG_anno_05h_fem.RDS", "./05h/anno_res_05h/regions_anno_05h_fem.RDS",
          "./05h/anno_res_05h/DMR_annStats_perc_05h_fem.RDS", "./05h/anno_res_05h/DMS_annStats_perc_05h_fem.RDS", 
          "./05h/anno_res_05h/CpG_annStats_perc_05h_fem.RDS", "./05h/anno_res_05h/regions_annStats_perc_05h_fem.RDS", 
          "./05h/anno_res_05h/DMR_annStats_num_05h_fem.RDS", "./05h/anno_res_05h/DMS_annStats_num_05h_fem.RDS", 
          "./05h/anno_res_05h/CpG_annStats_num_05h_fem.RDS", "./05h/anno_res_05h/regions_annStats_num_05h_fem.RDS",
          "./05h/anno_res_05h/DMR_TSS_05h_fem.RDS", "./05h/anno_res_05h/DMR_TSS_05h_fem_hyper.RDS", "./05h/anno_res_05h/DMR_TSS_05h_fem_hypo.RDS", 
          "./05h/anno_res_05h/DMS_TSS_05h_fem.RDS", "./05h/anno_res_05h/DMS_TSS_05h_fem_hyper.RDS", "./05h/anno_res_05h/DMS_TSS_05h_fem_hypo.RDS")

anno.func(DMR.diffmeth.05h.mal.gr.rename, DMR.diffmeth.05h.mal.gr.rename.hyper, DMR.diffmeth.05h.mal.gr.rename.hypo, 
          DMS.diffmeth.05h.mal.gr.rename, DMS.diffmeth.05h.mal.gr.rename.hyper, DMS.diffmeth.05h.mal.gr.rename.hypo, 
          CpG.05h.mal.gr.rename, regions.05h.mal.gr.rename,
          "./05h/anno_res_05h/DMR_anno_05h_mal.RDS", "./05h/anno_res_05h/DMR_anno_05h_mal_hyper.RDS", "./05h/anno_res_05h/DMR_anno_05h_mal_hypo.RDS", 
          "./05h/anno_res_05h/DMS_anno_05h_mal.RDS", "./05h/anno_res_05h/DMS_anno_05h_mal_hyper.RDS", "./05h/anno_res_05h/DMS_anno_05h_mal_hypo.RDS", 
          "./05h/anno_res_05h/CpG_anno_05h_mal.RDS", "./05h/anno_res_05h/regions_anno_05h_mal.RDS",
          "./05h/anno_res_05h/DMR_annStats_perc_05h_mal.RDS", "./05h/anno_res_05h/DMS_annStats_perc_05h_mal.RDS", 
          "./05h/anno_res_05h/CpG_annStats_perc_05h_mal.RDS", "./05h/anno_res_05h/regions_annStats_perc_05h_mal.RDS",
          "./05h/anno_res_05h/DMR_annStats_num_05h_mal.RDS", "./05h/anno_res_05h/DMS_annStats_num_05h_mal.RDS", 
          "./05h/anno_res_05h/CpG_annStats_num_05h_mal.RDS", "./05h/anno_res_05h/regions_annStats_num_05h_mal.RDS",
          "./05h/anno_res_05h/DMR_TSS_05h_mal.RDS", "./05h/anno_res_05h/DMR_TSS_05h_mal_hyper.RDS", "./05h/anno_res_05h/DMR_TSS_05h_mal_hypo.RDS", 
          "./05h/anno_res_05h/DMS_TSS_05h_mal.RDS", "./05h/anno_res_05h/DMS_TSS_05h_mal_hyper.RDS", "./05h/anno_res_05h/DMS_TSS_05h_mal_hypo.RDS")

anno.func(DMR.diffmeth.1h.gr.rename, DMR.diffmeth.1h.gr.rename.hyper, DMR.diffmeth.1h.gr.rename.hypo,
          DMS.diffmeth.1h.gr.rename, DMS.diffmeth.1h.gr.rename.hyper, DMS.diffmeth.1h.gr.rename.hypo,
          CpG.1h.gr.rename, regions.1h.gr.rename,
          "./1h/anno_res_1h/DMR_anno_1h.RDS", "./1h/anno_res_1h/DMR_anno_1h_hyper.RDS", "./1h/anno_res_1h/DMR_anno_1h_hypo.RDS",
          "./1h/anno_res_1h/DMS_anno_1h.RDS", "./1h/anno_res_1h/DMS_anno_1h_hyper.RDS", "./1h/anno_res_1h/DMS_anno_1h_hypo.RDS",
          "./1h/anno_res_1h/CpG_anno_1h.RDS", "./1h/anno_res_1h/regions_anno_1h.RDS",
          "./1h/anno_res_1h/DMR_annStats_perc_1h.RDS", "./1h/anno_res_1h/DMS_annStats_perc_1h.RDS",
          "./1h/anno_res_1h/CpG_annStats_perc_1h.RDS", "./1h/anno_res_1h/regions_annStats_perc_1h.RDS",
          "./1h/anno_res_1h/DMR_annStats_num_1h.RDS", "./1h/anno_res_1h/DMS_annStats_num_1h.RDS",
          "./1h/anno_res_1h/CpG_annStats_num_1h.RDS", "./1h/anno_res_1h/regions_annStats_num_1h.RDS",
          "./1h/anno_res_1h/DMR_TSS_1h.RDS", "./1h/anno_res_1h/DMR_TSS_1h_hyper.RDS", "./1h/anno_res_1h/DMR_TSS_1h_hypo.RDS",
          "./1h/anno_res_1h/DMS_TSS_1h.RDS", "./1h/anno_res_1h/DMS_TSS_1h_hyper.RDS", "./1h/anno_res_1h/DMS_TSS_1h_hypo.RDS")

anno.func(DMR.diffmeth.1h.fem.gr.rename, DMR.diffmeth.1h.fem.gr.rename.hyper, DMR.diffmeth.1h.fem.gr.rename.hypo,
          DMS.diffmeth.1h.fem.gr.rename, DMS.diffmeth.1h.fem.gr.rename.hyper, DMS.diffmeth.1h.fem.gr.rename.hypo,
          CpG.1h.fem.gr.rename, regions.1h.fem.gr.rename,
          "./1h/anno_res_1h/DMR_anno_1h_fem.RDS", "./1h/anno_res_1h/DMR_anno_1h_fem_hyper.RDS", "./1h/anno_res_1h/DMR_anno_1h_fem_hypo.RDS",
          "./1h/anno_res_1h/DMS_anno_1h_fem.RDS", "./1h/anno_res_1h/DMS_anno_1h_fem_hyper.RDS", "./1h/anno_res_1h/DMS_anno_1h_fem_hypo.RDS",
          "./1h/anno_res_1h/CpG_anno_1h_fem.RDS", "./1h/anno_res_1h/regions_anno_1h_fem.RDS",
          "./1h/anno_res_1h/DMR_annStats_perc_1h_fem.RDS", "./1h/anno_res_1h/DMS_annStats_perc_1h_fem.RDS",
          "./1h/anno_res_1h/CpG_annStats_perc_1h_fem.RDS", "./1h/anno_res_1h/regions_annStats_perc_1h_fem.RDS",
          "./1h/anno_res_1h/DMR_annStats_num_1h_fem.RDS", "./1h/anno_res_1h/DMS_annStats_num_1h_fem.RDS",
          "./1h/anno_res_1h/CpG_annStats_num_1h_fem.RDS", "./1h/anno_res_1h/regions_annStats_num_1h_fem.RDS",
          "./1h/anno_res_1h/DMR_TSS_1h_fem.RDS", "./1h/anno_res_1h/DMR_TSS_1h_fem_hyper.RDS", "./1h/anno_res_1h/DMR_TSS_1h_fem_hypo.RDS",
          "./1h/anno_res_1h/DMS_TSS_1h_fem.RDS", "./1h/anno_res_1h/DMS_TSS_1h_fem_hyper.RDS", "./1h/anno_res_1h/DMS_TSS_1h_fem_hypo.RDS")

anno.func(DMR.diffmeth.1h.mal.gr.rename, DMR.diffmeth.1h.mal.gr.rename.hyper, DMR.diffmeth.1h.mal.gr.rename.hypo,
          DMS.diffmeth.1h.mal.gr.rename, DMS.diffmeth.1h.mal.gr.rename.hyper, DMS.diffmeth.1h.mal.gr.rename.hypo,
          CpG.1h.mal.gr.rename, regions.1h.mal.gr.rename,
          "./1h/anno_res_1h/DMR_anno_1h_mal.RDS", "./1h/anno_res_1h/DMR_anno_1h_mal_hyper.RDS", "./1h/anno_res_1h/DMR_anno_1h_mal_hypo.RDS",
          "./1h/anno_res_1h/DMS_anno_1h_mal.RDS", "./1h/anno_res_1h/DMS_anno_1h_mal_hyper.RDS", "./1h/anno_res_1h/DMS_anno_1h_mal_hypo.RDS",
          "./1h/anno_res_1h/CpG_anno_1h_mal.RDS", "./1h/anno_res_1h/regions_anno_1h_mal.RDS",
          "./1h/anno_res_1h/DMR_annStats_perc_1h_mal.RDS", "./1h/anno_res_1h/DMS_annStats_perc_1h_mal.RDS",
          "./1h/anno_res_1h/CpG_annStats_perc_1h_mal.RDS", "./1h/anno_res_1h/regions_annStats_perc_1h_mal.RDS",
          "./1h/anno_res_1h/DMR_annStats_num_1h_mal.RDS", "./1h/anno_res_1h/DMS_annStats_num_1h_mal.RDS",
          "./1h/anno_res_1h/CpG_annStats_num_1h_mal.RDS", "./1h/anno_res_1h/regions_annStats_num_1h_mal.RDS",
          "./1h/anno_res_1h/DMR_TSS_1h_mal.RDS", "./1h/anno_res_1h/DMR_TSS_1h_mal_hyper.RDS", "./1h/anno_res_1h/DMR_TSS_1h_mal_hypo.RDS",
          "./1h/anno_res_1h/DMS_TSS_1h_mal.RDS", "./1h/anno_res_1h/DMS_TSS_1h_mal_hyper.RDS", "./1h/anno_res_1h/DMS_TSS_1h_mal_hypo.RDS")

anno.func(DMR.diffmeth.4h.gr.rename, DMR.diffmeth.4h.gr.rename.hyper, DMR.diffmeth.4h.gr.rename.hypo,
          DMS.diffmeth.4h.gr.rename, DMS.diffmeth.4h.gr.rename.hyper, DMS.diffmeth.4h.gr.rename.hypo,
          CpG.4h.gr.rename, regions.4h.gr.rename,
          "./4h/anno_res_4h/DMR_anno_4h.RDS", "./4h/anno_res_4h/DMR_anno_4h_hyper.RDS", "./4h/anno_res_4h/DMR_anno_4h_hypo.RDS",
          "./4h/anno_res_4h/DMS_anno_4h.RDS", "./4h/anno_res_4h/DMS_anno_4h_hyper.RDS", "./4h/anno_res_4h/DMS_anno_4h_hypo.RDS",
          "./4h/anno_res_4h/CpG_anno_4h.RDS", "./4h/anno_res_4h/regions_anno_4h.RDS",
          "./4h/anno_res_4h/DMR_annStats_perc_4h.RDS", "./4h/anno_res_4h/DMS_annStats_perc_4h.RDS",
          "./4h/anno_res_4h/CpG_annStats_perc_4h.RDS", "./4h/anno_res_4h/regions_annStats_perc_4h.RDS",
          "./4h/anno_res_4h/DMR_annStats_num_4h.RDS", "./4h/anno_res_4h/DMS_annStats_num_4h.RDS",
          "./4h/anno_res_4h/CpG_annStats_num_4h.RDS", "./4h/anno_res_4h/regions_annStats_num_4h.RDS",
          "./4h/anno_res_4h/DMR_TSS_4h.RDS", "./4h/anno_res_4h/DMR_TSS_4h_hyper.RDS", "./4h/anno_res_4h/DMR_TSS_4h_hypo.RDS",
          "./4h/anno_res_4h/DMS_TSS_4h.RDS", "./4h/anno_res_4h/DMS_TSS_4h_hyper.RDS", "./4h/anno_res_4h/DMS_TSS_4h_hypo.RDS")

anno.func(DMR.diffmeth.4h.fem.gr.rename, DMR.diffmeth.4h.fem.gr.rename.hyper, DMR.diffmeth.4h.fem.gr.rename.hypo,
          DMS.diffmeth.4h.fem.gr.rename, DMS.diffmeth.4h.fem.gr.rename.hyper, DMS.diffmeth.4h.fem.gr.rename.hypo,
          CpG.4h.fem.gr.rename, regions.4h.fem.gr.rename,
          "./4h/anno_res_4h/DMR_anno_4h_fem.RDS", "./4h/anno_res_4h/DMR_anno_4h_fem_hyper.RDS", "./4h/anno_res_4h/DMR_anno_4h_fem_hypo.RDS",
          "./4h/anno_res_4h/DMS_anno_4h_fem.RDS", "./4h/anno_res_4h/DMS_anno_4h_fem_hyper.RDS", "./4h/anno_res_4h/DMS_anno_4h_fem_hypo.RDS",
          "./4h/anno_res_4h/CpG_anno_4h_fem.RDS", "./4h/anno_res_4h/regions_anno_4h_fem.RDS",
          "./4h/anno_res_4h/DMR_annStats_perc_4h_fem.RDS", "./4h/anno_res_4h/DMS_annStats_perc_4h_fem.RDS",
          "./4h/anno_res_4h/CpG_annStats_perc_4h_fem.RDS", "./4h/anno_res_4h/regions_annStats_perc_4h_fem.RDS",
          "./4h/anno_res_4h/DMR_annStats_num_4h_fem.RDS", "./4h/anno_res_4h/DMS_annStats_num_4h_fem.RDS",
          "./4h/anno_res_4h/CpG_annStats_num_4h_fem.RDS", "./4h/anno_res_4h/regions_annStats_num_4h_fem.RDS",
          "./4h/anno_res_4h/DMR_TSS_4h_fem.RDS", "./4h/anno_res_4h/DMR_TSS_4h_fem_hyper.RDS", "./4h/anno_res_4h/DMR_TSS_4h_fem_hypo.RDS",
          "./4h/anno_res_4h/DMS_TSS_4h_fem.RDS", "./4h/anno_res_4h/DMS_TSS_4h_fem_hyper.RDS", "./4h/anno_res_4h/DMS_TSS_4h_fem_hypo.RDS")

anno.func(DMR.diffmeth.4h.mal.gr.rename, DMR.diffmeth.4h.mal.gr.rename.hyper, DMR.diffmeth.4h.mal.gr.rename.hypo,
          DMS.diffmeth.4h.mal.gr.rename, DMS.diffmeth.4h.mal.gr.rename.hyper, DMS.diffmeth.4h.mal.gr.rename.hypo,
          CpG.4h.mal.gr.rename, regions.4h.mal.gr.rename,
          "./4h/anno_res_4h/DMR_anno_4h_mal.RDS", "./4h/anno_res_4h/DMR_anno_4h_mal_hyper.RDS", "./4h/anno_res_4h/DMR_anno_4h_mal_hypo.RDS",
          "./4h/anno_res_4h/DMS_anno_4h_mal.RDS", "./4h/anno_res_4h/DMS_anno_4h_mal_hyper.RDS", "./4h/anno_res_4h/DMS_anno_4h_mal_hypo.RDS",
          "./4h/anno_res_4h/CpG_anno_4h_mal.RDS", "./4h/anno_res_4h/regions_anno_4h_mal.RDS",
          "./4h/anno_res_4h/DMR_annStats_perc_4h_mal.RDS", "./4h/anno_res_4h/DMS_annStats_perc_4h_mal.RDS",
          "./4h/anno_res_4h/CpG_annStats_perc_4h_mal.RDS", "./4h/anno_res_4h/regions_annStats_perc_4h_mal.RDS",
          "./4h/anno_res_4h/DMR_annStats_num_4h_mal.RDS", "./4h/anno_res_4h/DMS_annStats_num_4h_mal.RDS",
          "./4h/anno_res_4h/CpG_annStats_num_4h_mal.RDS", "./4h/anno_res_4h/regions_annStats_num_4h_mal.RDS",
          "./4h/anno_res_4h/DMR_TSS_4h_mal.RDS", "./4h/anno_res_4h/DMR_TSS_4h_mal_hyper.RDS", "./4h/anno_res_4h/DMR_TSS_4h_mal_hypo.RDS",
          "./4h/anno_res_4h/DMS_TSS_4h_mal.RDS", "./4h/anno_res_4h/DMS_TSS_4h_mal_hyper.RDS", "./4h/anno_res_4h/DMS_TSS_4h_mal_hypo.RDS")

anno.func(DMR.diffmeth.24h.gr.rename, DMR.diffmeth.24h.gr.rename.hyper, DMR.diffmeth.24h.gr.rename.hypo,
          DMS.diffmeth.24h.gr.rename, DMS.diffmeth.24h.gr.rename.hyper, DMS.diffmeth.24h.gr.rename.hypo,
          CpG.24h.gr.rename, regions.24h.gr.rename,
          "./24h/anno_res_24h/DMR_anno_24h.RDS", "./24h/anno_res_24h/DMR_anno_24h_hyper.RDS", "./24h/anno_res_24h/DMR_anno_24h_hypo.RDS",
          "./24h/anno_res_24h/DMS_anno_24h.RDS", "./24h/anno_res_24h/DMS_anno_24h_hyper.RDS", "./24h/anno_res_24h/DMS_anno_24h_hypo.RDS",
          "./24h/anno_res_24h/CpG_anno_24h.RDS", "./24h/anno_res_24h/regions_anno_24h.RDS",
          "./24h/anno_res_24h/DMR_annStats_perc_24h.RDS", "./24h/anno_res_24h/DMS_annStats_perc_24h.RDS",
          "./24h/anno_res_24h/CpG_annStats_perc_24h.RDS", "./24h/anno_res_24h/regions_annStats_perc_24h.RDS",
          "./24h/anno_res_24h/DMR_annStats_num_24h.RDS", "./24h/anno_res_24h/DMS_annStats_num_24h.RDS",
          "./24h/anno_res_24h/CpG_annStats_num_24h.RDS", "./24h/anno_res_24h/regions_annStats_num_24h.RDS",
          "./24h/anno_res_24h/DMR_TSS_24h.RDS", "./24h/anno_res_24h/DMR_TSS_24h_hyper.RDS", "./24h/anno_res_24h/DMR_TSS_24h_hypo.RDS",
          "./24h/anno_res_24h/DMS_TSS_24h.RDS", "./24h/anno_res_24h/DMS_TSS_24h_hyper.RDS", "./24h/anno_res_24h/DMS_TSS_24h_hypo.RDS")

anno.func(DMR.diffmeth.24h.fem.gr.rename, DMR.diffmeth.24h.fem.gr.rename.hyper, DMR.diffmeth.24h.fem.gr.rename.hypo,
          DMS.diffmeth.24h.fem.gr.rename, DMS.diffmeth.24h.fem.gr.rename.hyper, DMS.diffmeth.24h.fem.gr.rename.hypo,
          CpG.24h.fem.gr.rename, regions.24h.fem.gr.rename,
          "./24h/anno_res_24h/DMR_anno_24h_fem.RDS", "./24h/anno_res_24h/DMR_anno_24h_fem_hyper.RDS", "./24h/anno_res_24h/DMR_anno_24h_fem_hypo.RDS",
          "./24h/anno_res_24h/DMS_anno_24h_fem.RDS", "./24h/anno_res_24h/DMS_anno_24h_fem_hyper.RDS", "./24h/anno_res_24h/DMS_anno_24h_fem_hypo.RDS",
          "./24h/anno_res_24h/CpG_anno_24h_fem.RDS", "./24h/anno_res_24h/regions_anno_24h_fem.RDS",
          "./24h/anno_res_24h/DMR_annStats_perc_24h_fem.RDS", "./24h/anno_res_24h/DMS_annStats_perc_24h_fem.RDS",
          "./24h/anno_res_24h/CpG_annStats_perc_24h_fem.RDS", "./24h/anno_res_24h/regions_annStats_perc_24h_fem.RDS",
          "./24h/anno_res_24h/DMR_annStats_num_24h_fem.RDS", "./24h/anno_res_24h/DMS_annStats_num_24h_fem.RDS",
          "./24h/anno_res_24h/CpG_annStats_num_24h_fem.RDS", "./24h/anno_res_24h/regions_annStats_num_24h_fem.RDS",
          "./24h/anno_res_24h/DMR_TSS_24h_fem.RDS", "./24h/anno_res_24h/DMR_TSS_24h_fem_hyper.RDS", "./24h/anno_res_24h/DMR_TSS_24h_fem_hypo.RDS",
          "./24h/anno_res_24h/DMS_TSS_24h_fem.RDS", "./24h/anno_res_24h/DMS_TSS_24h_fem_hyper.RDS", "./24h/anno_res_24h/DMS_TSS_24h_fem_hypo.RDS")

anno.func(DMR.diffmeth.24h.mal.gr.rename, DMR.diffmeth.24h.mal.gr.rename.hyper, DMR.diffmeth.24h.mal.gr.rename.hypo,
          DMS.diffmeth.24h.mal.gr.rename, DMS.diffmeth.24h.mal.gr.rename.hyper, DMS.diffmeth.24h.mal.gr.rename.hypo,
          CpG.24h.mal.gr.rename, regions.24h.mal.gr.rename,
          "./24h/anno_res_24h/DMR_anno_24h_mal.RDS", "./24h/anno_res_24h/DMR_anno_24h_mal_hyper.RDS", "./24h/anno_res_24h/DMR_anno_24h_mal_hypo.RDS",
          "./24h/anno_res_24h/DMS_anno_24h_mal.RDS", "./24h/anno_res_24h/DMS_anno_24h_mal_hyper.RDS", "./24h/anno_res_24h/DMS_anno_24h_mal_hypo.RDS",
          "./24h/anno_res_24h/CpG_anno_24h_mal.RDS", "./24h/anno_res_24h/regions_anno_24h_mal.RDS",
          "./24h/anno_res_24h/DMR_annStats_perc_24h_mal.RDS", "./24h/anno_res_24h/DMS_annStats_perc_24h_mal.RDS",
          "./24h/anno_res_24h/CpG_annStats_perc_24h_mal.RDS", "./24h/anno_res_24h/regions_annStats_perc_24h_mal.RDS",
          "./24h/anno_res_24h/DMR_annStats_num_24h_mal.RDS", "./24h/anno_res_24h/DMS_annStats_num_24h_mal.RDS",
          "./24h/anno_res_24h/CpG_annStats_num_24h_mal.RDS", "./24h/anno_res_24h/regions_annStats_num_24h_mal.RDS",
          "./24h/anno_res_24h/DMR_TSS_24h_mal.RDS", "./24h/anno_res_24h/DMR_TSS_24h_mal_hyper.RDS", "./24h/anno_res_24h/DMR_TSS_24h_mal_hypo.RDS",
          "./24h/anno_res_24h/DMS_TSS_24h_mal.RDS", "./24h/anno_res_24h/DMS_TSS_24h_mal_hyper.RDS", "./24h/anno_res_24h/DMS_TSS_24h_mal_hypo.RDS")

anno.func(DMR.diffmeth.72h.gr.rename, DMR.diffmeth.72h.gr.rename.hyper, DMR.diffmeth.72h.gr.rename.hypo,
          DMS.diffmeth.72h.gr.rename, DMS.diffmeth.72h.gr.rename.hyper, DMS.diffmeth.72h.gr.rename.hypo,
          CpG.72h.gr.rename, regions.72h.gr.rename,
          "./72h/anno_res_72h/DMR_anno_72h.RDS", "./72h/anno_res_72h/DMR_anno_72h_hyper.RDS", "./72h/anno_res_72h/DMR_anno_72h_hypo.RDS",
          "./72h/anno_res_72h/DMS_anno_72h.RDS", "./72h/anno_res_72h/DMS_anno_72h_hyper.RDS", "./72h/anno_res_72h/DMS_anno_72h_hypo.RDS",
          "./72h/anno_res_72h/CpG_anno_72h.RDS", "./72h/anno_res_72h/regions_anno_72h.RDS",
          "./72h/anno_res_72h/DMR_annStats_perc_72h.RDS", "./72h/anno_res_72h/DMS_annStats_perc_72h.RDS",
          "./72h/anno_res_72h/CpG_annStats_perc_72h.RDS", "./72h/anno_res_72h/regions_annStats_perc_72h.RDS",
          "./72h/anno_res_72h/DMR_annStats_num_72h.RDS", "./72h/anno_res_72h/DMS_annStats_num_72h.RDS",
          "./72h/anno_res_72h/CpG_annStats_num_72h.RDS", "./72h/anno_res_72h/regions_annStats_num_72h.RDS",
          "./72h/anno_res_72h/DMR_TSS_72h.RDS", "./72h/anno_res_72h/DMR_TSS_72h_hyper.RDS", "./72h/anno_res_72h/DMR_TSS_72h_hypo.RDS",
          "./72h/anno_res_72h/DMS_TSS_72h.RDS", "./72h/anno_res_72h/DMS_TSS_72h_hyper.RDS", "./72h/anno_res_72h/DMS_TSS_72h_hypo.RDS")

anno.func(DMR.diffmeth.72h.fem.gr.rename, DMR.diffmeth.72h.fem.gr.rename.hyper, DMR.diffmeth.72h.fem.gr.rename.hypo,
          DMS.diffmeth.72h.fem.gr.rename, DMS.diffmeth.72h.fem.gr.rename.hyper, DMS.diffmeth.72h.fem.gr.rename.hypo,
          CpG.72h.fem.gr.rename, regions.72h.fem.gr.rename,
          "./72h/anno_res_72h/DMR_anno_72h_fem.RDS", "./72h/anno_res_72h/DMR_anno_72h_fem_hyper.RDS", "./72h/anno_res_72h/DMR_anno_72h_fem_hypo.RDS",
          "./72h/anno_res_72h/DMS_anno_72h_fem.RDS", "./72h/anno_res_72h/DMS_anno_72h_fem_hyper.RDS", "./72h/anno_res_72h/DMS_anno_72h_fem_hypo.RDS",
          "./72h/anno_res_72h/CpG_anno_72h_fem.RDS", "./72h/anno_res_72h/regions_anno_72h_fem.RDS",
          "./72h/anno_res_72h/DMR_annStats_perc_72h_fem.RDS", "./72h/anno_res_72h/DMS_annStats_perc_72h_fem.RDS",
          "./72h/anno_res_72h/CpG_annStats_perc_72h_fem.RDS", "./72h/anno_res_72h/regions_annStats_perc_72h_fem.RDS",
          "./72h/anno_res_72h/DMR_annStats_num_72h_fem.RDS", "./72h/anno_res_72h/DMS_annStats_num_72h_fem.RDS",
          "./72h/anno_res_72h/CpG_annStats_num_72h_fem.RDS", "./72h/anno_res_72h/regions_annStats_num_72h_fem.RDS",
          "./72h/anno_res_72h/DMR_TSS_72h_fem.RDS", "./72h/anno_res_72h/DMR_TSS_72h_fem_hyper.RDS", "./72h/anno_res_72h/DMR_TSS_72h_fem_hypo.RDS",
          "./72h/anno_res_72h/DMS_TSS_72h_fem.RDS", "./72h/anno_res_72h/DMS_TSS_72h_fem_hyper.RDS", "./72h/anno_res_72h/DMS_TSS_72h_fem_hypo.RDS")

anno.func(DMR.diffmeth.72h.mal.gr.rename, DMR.diffmeth.72h.mal.gr.rename.hyper, DMR.diffmeth.72h.mal.gr.rename.hypo,
          DMS.diffmeth.72h.mal.gr.rename, DMS.diffmeth.72h.mal.gr.rename.hyper, DMS.diffmeth.72h.mal.gr.rename.hypo,
          CpG.72h.mal.gr.rename, regions.72h.mal.gr.rename,
          "./72h/anno_res_72h/DMR_anno_72h_mal.RDS", "./72h/anno_res_72h/DMR_anno_72h_mal_hyper.RDS", "./72h/anno_res_72h/DMR_anno_72h_mal_hypo.RDS",
          "./72h/anno_res_72h/DMS_anno_72h_mal.RDS", "./72h/anno_res_72h/DMS_anno_72h_mal_hyper.RDS", "./72h/anno_res_72h/DMS_anno_72h_mal_hypo.RDS",
          "./72h/anno_res_72h/CpG_anno_72h_mal.RDS", "./72h/anno_res_72h/regions_anno_72h_mal.RDS",
          "./72h/anno_res_72h/DMR_annStats_perc_72h_mal.RDS", "./72h/anno_res_72h/DMS_annStats_perc_72h_mal.RDS",
          "./72h/anno_res_72h/CpG_annStats_perc_72h_mal.RDS", "./72h/anno_res_72h/regions_annStats_perc_72h_mal.RDS",
          "./72h/anno_res_72h/DMR_annStats_num_72h_mal.RDS", "./72h/anno_res_72h/DMS_annStats_num_72h_mal.RDS",
          "./72h/anno_res_72h/CpG_annStats_num_72h_mal.RDS", "./72h/anno_res_72h/regions_annStats_num_72h_mal.RDS",
          "./72h/anno_res_72h/DMR_TSS_72h_mal.RDS", "./72h/anno_res_72h/DMR_TSS_72h_mal_hyper.RDS", "./72h/anno_res_72h/DMR_TSS_72h_mal_hypo.RDS",
          "./72h/anno_res_72h/DMS_TSS_72h_mal.RDS", "./72h/anno_res_72h/DMS_TSS_72h_mal_hyper.RDS", "./72h/anno_res_72h/DMS_TSS_72h_mal_hypo.RDS")

anno.func(DMR.diffmeth.all.gr.rename, DMR.diffmeth.all.gr.rename.hyper, DMR.diffmeth.all.gr.rename.hypo,
          DMS.diffmeth.all.gr.rename, DMS.diffmeth.all.gr.rename.hyper, DMS.diffmeth.all.gr.rename.hypo,
          CpG.all.gr.rename, regions.all.gr.rename,
          "./all/anno_res_all/DMR_anno_all.RDS", "./all/anno_res_all/DMR_anno_all_hyper.RDS", "./all/anno_res_all/DMR_anno_all_hypo.RDS",
          "./all/anno_res_all/DMS_anno_all.RDS", "./all/anno_res_all/DMS_anno_all_hyper.RDS", "./all/anno_res_all/DMS_anno_all_hypo.RDS",
          "./all/anno_res_all/CpG_anno_all.RDS", "./all/anno_res_all/regions_anno_all.RDS",
          "./all/anno_res_all/DMR_annStats_perc_all.RDS", "./all/anno_res_all/DMS_annStats_perc_all.RDS",
          "./all/anno_res_all/CpG_annStats_perc_all.RDS", "./all/anno_res_all/region_annStats_perc_all.RDS",
          "./all/anno_res_all/DMR_annStats_num_all.RDS", "./all/anno_res_all/DMS_annStats_num_all.RDS",
          "./all/anno_res_all/CpG_annStats_num_all.RDS", "./all/anno_res_all/regions_annStats_num_all.RDS",
          "./all/anno_res_all/DMR_TSS_all.RDS", "./all/anno_res_all/DMR_TSS_all_hyper.RDS", "./all/anno_res_all/DMR_TSS_all_hypo.RDS",
          "./all/anno_res_all/DMS_TSS_all.RDS", "./all/anno_res_all/DMS_TSS_all_hyper.RDS", "./all/anno_res_all/DMS_TSS_all_hypo.RDS")

anno.func(DMR.diffmeth.all.fem.gr.rename, DMR.diffmeth.all.fem.gr.rename.hyper, DMR.diffmeth.all.fem.gr.rename.hypo,
          DMS.diffmeth.all.fem.gr.rename, DMS.diffmeth.all.fem.gr.rename.hyper, DMS.diffmeth.all.fem.gr.rename.hypo,
          CpG.all.fem.gr.rename, regions.all.fem.gr.rename,
          "./all/anno_res_all/DMR_anno_all_fem.RDS", "./all/anno_res_all/DMR_anno_all_fem_hyper.RDS", "./all/anno_res_all/DMR_anno_all_fem_hypo.RDS",
          "./all/anno_res_all/DMS_anno_all_fem.RDS", "./all/anno_res_all/DMS_anno_all_fem_hyper.RDS", "./all/anno_res_all/DMS_anno_all_fem_hypo.RDS",
          "./all/anno_res_all/CpG_anno_all_fem.RDS", "./all/anno_res_all/regions_anno_all_fem.RDS",
          "./all/anno_res_all/DMR_annStats_perc_all_fem.RDS", "./all/anno_res_all/DMS_annStats_perc_all_fem.RDS",
          "./all/anno_res_all/CpG_annStats_perc_all_fem.RDS", "./all/anno_res_all/regions_annStats_perc_all_fem.RDS",
          "./all/anno_res_all/DMR_annStats_num_all_fem.RDS", "./all/anno_res_all/DMS_annStats_num_all_fem.RDS",
          "./all/anno_res_all/CpG_annStats_num_all_fem.RDS", "./all/anno_res_all/regions_annStats_num_all_fem.RDS",
          "./all/anno_res_all/DMR_TSS_all_fem.RDS", "./all/anno_res_all/DMR_TSS_all_fem_hyper.RDS", "./all/anno_res_all/DMR_TSS_all_fem_hypo.RDS",
          "./all/anno_res_all/DMS_TSS_all_fem.RDS", "./all/anno_res_all/DMS_TSS_all_fem_hyper.RDS", "./all/anno_res_all/DMS_TSS_all_fem_hypo.RDS")

anno.func(DMR.diffmeth.all.mal.gr.rename, DMR.diffmeth.all.mal.gr.rename.hyper, DMR.diffmeth.all.mal.gr.rename.hypo,
          DMS.diffmeth.all.mal.gr.rename, DMS.diffmeth.all.mal.gr.rename.hyper, DMS.diffmeth.all.mal.gr.rename.hypo,
          CpG.all.mal.gr.rename, regions.all.mal.gr.rename,
          "./all/anno_res_all/DMR_anno_all_mal.RDS", "./all/anno_res_all/DMR_anno_all_mal_hyper.RDS", "./all/anno_res_all/DMR_anno_all_mal_hypo.RDS",
          "./all/anno_res_all/DMS_anno_all_mal.RDS", "./all/anno_res_all/DMS_anno_all_mal_hyper.RDS", "./all/anno_res_all/DMS_anno_all_mal_hypo.RDS",
          "./all/anno_res_all/CpG_anno_all_mal.RDS", "./all/anno_res_all/regions_anno_all_mal.RDS",
          "./all/anno_res_all/DMR_annStats_perc_all_mal.RDS", "./all/anno_res_all/DMS_annStats_perc_all_mal.RDS",
          "./all/anno_res_all/CpG_annStats_perc_all_mal.RDS", "./all/anno_res_all/regions_annStats_perc_all_mal.RDS",
          "./all/anno_res_all/DMR_annStats_num_all_mal.RDS", "./all/anno_res_all/DMS_annStats_num_all_mal.RDS",
          "./all/anno_res_all/CpG_annStats_num_all_mal.RDS", "./all/anno_res_all/regions_annStats_num_all_mal.RDS",
          "./all/anno_res_all/DMR_TSS_all_mal.RDS", "./all/anno_res_all/DMR_TSS_all_mal_hyper.RDS", "./all/anno_res_all/DMR_TSS_all_mal_hypo.RDS",
          "./all/anno_res_all/DMS_TSS_all_mal.RDS", "./all/anno_res_all/DMS_TSS_all_mal_hyper.RDS", "./all/anno_res_all/DMS_TSS_all_mal_hypo.RDS")
