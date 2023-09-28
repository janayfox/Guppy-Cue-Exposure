#load packages 
library("S4Vectors", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")
library("IRanges", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")
library("GenomicRanges", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")
library("methylKit", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")
library("data.table", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")

setwd("/scratch/janayfox/guppyWGBS/")

##run with no over dispersion 
DMSmeth5X_21L <- readRDS("./newMethylKitRes/DMSmeth5X_21L.RDS")
DMSmeth5X_25L <- readRDS("./newMethylKitRes/DMSmeth5X_25L.RDS")

#no overdispersion correction
nod.DMSmyDiff5X_21L <- calculateDiffMeth(DMSmeth5X_21L, mc.cores=2, covariates=covariates, test="Chisq", save.db = TRUE, suffix = "myDiff")
nod.DMSmyDiff5X_25L <- calculateDiffMeth(DMSmeth5X_25L, mc.cores=2, covariates=covariates, test="Chisq", save.db = TRUE, suffix = "myDiff")

# Call significant methylation
nod.DMSdiffMeth5X_21L <- getMethylDiff(nod.DMSmyDiff5X_21L, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "diffMeth")
nod.DMSdiffMeth5X_25L <- getMethylDiff(nod.DMSmyDiff5X_25L, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "diffMeth")

# Check number of significant DMS
nod.DMSdiffMeth5X_21L
nod.DMSdiffMeth5X_25L

# Get meth per chromosome
nod.DMSdiffMethChr5X_21L <- diffMethPerChr(nod.DMSmyDiff5X_21L, plot=FALSE, qvalue.cutoff=0.0125, meth.cutoff=15, save.db = TRUE, suffix = "chr")
nod.DMSdiffMethChr5X_21L
nod.DMSdiffMethChr5X_25L <- diffMethPerChr(nod.DMSmyDiff5X_25L, plot=FALSE, qvalue.cutoff=0.0125, meth.cutoff=15, save.db = TRUE, suffix = "chr")
nod.DMSdiffMethChr5X_25L

saveRDS(nod.DMSmyDiff5X_21L, file = "./nod_DMSmyDiff5X_21L.RDS")
saveRDS(nod.DMSmyDiff5X_25L, file = "./nod_DMSmyDiff5X_25L.RDS")
saveRDS(nod.DMSdiffMeth5X_21L, file = "./nod_DMSdiffmeth5X_21L.RDS")
saveRDS(nod.DMSdiffMeth5X_25L, file = "./nod_DMSdiffMeth5X_25L.RDS")
saveRDS(nod.DMSdiffMethChr5X_21L, file = "./nod_DMSdiffMethChr5x_21L.RDS")
saveRDS(nod.DMSdiffMethChr5X_25L, file = "./nod_DMSdiffMethChr5X_25L.RDS")

saveRDS(getData(nod.DMSdiffMeth5X_21L), file = "./nod_DMSdiffMeth5X_21Ldata.RDS")
saveRDS(getData(nod.DMSdiffMeth5X_25L), file = "./nod_DMSdiffMeth5X_25Ldata.RDS")