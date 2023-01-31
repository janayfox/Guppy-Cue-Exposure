#############################################
### Goal: Get data from methylkit files
### Author: Janay Fox
### R script
#############################################

## Set up ##
#install packages 
#install.packages("S4Vectors")
#install.packages("IRanges")
#install.packages("GenomicRanges")
#install.packages("methylKit")

#load packages
library("IRanges", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")
library("GenomicRanges", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")
library("methylKit", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")
library("methylKit", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")

setwd("/scratch/janayfox/guppyWGBS/")
load(file="./backupRData/05_methylCallingDMS-backup.RData")

## Get data ## 
DMS.meth.2L.10X.data <- getData(DMS.meth.2L.10X)
DMS.meth.10L.10X.data <- getData(DMS.meth.10L.10X)
DMS.meth.2L.5X.data <- getData(DMS.meth.2L.5X)
DMS.meth.10L.5X.data <- getData(DMS.meth.10L.5X)

DMS.myDiff.2L.10X.data <- getData(DMS.myDiff.2L.10X)
DMS.myDiff.10L.10X.data <- getData(DMS.myDiff.10L.10X)
DMS.myDiff.2L.5X.data <- getData(DMS.myDiff.2L.5X)
DMS.myDiff.10L.5X.data <- getData(DMS.myDiff.10L.5X)

DMS.diffMeth.2L.10X.data <- getData(DMS.diffMeth.2L.10X)
DMS.diffMeth.10L.10X.data <- getData(DMS.diffMeth.10L.10X)
DMS.diffMeth.2L.5X.data <- getData(DMS.diffMeth.2L.5X)
DMS.diffMeth.10L.5X.data <- getData(DMS.diffMeth.10L.5X)

DMR.meth.2L10X.data <- getData(DMR.meth.2L10X)
DMR.meth.10L10X.data <- getData(DMR.meth.10L10X)
DMR.meth.2L5X.data <- getData(DMR.meth.2L5X)
DMR.meth.10L5X.data <- getData(DMR.meth.10L5X)

DMR.myDiff.2L.10X.data <- getData(DMR.myDiff.2L.10X)
DMR.myDiff.10L.10X.data <- getData(DMR.myDiff.10L.10X)
DMR.myDiff.2L.5X.data <- getData(DMR.myDiff.2L.5X)
DMR.myDiff.10L.5X.data <- getData(DMR.myDiff.10L.5X)

DMR.diffMeth.2L.10X.data <- getData(DMR.diffMeth.2L.10X)
DMR.diffMeth.10L.10X.data <- getData(DMR.diffMeth.10L.10X)
DMR.diffMeth.2L.5X.data <- getData(DMR.diffMeth.2L.5X)
DMR.diffMeth.10L.5X.data <- getData(DMR.diffMeth.10L.5X)

## Save RDS's ## 
saveRDS(DMS.meth.2L.10X.data, file ="./backupRData/DMS_meth_2l10X.rds")
saveRDS(DMS.meth.10L.10X.data, file ="./backupRData/DMS_meth_10l10X.rds")
saveRDS(DMS.meth.2L.5X.data, file ="./backupRData/DMS_meth_2l5X.rds")
saveRDS(DMS.meth.10L.5X.data, file ="./backupRData/DMS_meth_10l5X.rds")

saveRDS(DMS.myDiff.2L.10X.data, file ="./backupRData/DMS_myDiff_2l10X.rds")
saveRDS(DMS.myDiff.10L.10X.data, file ="./backupRData/DMS_myDiff_10l10X.rds")
saveRDS(DMS.myDiff.2L.5X.data, file ="./backupRData/DMS_myDiff_2l5X.rds")
saveRDS(DMS.myDiff.10L.5X.data, file ="./backupRData/DMS_myDiff_10l5X.rds")

saveRDS(DMS.diffMeth.2L.10X.data, file ="./backupRData/DMS_diffMeth_2l10X.rds")
saveRDS(DMS.diffMeth.10L.10X.data, file ="./backupRData/DMS_diffMeth_10l10X.rds")
saveRDS(DMS.diffMeth.2L.5X.data, file ="./backupRData/DMS_diffMeth_2l5X.rds")
saveRDS(DMS.diffMeth.10L.5X.data, file ="./backupRData/DMS_diffMeth_10l5X.rds")

saveRDS(DMR.meth.2L.10X.data, file ="./backupRData/DMR_meth_2l10X.rds")
saveRDS(DMR.meth.10L.10X.data, file ="./backupRData/DMR_meth_10l10X.rds")
saveRDS(DMR.meth.2L.5X.data, file ="./backupRData/DMR_meth_2l5X.rds")
saveRDS(DMR.meth.10L.5X.data, file ="./backupRData/DMR_meth_10l5X.rds")

saveRDS(DMR.myDiff.2L.10X.data, file ="./backupRData/DMR_myDiff_2l10X.rds")
saveRDS(DMR.myDiff.10L.10X.data, file ="./backupRData/DMR_myDiff_10l10X.rds")
saveRDS(DMR.myDiff.2L.5X.data, file ="./backupRData/DMR_myDiff_2l5X.rds")
saveRDS(DMR.myDiff.10L.5X.data, file ="./backupRData/DMR_myDiff_10l5X.rds")

saveRDS(DMR.diffMeth.2L.10X.data, file ="./backupRData/DMR_diffMeth_2l10X.rds")
saveRDS(DMR.diffMeth.10L.10X.data, file ="./backupRData/DMR_diffMeth_10l10X.rds")
saveRDS(DMR.diffMeth.2L.5X.data, file ="./backupRData/DMR_diffMeth_2l5X.rds")
saveRDS(DMR.diffMeth.10L.5X.data, file ="./backupRData/DMR_diffMeth_10l5X.rds")

## Save workspace image for later loading ##
# save.image(file = ".RData")
# save.image(file = "./backupRData/07_getData-backup.RData")
# 
# q(save="yes")