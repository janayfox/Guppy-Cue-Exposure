##############################################################
### Goal: Get data for CpGs, DMSs, and DMRs
### Author: Janay Fox
### R script
#############################################################

## Set up ##
#install packages 
#install.packages("S4Vectors")
#install.packages("IRanges")
#install.packages("methylKit")

#load packages
library("S4Vectors", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")
library("IRanges", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")
library("methylKit", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")

setwd("/scratch/janayfox/guppyWGBS/")
load(file=".RData")

#get data 
data.CpG <- getData(meth)
data.DMS <- getData(DMS.diffMeth)
data.DMR <- getData(DMR.diffMeth)

#perc.meth <- percMethylation(meth)
#would this help with anything?? 

#save data 
saveRDS(data.CpG, "../guppyWGBS/data_CpG.rds")
saveRDS(data.DMS, "../guppyWGBS/data_DMS.rds")
saveRDS(data.DMR, "../guppyWGBS/data_DMR.rds")

## get CpG info - not using right now
# #find number of CpGs
# nrow(meth.2L.5X)
# nrow(meth.2L.10X)
# 
# #get CpG data
# data.CpG.2L.5X <- getData(meth.2L.5X)
# data.CpG.2L.10X <- getData(meth.2L.10X)
# 
# #save data
# saveRDS(data.CpG.2L.5X, "../guppyWGBS/data_CpG_2L_5X.rds")
# saveRDS(data.CpG.2L.10X, "../guppyWGBS/data_CpG_2L_10X.rds")
