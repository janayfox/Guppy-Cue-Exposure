################################################################################
### Goal: Get methylation tables and clean data for cluster analysis
### Author: Janay Fox
### R script
###############################################################################

#load packages 
library("S4Vectors", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")
library("IRanges", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")
library("GenomicRanges", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")
library("methylKit", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")
library("vegan", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")
library("goeveg", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")
library("tibble", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")
library("tidyr", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")
library("data.table", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")
library("dplyr", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")

setwd("/scratch/janayfox/guppyWGBS/")

#read in meth files
DMSmeth5X_25L <- readRDS("./newMethylKitRes/DMSmeth5X_25L.RDS")
DMSmeth5X_21L <- readRDS("./newMethylKitRes/DMSmeth5X_21L.RDS")
DMSdiffMeth5X_21L <- readRDS("./newMethylKitRes/nod_DMSdiffMeth5X_21Ldata.RDS")
DMSdiffMeth5X_25L <- readRDS("./newMethylKitRes/nod_DMSdiffMeth5X_25Ldata.RDS")

#read in metadata file
metaData <- read.csv("./clean_developmental_metadata.csv")
  
## Get percent methylation and clean data for cluster analysis##
#calculate percent methylation
percDMSmeth5X_21L <- percMethylation(DMSmeth5X_21L, rowids = TRUE)
percDMSmeth5X_25L <- percMethylation(DMSmeth5X_25L, rowids = TRUE)

#convert to dataframe
percDMSmeth5X_21L <- as.data.frame(percDMSmeth5X_21L)
percDMSmeth5X_25L <- as.data.frame(percDMSmeth5X_25L)

# #remove methylation sites with low variation (sd < 0.3)
# percDMSmeth5X_21L$sd <- apply(percDMSmeth5X_21L, 1, sd, na.rm = TRUE) #calculate sd
# percDMSmeth5X_25L$sd <- apply(percDMSmeth5X_25L, 1, sd, na.rm = TRUE)
# 
# percDMSmeth5X_21L <- subset(percDMSmeth5X_21L, sd >= 30) #subset sd >= 30 (would be equal to 0.3 percentage) 
# percDMSmeth5X_25L <- subset(percDMSmeth5X_25L, sd >= 30) #subset sd >= 30 (would be equal to 0.3 percentage)
# 
# percDMSmeth5X_21L <- subset(percDMSmeth5X_21L, select = -c(sd)) #remove sd column 
# percDMSmeth5X_25L <- subset(percDMSmeth5X_25L, select = -c(sd)) 

#transpose 
t.percDMSmeth5X_21L <- transpose(percDMSmeth5X_21L)
rownames(t.percDMSmeth5X_21L) <- colnames(percDMSmeth5X_21L)
colnames(t.percDMSmeth5X_21L) <- rownames(percDMSmeth5X_21L)

t.percDMSmeth5X_25L <- transpose(percDMSmeth5X_25L)
rownames(t.percDMSmeth5X_25L) <- colnames(percDMSmeth5X_25L)
colnames(t.percDMSmeth5X_25L) <- rownames(percDMSmeth5X_25L)

#format chromosome position to be the same in DMS as CpGs
DMSdiffMeth5X_21L <- unite(DMSdiffMeth5X_21L, col = "DMS", c("chr", "start", "end"), sep = ".")
DMS <- DMSdiffMeth5X_21L$DMS #extract DMS as list 

DMSdiffMeth5X_25L <- unite(DMSdiffMeth5X_25L, col = "DMS", c("chr", "start", "end"), sep = ".")
DMS <- DMSdiffMeth5X_25L$DMS #extract DMS as list 

#select only DMS from CpG
DMS_21L.dat <- select(t.percDMSmeth5X_21L, DMS)
DMS_25L.dat <- select(t.percDMSmeth5X_25L, DMS)

#convert row names to column for fish ID
t.percDMSmeth5X_21L<- tibble::rownames_to_column(t.percDMSmeth5X_21L, var = "fish_ID")
t.percDMSmeth5X_25L <- tibble::rownames_to_column(t.percDMSmeth5X_25L, var = "fish_ID")

DMS_21L.dat <- tibble::rownames_to_column(DMS_21L.dat, var = "fish_ID")
DMS_25L.dat <- tibble::rownames_to_column(DMS_25L.dat, var = "fish_ID")

#check dataframe 
print(t.percDMSmeth5X_21L[1:10,1:10])
print(t.percDMSmeth5X_25L[1:10,1:10])

print(DMS_21L.dat[1:40,1:10])
print(DMS_25L.dat[1:10,1:10])

#remove columns with NAs
t.percDMSmeth5X_21L <- t.percDMSmeth5X_21L[ , colSums(is.na(t.percDMSmeth5X_21L))==0]
t.percDMSmeth5X_25L <- t.percDMSmeth5X_25L[ , colSums(is.na(t.percDMSmeth5X_25L))==0]

DMS_21L.dat <- DMS_21L.dat[ , colSums(is.na(DMS_21L.dat))==0]
DMS_25L.dat <- DMS_25L.dat[ , colSums(is.na(DMS_25L.dat))==0]

#check number of columns 
ncol.21L <- ncol(t.percDMSmeth5X_21L)
ncol.21L
ncol.25L <- ncol(t.percDMSmeth5X_25L)
ncol.25L

ncol.DMS.21L <- ncol(DMS_21L.dat)
ncol.DMS.21L
ncol.DMS.25L <- ncol(DMS_25L.dat)
ncol.DMS.25L

#combine with metadata file 
all.data_21L <- merge(t.percDMSmeth5X_21L, metaData, by = "fish_ID")
all.data_25L <- merge(t.percDMSmeth5X_25L, metaData, by = "fish_ID")

DMS.all.data_21L <- merge(DMS_21L.dat, metaData, by = "fish_ID")
DMS.all.data_25L <- merge(DMS_25L.dat, metaData, by = "fish_ID")

#get number of columns 
ncol.21L.alldat <- ncol(all.data_21L)
ncol.21L.alldat
ncol.25L.alldat <- ncol(all.data_25L)
ncol.25L.alldat

ncol.21L.alldatDMS <- ncol(DMS.all.data_21L)
ncol.21L.alldatDMS
ncol.25L.alldatDMS <- ncol(DMS.all.data_25L)
ncol.25L.alldatDMS

#check dataframe 
print(all.data_21L[1:10,(ncol.21L.alldat-10):ncol.21L.alldat])
print(all.data_25L[1:10,(ncol.25L.alldat-10):ncol.25L.alldat])

## Save RDS files ##
saveRDS(percDMSmeth5X_21L, "./percDMSmeth5X_21L.RDS")
saveRDS(percDMSmeth5X_25L, "./percDMSmeth5X_25L.RDS")

saveRDS(all.data_21L, "./methAllData21L.RDS")
saveRDS(all.data_25L, "./methAllData25L.RDS")

saveRDS(DMS.all.data_21L, "./DMSmethAllData21L.RDS")
saveRDS(DMS.all.data_25L, "./DMSmethAllData25L.RDS")

