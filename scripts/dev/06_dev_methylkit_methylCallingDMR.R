###########################################################################################
### Goal: Remove sex chr and unplaced scaffolds and methylation calling for DMRs 
### Author: Janay Fox
### R script
###########################################################################################

## Set up  ##
#install packages 
#install.packages("S4Vectors")
#install.packages("IRanges")
#install.packages("GenomicRanges")
#install.packages("methylKit")

#load packages
library("S4Vectors", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")
library("IRanges", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")
library("GenomicRanges", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")
library("methylKit", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")

setwd("/scratch/janayfox/guppyWGBS/")

#read in R objects 
myobj.3X <- readRDS("./myObj3X.RDS")

#read in tank covariate data
covariates <- data.frame(tank=c("AC2","AC2","AC2","AC3","AC3","AC3",
                             "AC3","AC4","AC4","AC4","AC4","AC4","AC4",
                             "AC5", "AC5","AC5","AC5","AC5","AC5","AC5",
                             "AC5","AC6","AC6","AC6","AC6","AC6", "AC6",
                             "AC7","AC7","AC7","AC7","AC7","AC7","AC7",
                             "AC7","AC7","C2","C2", "C2","C2","C3","C3","C3",
                             "C3","C3","C3","C4","C4","C4","C4",
                             "C4","C4","C4","C5","C5","C5","C5",
                             "C5","C5","C5","C5","C6","C6","C6",
                             "C6","C6","C6","C6","C7","C7","C7",
                             "C7","C7","C7","C7","C7"), stringsAsFactors = TRUE)

## Find DMRs ##
#tile into 100 bp windows with min coverage 10X and 5X
tiles.10X <- tileMethylCounts(myobj.3X, win.size = 100, step.size = 100, cov.bases = 10, save.db = TRUE, suffix = "tiles10X")
tiles.5X <- tileMethylCounts(myobj.3X, win.size = 100, step.size = 100, cov.bases = 5, save.db = TRUE, suffix = "tiles5X")

#check number of tiles 
tiles.10X
tiles.5X

#unite calls for 60% of samples
DMRmeth10X <- unite(tiles.10X, min.per.group=20L, save.db = TRUE, suffix = "DMRunite10X")
DMRmeth5X <- unite(tiles.5X, min.per.group=20L, save.db = TRUE, suffix = "DMRunite5X")

#check number of regions retained 
DMRmeth10X
DMRmeth5X

#calculate differential methylation 
DMRmyDiff10X <- calculateDiffMeth(DMRmeth10X, mc.cores=1, covariates=covariates, save.db = TRUE, suffix = "myDiff")
DMRmyDiff5X <- calculateDiffMeth(DMRmeth5X, mc.cores=1, covariates=covariates, save.db = TRUE, suffix = "myDiff")
 
#call significant methylation
DMRdiffMeth10X <- getMethylDiff(DMRmyDiff10X, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "diffMeth")
DMRdiffMeth5X <- getMethylDiff(DMRmyDiff5X, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "diffMeth")

#check number of DMRs
DMRdiffMeth10X
DMRdiffMeth5X

#get meth per chromosome
DMRdiffMethChr10X <- diffMethPerChr(DMRmyDiff10X, plot=FALSE,qvalue.cutoff=0.0125, meth.cutoff=15, save.db =  TRUE, suffix = "chr")
DMRdiffMethChr10X
DMRdiffMethChr5X <- diffMethPerChr(DMRmyDiff5X, plot=FALSE,qvalue.cutoff=0.0125, meth.cutoff=15, save.db =  TRUE, suffix = "chr")
DMRdiffMethChr5X

## Save R objects ##
saveRDS(DMRmeth10X, file = "./DMRmeth10X.RDs")
saveRDS(DMRmeth5X, file = "./DMRmeth5X.RDs")

saveRDS(DMRmyDiff10X, file = "./DMRmyDiff10X.RDs")
saveRDS(DMRmyDiff5X, file = "./DMRmyDiff5X.RDs")

saveRDS(DMRdiffMeth10X, file = "./DMRdiffMeth10X.RDS")
saveRDS(DMRdiffMeth5X, file = "./DMRdiffMeth5X.RDS")

saveRDS(DMRdiffMethChr10X, file = "./DMRdiffMethChr10X.RDS")
saveRDS(DMRdiffMethChr5X, file = "./DMRdiffMethChr5X.RDS")

## Get data ##
getData(DMRmeth10X) %>% saveRDS(file = "./DMRmeth10Xdata.RDS")
getData(DMRmeth5X) %>% saveRDS(file = "./DMRmeth5Xdata.RDS")

getData(DMRmyDiff10X) %>% saveRDS(file = "./DMRmyDiff10X.RDS")
getData(DMRmyDiff5X) %>% saveRDS(file = "./DMRmyDiff5X.RDS")

getData(DMRdiffMeth10X) %>% saveRDS(file = "./DMRdiffMeth10X.RDS")
getData(DMRdiffMeth5X) %>% saveRDS(file = "./DMRdiffMeth5X.RDS")






