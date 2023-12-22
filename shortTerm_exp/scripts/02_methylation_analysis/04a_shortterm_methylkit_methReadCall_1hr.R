#####################################################################################################################
### Goal: Read alignment files into methylKit and filter to create tabix files of filtered cytosine methylation
### Author: Janay Fox
### R script
#####################################################################################################################

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
library("genomation", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")

setwd("/scratch/janayfox/guppyWGBS/methylKit/st/perc20/1h")

# ## Prepare tabix files
# #create lists of file locations
# file.list.1h = list("../../../../mergedCov/st/ST2AC14F.CpG_merged.cov",
#                     "../../../../mergedCov/st/ST2AC14M.CpG_merged.cov",
#                     "../../../../mergedCov/st/ST2AC2F.CpG_merged.cov",
#                     "../../../../mergedCov/st/ST2AC2M.CpG_merged.cov",
#                     "../../../../mergedCov/st/ST2AC9F.CpG_merged.cov",
#                     "../../../../mergedCov/st/ST2AC9M.CpG_merged.cov",
#                     "../../../../mergedCov/st/ST2C14F.CpG_merged.cov",
#                     "../../../../mergedCov/st/ST2C14M.CpG_merged.cov",
#                     "../../../../mergedCov/st/ST2C3F.CpG_merged.cov",
#                     "../../../../mergedCov/st/ST2C3M.CpG_merged.cov",
#                     "../../../../mergedCov/st/ST2C9F.CpG_merged.cov",
#                     "../../../../mergedCov/st/ST2C9M.CpG_merged.cov")

# file.list.1h.fem = list("../../../../mergedCov/st/ST2AC14F.CpG_merged.cov",
#                         "../../../../mergedCov/st/ST2AC2F.CpG_merged.cov",
#                         "../../../../mergedCov/st/ST2AC9F.CpG_merged.cov",
#                         "../../../../mergedCov/st/ST2C14F.CpG_merged.cov",
#                         "../../../../mergedCov/st/ST2C3F.CpG_merged.cov",
#                         "../../../../mergedCov/st/ST2C9F.CpG_merged.cov")

# file.list.1h.mal = list("../../../../mergedCov/st/ST2AC14M.CpG_merged.cov",
#                         "../../../../mergedCov/st/ST2AC2M.CpG_merged.cov",
#                         "../../../../mergedCov/st/ST2AC9M.CpG_merged.cov",
#                         "../../../../mergedCov/st/ST2C14M.CpG_merged.cov",
#                         "../../../../mergedCov/st/ST2C3M.CpG_merged.cov",
#                         "../../../../mergedCov/st/ST2C9M.CpG_merged.cov")

# #create tabix file
# myobj.1h=methRead(file.list.1h,
#                   sample.id=list("ST2AC14F","ST2AC14M","ST2AC2F","ST2AC2M","ST2AC9F","ST2AC9M",
#                                  "ST2C14F","ST2C14M","ST2C3F","ST2C3M","ST2C9F", "ST2C9M"),
#                   assembly="guppyWGBS_shortterm",
#                   pipeline="bismarkCoverage",
#                   treatment=c(1,1,1,1,1,1,
#                               0,0,0,0,0,0),
#                   context="CpG",
#                   mincov = 1,
#                   dbtype = "tabix",
#                   dbdir = "shortterm_1h_DB_merged"
# )

# myobj.1h.fem=methRead(file.list.1h.fem,
#                       sample.id=list("ST2AC14F","ST2AC2F","ST2AC9F",
#                                      "ST2C14F","ST2C3F","ST2C9F"),
#                       assembly="guppyWGBS_shortterm",
#                       pipeline="bismarkCoverage",
#                       treatment=c(1,1,1,
#                                   0,0,0),
#                       context="CpG",
#                       mincov = 1,
#                       dbtype = "tabix",
#                       dbdir = "shortterm_1hF_DB_merged"
# )

# myobj.1h.mal=methRead(file.list.1h.mal,
#                       sample.id=list("ST2AC14M","ST2AC2M","ST2AC9M",
#                                      "ST2C14M","ST2C3M","ST2C9M"),
#                       assembly="guppyWGBS_shortterm",
#                       pipeline="bismarkCoverage",
#                       treatment=c(1,1,1,
#                                   0,0,0),
#                       context="CpG",
#                       mincov = 1,
#                       dbtype = "tabix",
#                       dbdir = "shortterm_1hM_DB_merged"
# )

# #get coverage stats 
# getCoverageStats(myobj.1h[[2]], both.strands = FALSE)
# getCoverageStats(myobj.1h.fem[[2]], both.strands = FALSE)
# getCoverageStats(myobj.1h.mal[[2]], both.strands = FALSE)

# #filter out sites in the 99.9th percentile of coverage (PCR bias) 
# myobj.1h.5X=filterByCoverage(myobj.1h,lo.count=5,lo.perc=NULL,
#                                hi.count=NULL, hi.perc=99.9, suffix = "1h_5X_merged")
# myobj.1h.fem.5X=filterByCoverage(myobj.1h.fem,lo.count=5,lo.perc=NULL,
#                                    hi.count=NULL, hi.perc=99.9, suffix = "1h_5X_fem_merged")
# myobj.1h.mal.5X=filterByCoverage(myobj.1h.mal,lo.count=5,lo.perc=NULL,
#                                    hi.count=NULL, hi.perc=99.9, suffix = "1h_5X_mal_merged")

# #normalize by median coverage
# norm.myobj.1h.5X=normalizeCoverage(myobj.1h.5X, method="median")
# norm.myobj.1h.fem.5X=normalizeCoverage(myobj.1h.fem.5X, method="median")
# norm.myobj.1h.mal.5X=normalizeCoverage(myobj.1h.mal.5X, method="median")

# #prepare GRanges object for chromosomes to keep 
# #to remove unplaced scaffolds and sex chromosomes
# keep.chr.noXY <- GRanges(seqnames = c("NC_024331.1", "NC_024332.1", "NC_024333.1", "NC_024334.1", "NC_024335.1",
#                                       "NC_024336.1", "NC_024337.1", "NC_024338.1", "NC_024339.1", "NC_024340.1",
#                                       "NC_024341.1", "NC_024343.1", "NC_024344.1", "NC_024345.1", "NC_024346.1",
#                                       "NC_024347.1", "NC_024348.1", "NC_024349.1", "NC_024350.1", "NC_024351.1",
#                                       "NC_024352.1", "NC_024353.1"),
#                          ranges=IRanges(start = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
#                                         end = c(34115677,46286544,35265442,31497199,33908744,
#                                                 31529174,31413364,27946405,34117797,32819797,
#                                                 28875558,33524197,28338960,30644713,33199742,
#                                                 30788009,22026651,28470737,26385442,25773841,
#                                                 25248790,18084596)),
#                          strand="*")

# seqlengths(keep.chr.noXY)=c(34115677,46286544,35265442,31497199,33908744,
#                             31529174,31413364,27946405,34117797,32819797,
#                             28875558,33524197,28338960,30644713,33199742,
#                             30788009,22026651,28470737,26385442,25773841,
#                             25248790,18084596)

# #to remove unplaced scaffolds
# keep.chr.allchr <- GRanges(seqnames = c("NC_024331.1", "NC_024332.1", "NC_024333.1", "NC_024334.1", "NC_024335.1",
#                                         "NC_024336.1", "NC_024337.1", "NC_024338.1", "NC_024339.1", "NC_024340.1",
#                                         "NC_024341.1", "NC_024342.1", "NC_024343.1", "NC_024344.1", "NC_024345.1", "NC_024346.1",
#                                         "NC_024347.1", "NC_024348.1", "NC_024349.1", "NC_024350.1", "NC_024351.1",
#                                         "NC_024352.1", "NC_024353.1"),
#                            ranges=IRanges(start = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
#                                           end = c(34115677,46286544,35265442,31497199,33908744,
#                                                   31529174,31413364,27946405,34117797,32819797,
#                                                   28875558,26439574,33524197,28338960,30644713,33199742,
#                                                   30788009,22026651,28470737,26385442,25773841,
#                                                   25248790,18084596)),
#                            strand="*")

# seqlengths(keep.chr.allchr)=c(34115677,46286544,35265442,31497199,33908744,
#                               31529174,31413364,27946405,34117797,32819797,
#                               28875558,26439574,33524197,28338960,30644713,33199742,
#                               30788009,22026651,28470737,26385442,25773841,
#                               25248790,18084596)

# #remove sex chr (LG12) and unplaced scaffolds
# myobj.1h.subset <- selectByOverlap(norm.myobj.1h.5X, keep.chr.noXY)
# myobj.1h.fem.subset <- selectByOverlap(norm.myobj.1h.fem.5X, keep.chr.allchr)
# myobj.1h.mal.subset <- selectByOverlap(norm.myobj.1h.mal.5X, keep.chr.allchr)

# ##Find DMS##
# #unite sites 
# DMS.meth.1h.5X=unite(myobj.1h.subset, min.per.group = 6L, destrand=FALSE, save.db = TRUE, suffix = "DMS_unite_1h")
# DMS.meth.1h.fem.5X=unite(myobj.1h.fem.subset, min.per.group = 3L, destrand=FALSE, save.db = TRUE, suffix = "DMS_unite_1h_fem")
# DMS.meth.1h.mal.5X=unite(myobj.1h.mal.subset, min.per.group = 3L, destrand=FALSE, save.db = TRUE, suffix = "DMS_unite_1h_mal")

# #convert to non DB object 
# DMS.meth.1h.5X <- as(DMS.meth.1h.5X, "methylBase")
# DMS.meth.1h.fem.5X <- as(DMS.meth.1h.fem.5X, "methylBase")
# DMS.meth.1h.mal.5X <- as(DMS.meth.1h.mal.5X, "methylBase")

# #filter out low variation sites 
# pm.1h.5x <- percMethylation(DMS.meth.1h.5X) #get percent methylation matrix
# sds.1h.5x <- matrixStats::rowSds(pm.1h.5x) #calculate standard deviation of CpGs 
# DMS.meth.1h.5X <- DMS.meth.1h.5X[sds.1h.5x > 2]

# pm.1h.fem.5x <- percMethylation(DMS.meth.1h.fem.5X) #get percent methylation matrix
# sds.1h.fem.5x <- matrixStats::rowSds(pm.1h.fem.5x) #calculate standard deviation of CpGs 
# DMS.meth.1h.fem.5X <- DMS.meth.1h.fem.5X[sds.1h.fem.5x > 2]

# pm.1h.mal.5x <- percMethylation(DMS.meth.1h.mal.5X) #get percent methylation matrix
# sds.1h.mal.5x <- matrixStats::rowSds(pm.1h.mal.5x) #calculate standard deviation of CpGs 
# DMS.meth.1h.mal.5X <- DMS.meth.1h.mal.5X[sds.1h.mal.5x > 2]

# #filter out SNPs
# snp <- read.csv("../../../../BS-SNPer/shortterm_CT_SNP_edit.csv") #read in snps
# snp.granges <- makeGRangesFromDataFrame(snp, ignore.strand = TRUE) #convert to granges 

# DMS.meth.1h.5X <- DMS.meth.1h.5X[!as(DMS.meth.1h.5X, "GRanges") %over% snp.granges, ] #select CpGs that do not overlap
# DMS.meth.1h.fem.5X <- DMS.meth.1h.fem.5X[!as(DMS.meth.1h.fem.5X, "GRanges") %over% snp.granges, ] #select CpGs that do not overlap
# DMS.meth.1h.mal.5X <- DMS.meth.1h.mal.5X[!as(DMS.meth.1h.mal.5X, "GRanges") %over% snp.granges, ] #select CpGs that do not overlap

# # Check number of CpGs 
# DMS.meth.1h.5X
# DMS.meth.1h.fem.5X
# DMS.meth.1h.mal.5X

# #enter covariates 
covariates.1h <- data.frame(tank=c("AC14","AC14","AC2","AC2","AC9","AC9",
                                   "C14","C14","C3","C3","C9","C9"), 
                            sex=c("F","M","F","M","F","M",
                                  "F","M","F","M","F","M"), 
                            stringsAsFactors = TRUE)

# #calculate differential methylation
# DMS.myDiff.1h.5X <- calculateDiffMeth(DMS.meth.1h.5X, covariates=covariates.1h, mc.cores=2, test="Chisq", save.db = TRUE, suffix = "DMS_myDiff_1h")
# DMS.myDiff.1h.fem.5X <- calculateDiffMeth(DMS.meth.1h.fem.5X, mc.cores=2, test="Chisq", save.db = TRUE, suffix = "DMS_myDiff_1h_fem")
# DMS.myDiff.1h.mal.5X <- calculateDiffMeth(DMS.meth.1h.mal.5X, mc.cores=2, test="Chisq", save.db = TRUE, suffix = "DMS_myDiff_1h_mal")

# #call significant methylation
# DMS.diffMeth.1h.5X <- getMethylDiff(DMS.myDiff.1h.5X, difference = 20, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_1h")
# DMS.diffMeth.1h.fem.5X <- getMethylDiff(DMS.myDiff.1h.fem.5X, difference = 20, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_1h_fem")
# DMS.diffMeth.1h.mal.5X <- getMethylDiff(DMS.myDiff.1h.mal.5X, difference = 20, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_1h_mal")

# #check number of significant DMS
# DMS.diffMeth.1h.5X
# DMS.diffMeth.1h.fem.5X
# DMS.diffMeth.1h.mal.5X

# # Get meth per chromosome
# DMS.diffMethChr.1h.5X <- diffMethPerChr(DMS.myDiff.1h.5X, plot=FALSE, qvalue.cutoff=0.0125, meth.cutoff=20, save.db = TRUE, suffix = "DMS_chr_1h")
# DMS.diffMethChr.1h.5X
# DMS.diffMethChr.1h.fem.5X <- diffMethPerChr(DMS.myDiff.1h.fem.5X, plot=FALSE, qvalue.cutoff=0.0125, meth.cutoff=20, save.db = TRUE, suffix = "DMS_chr_1h_fem")
# DMS.diffMethChr.1h.fem.5X
# DMS.diffMethChr.1h.mal.5X <- diffMethPerChr(DMS.myDiff.1h.mal.5X, plot=FALSE, qvalue.cutoff=0.0125, meth.cutoff=20, save.db = TRUE, suffix = "DMS_chr_1h_mal")
# DMS.diffMethChr.1h.mal.5X

DMS.meth.1h.5X <- readRDS("./DMS_res_1h/DMSmeth_1h_5X.RDS")
DMS.meth.1h.fem.5X <- readRDS("./DMS_res_1h/DMSmeth_1h_fem_5X.RDS")
DMS.meth.1h.mal.5X <- readRDS("./DMS_res_1h/DMSmeth_1h_mal_5X.RDS")

##Find DMRs ##
#tile into 100 bp windows with min coverage 5X
DMR.meth.1h.5X <- tileMethylCounts(DMS.meth.1h.5X, win.size = 100, step.size = 100, cov.bases = 5)
DMR.meth.1h.fem.5X <- tileMethylCounts(DMS.meth.1h.fem.5X, win.size = 100, step.size = 100, cov.bases = 5)
DMR.meth.1h.mal.5X <- tileMethylCounts(DMS.meth.1h.mal.5X, win.size = 100, step.size = 100, cov.bases = 5)

#check number of regions retained 
DMR.meth.1h.5X
DMR.meth.1h.fem.5X
DMR.meth.1h.mal.5X

#calculate differential methylation 
DMR.myDiff.1h.5X <- calculateDiffMeth(DMR.meth.1h.5X, mc.cores=2, test="Chisq", covariates=covariates.1h, save.db = TRUE, suffix = "DMR_myDiff_1h_5X")
DMR.myDiff.1h.fem.5X <- calculateDiffMeth(DMR.meth.1h.fem.5X, mc.cores=2, test="Chisq", save.db = TRUE, suffix = "DMR_myDiff_1h_fem_5X")
DMR.myDiff.1h.mal.5X <- calculateDiffMeth(DMR.meth.1h.mal.5X, mc.cores=2, test="Chisq", save.db = TRUE, suffix = "DMR_myDiff_1h_mal_5X")

#call significant methylation
DMR.diffMeth.1h.5X <- getMethylDiff(DMR.myDiff.1h.5X, difference = 20, qvalue = 0.0125, save.db = TRUE, suffix = "DMR_diffMeth_1h_5X")
DMR.diffMeth.1h.fem.5X <- getMethylDiff(DMR.myDiff.1h.fem.5X, difference = 20, qvalue = 0.0125, save.db = TRUE, suffix = "DMR_diffMeth_1h_fem_5X")
DMR.diffMeth.1h.mal.5X <- getMethylDiff(DMR.myDiff.1h.mal.5X, difference = 20, qvalue = 0.0125, save.db = TRUE, suffix = "DMR_diffMeth_1h_mal_5X")

#check number of DMRs 
DMR.diffMeth.1h.5X
DMR.diffMeth.1h.fem.5X
DMR.diffMeth.1h.mal.5X

#get meth per chromosome 
DMR.diffMeth.1h.5X.chr <- diffMethPerChr(DMR.diffMeth.1h.5X, plot = FALSE,qvalue.cutoff=0.0125, meth.cutoff=20, save.db =  TRUE, suffix = "chrDMR_1h_5X")
DMR.diffMeth.1h.5X.chr
DMR.diffMeth.1h.fem.5X.chr <- diffMethPerChr(DMR.diffMeth.1h.fem.5X, plot = FALSE,qvalue.cutoff=0.0125, meth.cutoff=20, save.db =  TRUE, suffix = "chrDMR_1h_fem_5X")
DMR.diffMeth.1h.fem.5X.chr
DMR.diffMeth.1h.mal.5X.chr <- diffMethPerChr(DMR.diffMeth.1h.mal.5X, plot = FALSE,qvalue.cutoff=0.0125, meth.cutoff=20, save.db =  TRUE, suffix = "chrDMR_1h_mal_5X")
DMR.diffMeth.1h.mal.5X.chr

## Save R objects ##
# saveRDS(myobj.1h.subset, file = "./DMS_res/myobj_1h_5X.RDS")
# saveRDS(myobj.1h.fem.subset, file = "./DMS_res/myobj_1h_fem_5X.RDS")
# saveRDS(myobj.1h.mal.subset, file = "./DMS_res/myobj_1h_mal_5X.RDS")

# saveRDS(DMS.meth.1h.5X, file = "./DMS_res/DMSmeth_1h_5X.RDS")
# saveRDS(DMS.meth.1h.fem.5X, file = "./DMS_res/DMSmeth_1h_fem_5X.RDS")
# saveRDS(DMS.meth.1h.mal.5X, file = "./DMS_res/DMSmeth_1h_mal_5X.RDS")

# saveRDS(DMS.myDiff.1h.5X, file = "./DMS_res/DMSmydiff_1h_5X.RDS")
# saveRDS(DMS.myDiff.1h.fem.5X, file = "./DMS_res/DMSmydiff_1h_fem_5X.RDS")
# saveRDS(DMS.myDiff.1h.mal.5X, file = "./DMS_res/DMSmydiff_1h_mal_5X.RDS")

# saveRDS(DMS.diffMeth.1h.5X, file = "./DMS_res/DMSdiffmeth_1h_5X.RDS")
# saveRDS(DMS.diffMeth.1h.fem.5X, file = "./DMS_res/DMSdiffmeth_1h_fem_5X.RDS")
# saveRDS(DMS.diffMeth.1h.mal.5X, file = "./DMS_res/DMSdiffmeth_1h_mal_5X.RDS")

# saveRDS(getData(DMS.meth.1h.5X), file = "./DMS_res/meth_1h_5X_data.RDS")
# saveRDS(getData(DMS.meth.1h.fem.5X), file = "./DMS_res/meth_1h_fem_5X_data.RDS")
# saveRDS(getData(DMS.meth.1h.mal.5X), file = "./DMS_res/meth_1h_mal_5X_data.RDS")

# saveRDS(getData(DMS.diffMeth.1h.5X), file = "./DMS_res/DMSdiffmeth_1h_5X_data.RDS")
# saveRDS(getData(DMS.diffMeth.1h.fem.5X), file = "./DMS_res/DMSdiffmeth_1h_fem_5X_data.RDS")
# saveRDS(getData(DMS.diffMeth.1h.mal.5X), file = "./DMS_res/DMSdiffmeth_1h_mal_5X_data.RDS")

# saveRDS(getData(DMS.myDiff.1h.5X), file = "./DMS_res/DMSmydiff_1h_5X_data.RDS")
# saveRDS(getData(DMS.myDiff.1h.fem.5X), file = "./DMS_res/DMSmydiff_1h_fem_5X_data.RDS")
# saveRDS(getData(DMS.myDiff.1h.mal.5X), file = "./DMS_res/DMSmydiff_1h_mal_5X_data.RDS")

saveRDS(DMR.meth.1h.5X, file = "./DMR_res_1h_updated/DMRmeth_1h_5X.RDS")
saveRDS(DMR.meth.1h.fem.5X, file = "./DMR_res_1h_updated/DMRmeth_1h_fem_5X.RDS")
saveRDS(DMR.meth.1h.mal.5X, file = "./DMR_res_1h_updated/DMRmeth_1h_mal_5X.RDS")

saveRDS(getData(DMR.meth.1h.5X), file = "./DMR_res_1h_updated/DMRmeth_1h_5X_data.RDS")
saveRDS(getData(DMR.meth.1h.fem.5X), file = "./DMR_res_1h_updated/DMRmeth_1h_fem_5X_data.RDS")
saveRDS(getData(DMR.meth.1h.mal.5X), file = "./DMR_res_1h_updated/DMRmeth_1h_mal_5X_data.RDS")

saveRDS(DMR.myDiff.1h.5X, file = "./DMR_res_1h_updated/DMRmydiff_1h_5X.RDS")
saveRDS(DMR.myDiff.1h.fem.5X, file = "./DMR_res_1h_updated/DMRmydiff_1h_fem_5X.RDS")
saveRDS(DMR.myDiff.1h.mal.5X, file = "./DMR_res_1h_updated/DMRmydiff_1h_mal_5X.RDS")

saveRDS(getData(DMR.myDiff.1h.5X), file = "./DMR_res_1h_updated/DMRmydiff_1h_5X_data.RDS")
saveRDS(getData(DMR.myDiff.1h.fem.5X), file = "./DMR_res_1h_updated/DMRmydiff_1h_fem_5X_data.RDS")
saveRDS(getData(DMR.myDiff.1h.mal.5X), file = "./DMR_res_1h_updated/DMRmydiff_1h_mal_5X_data.RDS")

saveRDS(DMR.diffMeth.1h.5X, file = "./DMR_res_1h_updated/DMRdiffmeth_1h_5X.RDS")
saveRDS(DMR.diffMeth.1h.fem.5X, file = "./DMR_res_1h_updated/DMRdiffmeth_1h_fem_5X.RDS")
saveRDS(DMR.diffMeth.1h.mal.5X, file = "./DMR_res_1h_updated/DMRdiffmeth_1h_mal_5X.RDS")

saveRDS(getData(DMR.diffMeth.1h.5X), file = "./DMR_res_1h_updated/DMRdiffmeth_1h_5X_data.RDS")
saveRDS(getData(DMR.diffMeth.1h.fem.5X), file = "./DMR_res_1h_updated/DMRdiffmeth_1h_fem_5X_data.RDS")
saveRDS(getData(DMR.diffMeth.1h.mal.5X), file = "./DMR_res_1h_updated/DMRdiffmeth_1h_mal_5X_data.RDS")

saveRDS(DMR.diffMeth.1h.5X.chr, file = "./DMR_res_1h_updated/DMRdiffmethchr_1h_5X.RDS")
saveRDS(DMR.diffMeth.1h.fem.5X.chr, file = "./DMR_res_1h_updated/DMRdiffmethchr_1h_fem_5X.RDS")
saveRDS(DMR.diffMeth.1h.mal.5X.chr, file = "./DMR_res_1h_updated/DMRdiffmethchr_1h_mal_5X.RDS")
