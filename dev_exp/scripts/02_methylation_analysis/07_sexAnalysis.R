#####################################################################################################################
### Goal: Run methylation analysis to look for sex differences 
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
library("data.table", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")

setwd("/scratch/janayfox/guppyWGBS/methylKit/sexDiff")

file.list.dev = list("../../mergedCov/dev/DC2F1.CpG_merged.cov",
                     "../../mergedCov/dev/DC2F2.CpG_merged.cov",
                     "../../mergedCov/dev/DC3F1.CpG_merged.cov",
                     "../../mergedCov/dev/DC3F2.CpG_merged.cov",
                     "../../mergedCov/dev/DC3F3.CpG_merged.cov",
                     "../../mergedCov/dev/DC4F1.CpG_merged.cov",
                     "../../mergedCov/dev/DC4F2.CpG_merged.cov",
                     "../../mergedCov/dev/DC4F3.CpG_merged.cov",
                     "../../mergedCov/dev/DC4F4.CpG_merged.cov",
                     "../../mergedCov/dev/DC4F5.CpG_merged.cov",
                     "../../mergedCov/dev/DC5F1.CpG_merged.cov",
                     "../../mergedCov/dev/DC5F2.CpG_merged.cov",
                     "../../mergedCov/dev/DC5F3.CpG_merged.cov",
                     "../../mergedCov/dev/DC5F4.CpG_merged.cov",
                     "../../mergedCov/dev/DC5F5.CpG_merged.cov",
                     "../../mergedCov/dev/DC6F1.CpG_merged.cov",
                     "../../mergedCov/dev/DC6F2.CpG_merged.cov",
                     "../../mergedCov/dev/DC6F3.CpG_merged.cov",
                     "../../mergedCov/dev/DC6F4.CpG_merged.cov",
                     "../../mergedCov/dev/DC7F1.CpG_merged.cov",
                     "../../mergedCov/dev/DC7F2.CpG_merged.cov",
                     "../../mergedCov/dev/DC7F3.CpG_merged.cov",
                     "../../mergedCov/dev/DC2M1.CpG_merged.cov",
                     "../../mergedCov/dev/DC2M2.CpG_merged.cov",
                     "../../mergedCov/dev/DC3M1.CpG_merged.cov",
                     "../../mergedCov/dev/DC3M2.CpG_merged.cov",
                     "../../mergedCov/dev/DC3M4.CpG_merged.cov",
                     "../../mergedCov/dev/DC4M1.CpG_merged.cov",
                     "../../mergedCov/dev/DC4M2.CpG_merged.cov",
                     "../../mergedCov/dev/DC5M1.CpG_merged.cov",
                     "../../mergedCov/dev/DC5M2.CpG_merged.cov",
                     "../../mergedCov/dev/DC5M3.CpG_merged.cov",
                     "../../mergedCov/dev/DC6M1.CpG_merged.cov",
                     "../../mergedCov/dev/DC6M2.CpG_merged.cov",
                     "../../mergedCov/dev/DC6M4.CpG_merged.cov",
                     "../../mergedCov/dev/DC7M1.CpG_merged.cov",
                     "../../mergedCov/dev/DC7M2.CpG_merged.cov",
                     "../../mergedCov/dev/DC7M3.CpG_merged.cov",
                     "../../mergedCov/dev/DC7M4.CpG_merged.cov",
                     "../../mergedCov/dev/DC7M5.CpG_merged.cov")


file.list.st = list("../../mergedCov/st/ST2C10F.CpG_merged.cov",
                         "../../mergedCov/st/ST2C11F.CpG_merged.cov",
                         "../../mergedCov/st/ST2C12F.CpG_merged.cov",
                         "../../mergedCov/st/ST2C13F.CpG_merged.cov",
                         "../../mergedCov/st/ST2C14F.CpG_merged.cov",
                         "../../mergedCov/st/ST2C15F.CpG_merged.cov",
                         "../../mergedCov/st/ST2C1F.CpG_merged.cov",
                         "../../mergedCov/st/ST2C2F.CpG_merged.cov",
                         "../../mergedCov/st/ST2C3F.CpG_merged.cov",
                         "../../mergedCov/st/ST2C4F.CpG_merged.cov",
                         "../../mergedCov/st/ST2C5F.CpG_merged.cov",
                         "../../mergedCov/st/ST2C6F.CpG_merged.cov",
                         "../../mergedCov/st/ST2C7F.CpG_merged.cov",
                         "../../mergedCov/st/ST2C8F.CpG_merged.cov",
                         "../../mergedCov/st/ST2C9F.CpG_merged.cov",
                         "../../mergedCov/st/ST2C10M.CpG_merged.cov",
                         "../../mergedCov/st/ST2C11M.CpG_merged.cov",
                         "../../mergedCov/st/ST2C12M.CpG_merged.cov",
                         "../../mergedCov/st/ST2C13M.CpG_merged.cov",
                         "../../mergedCov/st/ST2C14M.CpG_merged.cov",
                         "../../mergedCov/st/ST2C15M.CpG_merged.cov",
                         "../../mergedCov/st/ST2C1M.CpG_merged.cov",
                         "../../mergedCov/st/ST2C2M.CpG_merged.cov",
                         "../../mergedCov/st/ST2C3M.CpG_merged.cov",
                         "../../mergedCov/st/ST2C4M.CpG_merged.cov",
                         "../../mergedCov/st/ST2C5M.CpG_merged.cov",
                         "../../mergedCov/st/ST2C6M.CpG_merged.cov",
                         "../../mergedCov/st/ST2C7M.CpG_merged.cov",
                         "../../mergedCov/st/ST2C8M.CpG_merged.cov",
                         "../../mergedCov/st/ST2C9M.CpG_merged.cov")

file.list.all = list("../../mergedCov/dev/DC2F1.CpG_merged.cov",
                     "../../mergedCov/dev/DC2F2.CpG_merged.cov",
                     "../../mergedCov/dev/DC3F1.CpG_merged.cov",
                     "../../mergedCov/dev/DC3F2.CpG_merged.cov",
                     "../../mergedCov/dev/DC3F3.CpG_merged.cov",
                     "../../mergedCov/dev/DC4F1.CpG_merged.cov",
                     "../../mergedCov/dev/DC4F2.CpG_merged.cov",
                     "../../mergedCov/dev/DC4F3.CpG_merged.cov",
                     "../../mergedCov/dev/DC4F4.CpG_merged.cov",
                     "../../mergedCov/dev/DC4F5.CpG_merged.cov",
                     "../../mergedCov/dev/DC5F1.CpG_merged.cov",
                     "../../mergedCov/dev/DC5F2.CpG_merged.cov",
                     "../../mergedCov/dev/DC5F3.CpG_merged.cov",
                     "../../mergedCov/dev/DC5F4.CpG_merged.cov",
                     "../../mergedCov/dev/DC5F5.CpG_merged.cov",
                     "../../mergedCov/dev/DC6F1.CpG_merged.cov",
                     "../../mergedCov/dev/DC6F2.CpG_merged.cov",
                     "../../mergedCov/dev/DC6F3.CpG_merged.cov",
                     "../../mergedCov/dev/DC6F4.CpG_merged.cov",
                     "../../mergedCov/dev/DC7F1.CpG_merged.cov",
                     "../../mergedCov/dev/DC7F2.CpG_merged.cov",
                     "../../mergedCov/dev/DC7F3.CpG_merged.cov",
                     "../../mergedCov/dev/DC2M1.CpG_merged.cov",
                     "../../mergedCov/dev/DC2M2.CpG_merged.cov",
                     "../../mergedCov/dev/DC3M1.CpG_merged.cov",
                     "../../mergedCov/dev/DC3M2.CpG_merged.cov",
                     "../../mergedCov/dev/DC3M4.CpG_merged.cov",
                     "../../mergedCov/dev/DC4M1.CpG_merged.cov",
                     "../../mergedCov/dev/DC4M2.CpG_merged.cov",
                     "../../mergedCov/dev/DC5M1.CpG_merged.cov",
                     "../../mergedCov/dev/DC5M2.CpG_merged.cov",
                     "../../mergedCov/dev/DC5M3.CpG_merged.cov",
                     "../../mergedCov/dev/DC6M1.CpG_merged.cov",
                     "../../mergedCov/dev/DC6M2.CpG_merged.cov",
                     "../../mergedCov/dev/DC6M4.CpG_merged.cov",
                     "../../mergedCov/dev/DC7M1.CpG_merged.cov",
                     "../../mergedCov/dev/DC7M2.CpG_merged.cov",
                     "../../mergedCov/dev/DC7M3.CpG_merged.cov",
                     "../../mergedCov/dev/DC7M4.CpG_merged.cov",
                     "../../mergedCov/dev/DC7M5.CpG_merged.cov",
                     "../../mergedCov/st/ST2C10F.CpG_merged.cov",
                    "../../mergedCov/st/ST2C11F.CpG_merged.cov",
                    "../../mergedCov/st/ST2C12F.CpG_merged.cov",
                    "../../mergedCov/st/ST2C13F.CpG_merged.cov",
                    "../../mergedCov/st/ST2C14F.CpG_merged.cov",
                    "../../mergedCov/st/ST2C15F.CpG_merged.cov",
                    "../../mergedCov/st/ST2C1F.CpG_merged.cov",
                    "../../mergedCov/st/ST2C2F.CpG_merged.cov",
                    "../../mergedCov/st/ST2C3F.CpG_merged.cov",
                    "../../mergedCov/st/ST2C4F.CpG_merged.cov",
                    "../../mergedCov/st/ST2C5F.CpG_merged.cov",
                    "../../mergedCov/st/ST2C6F.CpG_merged.cov",
                    "../../mergedCov/st/ST2C7F.CpG_merged.cov",
                    "../../mergedCov/st/ST2C8F.CpG_merged.cov",
                    "../../mergedCov/st/ST2C9F.CpG_merged.cov",
                    "../../mergedCov/st/ST2C10M.CpG_merged.cov",
                    "../../mergedCov/st/ST2C11M.CpG_merged.cov",
                    "../../mergedCov/st/ST2C12M.CpG_merged.cov",
                    "../../mergedCov/st/ST2C13M.CpG_merged.cov",
                    "../../mergedCov/st/ST2C14M.CpG_merged.cov",
                    "../../mergedCov/st/ST2C15M.CpG_merged.cov",
                    "../../mergedCov/st/ST2C1M.CpG_merged.cov",
                    "../../mergedCov/st/ST2C2M.CpG_merged.cov",
                    "../../mergedCov/st/ST2C3M.CpG_merged.cov",
                    "../../mergedCov/st/ST2C4M.CpG_merged.cov",
                    "../../mergedCov/st/ST2C5M.CpG_merged.cov",
                    "../../mergedCov/st/ST2C6M.CpG_merged.cov",
                    "../../mergedCov/st/ST2C7M.CpG_merged.cov",
                    "../../mergedCov/st/ST2C8M.CpG_merged.cov",
                    "../../mergedCov/st/ST2C9M.CpG_merged.cov")

myobj.dev=methRead(file.list.dev.fem,
               sample.id=list("DC2F1","DC2F2","DC3F1","DC3F2","DC3F3",
                              "DC4F1","DC4F2","DC4F3","DC4F4","DC4F5",
                              "DC5F1","DC5F2","DC5F3","DC5F4","DC5F5",
                              "DC6F1","DC6F2","DC6F3", "DC6F4","DC7F1",
                              "DC7F2","DC7F3",
                              "DC2M1","DC2M2","DC3M1","DC3M2","DC3M4",
                              "DC4M1","DC4M2","DC5M1","DC5M2","DC5M3",
                              "DC6M1","DC6M2","DC6M4","DC7M1","DC7M2",
                              "DC7M3","DC7M4","DC7M5"),
               assembly="guppyWGBS_sex_dev",
               pipeline="bismarkCoverage",
               treatment=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
               context="CpG",
              # dbtype = "tabix",
               dbdir = "guppy_sex_dev",
               mincov = 3
)

myobj.st=methRead(file.list.all.fem,
                       sample.id=list("ST2C10F","ST2C11F","ST2C12F","ST2C13F","ST2C14F",
                                      "ST2C15F","ST2C1F","ST2C2F","ST2C3F","ST2C4F",
                                      "ST2C5F","ST2C6F","ST2C7F","ST2C8F","ST2C9F",
                                      "ST2C10M","ST2C11M","ST2C12M","ST2C13M","ST2C14M",
                                      "ST2C15M","ST2C1M","ST2C2M","ST2C3M","ST2C4M",
                                      "ST2C5M","ST2C6M","ST2C7M","ST2C8M", "ST2C9M"),
                       assembly="guppyWGBS_sex_st",
                       pipeline="bismarkCoverage",
                       treatment=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                       context="CpG",
                       mincov = 3,
                      # dbtype = "tabix",
                       dbdir = "shortterm_allfem_DB"
)

myobj.all=methRead(file.list.dev.fem,
               sample.id=list("DC2F1","DC2F2","DC3F1","DC3F2","DC3F3",
                              "DC4F1","DC4F2","DC4F3","DC4F4","DC4F5",
                              "DC5F1","DC5F2","DC5F3","DC5F4","DC5F5",
                              "DC6F1","DC6F2","DC6F3", "DC6F4","DC7F1",
                              "DC7F2","DC7F3","ST2C10F","ST2C11F","ST2C12F",
                              "ST2C13F","ST2C14F","ST2C15F","ST2C1F","ST2C2F",
                              "ST2C3F","ST2C4F","ST2C5F","ST2C6F","ST2C7F",
                              "ST2C8F","ST2C9F",
                              "DC2M1","DC2M2","DC3M1","DC3M2","DC3M4",
                              "DC4M1","DC4M2","DC5M1","DC5M2","DC5M3",
                              "DC6M1","DC6M2","DC6M4","DC7M1","DC7M2",
                              "DC7M3","DC7M4","DC7M5","ST2C10M","ST2C11M",
                              "ST2C12M","ST2C13M","ST2C14M","ST2C15M","ST2C1M",
                              "ST2C2M","ST2C3M","ST2C4M","ST2C5M","ST2C6M",
                              "ST2C7M","ST2C8M", "ST2C9M"),
               assembly="guppyWGBS_sex_all",
               pipeline="bismarkCoverage",
               treatment=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
               context="CpG",
              # dbtype = "tabix",
               dbdir = "guppy_sex_all",
               mincov = 3
)

#filter out sites in the 99.9th percentile of coverage (PCR bias) and 5x coverage
myobj.dev.3X=filterByCoverage(myobj.dev,lo.count=3,lo.perc=NULL,
                          hi.count=NULL, hi.perc=99.9)

myobj.st.3X=filterByCoverage(myobj.st,lo.count=3,lo.perc=NULL,
                          hi.count=NULL, hi.perc=99.9)

myobj.all.3X=filterByCoverage(myobj.all,lo.count=3,lo.perc=NULL,
                          hi.count=NULL, hi.perc=99.9)

myobj.dev.5X=filterByCoverage(myobj.dev,lo.count=5,lo.perc=NULL,
                          hi.count=NULL, hi.perc=99.9)

myobj.st.5X=filterByCoverage(myobj.st,lo.count=5,lo.perc=NULL,
                          hi.count=NULL, hi.perc=99.9)

myobj.all.5X=filterByCoverage(myobj.all,lo.count=5,lo.perc=NULL,
                          hi.count=NULL, hi.perc=99.9)

#normalize by median coverage
norm.myobj.dev.3X=normalizeCoverage(myobj.dev.3X, method="median")
norm.myobj.st.3X=normalizeCoverage(myobj.st.3X, method="median")
norm.myobj.all.3X=normalizeCoverage(myobj.all.3X, method="median")

norm.myobj.dev.5X=normalizeCoverage(myobj.dev.5X, method="median")
norm.myobj.st.5X=normalizeCoverage(myobj.st.5X, method="median")
norm.myobj.all.5X=normalizeCoverage(myobj.all.5X, method="median")

#to remove unplaced scaffolds
keep.chr.allchr <- GRanges(seqnames = c("NC_024331.1", "NC_024332.1", "NC_024333.1", "NC_024334.1", "NC_024335.1",
                                        "NC_024336.1", "NC_024337.1", "NC_024338.1", "NC_024339.1", "NC_024340.1",
                                        "NC_024341.1", "NC_024342.1", "NC_024343.1", "NC_024344.1", "NC_024345.1", "NC_024346.1",
                                        "NC_024347.1", "NC_024348.1", "NC_024349.1", "NC_024350.1", "NC_024351.1",
                                        "NC_024352.1", "NC_024353.1"),
                           ranges=IRanges(start = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
                                          end = c(34115677,46286544,35265442,31497199,33908744,
                                                  31529174,31413364,27946405,34117797,32819797,
                                                  28875558,26439574,33524197,28338960,30644713,33199742,
                                                  30788009,22026651,28470737,26385442,25773841,
                                                  25248790,18084596)),
                           strand="*")

seqlengths(keep.chr.allchr)=c(34115677,46286544,35265442,31497199,33908744,
                              31529174,31413364,27946405,34117797,32819797,
                              28875558,26439574,33524197,28338960,30644713,33199742,
                              30788009,22026651,28470737,26385442,25773841,
                              25248790,18084596)

myobj3X.dev.subset <- selectByOverlap(norm.myobj.dev.3X, keep.chr.noXY)
myobj3X.st.subset <- selectByOverlap(norm.myobj.st.3X, keep.chr.noXY)
myobj3X.all.subset <- selectByOverlap(norm.myobj.all.3X, keep.chr.noXY)

myobj5X.dev.subset <- selectByOverlap(norm.myobj.dev.5X, keep.chr.noXY)
myobj5X.st.subset <- selectByOverlap(norm.myobj.st.5X, keep.chr.noXY)
myobj5X.all.subset <- selectByOverlap(norm.myobj.all.5X, keep.chr.noXY)

## Find DMSs ## 
# Unite methylation calls
DMSmeth5X.dev <- unite(myobj5X.dev.subset, destrand=FALSE, save.db = FALSE)
DMSmeth5X.st <- unite(myobj5X.st.subset, destrand=FALSE, save.db = FALSE)
DMSmeth5X.all <- unite(myobj5X.all.subset, destrand=FALSE, save.db = FALSE)

#filter out low variation sites 
pm <- percMethylation(DMSmeth5X.dev) #get percent methylation matrix
sds <- matrixStats::rowSds(pm, na.rm = TRUE) #calculate standard deviation of CpGs 
DMSmeth5X.dev <- DMSmeth5X.dev[sds > 2]

pm <- percMethylation(DMSmeth5X.st) #get percent methylation matrix
sds <- matrixStats::rowSds(pm, na.rm = TRUE) #calculate standard deviation of CpGs 
DMSmeth5X.st <- DMSmeth5X.st[sds > 2]

pm <- percMethylation(DMSmeth5X.all) #get percent methylation matrix
sds <- matrixStats::rowSds(pm, na.rm = TRUE) #calculate standard deviation of CpGs 
DMSmeth5X.all <- DMSmeth5X.all[sds > 2]

#filter out SNPs
snp <- read.csv("../../BS-SNPer/dev_CT_SNP_edit.csv") #read in snps
snp.granges <- makeGRangesFromDataFrame(snp, ignore.strand = TRUE) #convert to granges 
DMSmeth5X.dev <- DMSmeth5X.dev[!as(DMSmeth5X.dev, "GRanges") %over% snp.granges, ] #select CpGs that do not overlap
DMSmeth5X.st <- DMSmeth5X.st[!as(DMSmeth5X.st, "GRanges") %over% snp.granges, ] #select CpGs that do not overlap
DMSmeth5X.all <- DMSmeth5X.all[!as(DMSmeth5X.all, "GRanges") %over% snp.granges, ] #select CpGs that do not overlap

#check number of cpgs
DMSmeth5X.dev
DMSmeth5X.st
DMSmeth5X.all

#read in covariates
covariates.pop <- data.frame(pop = c("paria","paria","paria","paria","paria",
                                    "paria","paria","paria","paria","paria",
                                    "paria","paria","paria","paria","paria",
                                    "paria","paria","paria","paria","paria",
                                    "paria","paria",
                                    "aripo","aripo","aripo","aripo","aripo",
                                    "aripo","aripo","aripo","aripo","aripo",
                                    "aripo","aripo","aripo","aripo","aripo",
                                    "paria","paria","paria","paria","paria",
                                    "paria","paria","paria","paria","paria",
                                    "paria","paria","paria","paria","paria",
                                    "paria","paria","paria",
                                    "aripo","aripo","aripo","aripo","aripo",
                                    "aripo","aripo","aripo","aripo","aripo",
                                    "aripo","aripo","aripo","aripo","aripo"),
                         stringsAsFactors = TRUE)

# Calculate differential methylation
DMSmyDiff5X.dev <- calculateDiffMeth(DMSmeth5X.dev, mc.cores=2, test="Chisq")
DMSmyDiff5X.st <- calculateDiffMeth(DMSmeth5X.st, mc.cores=2, test="Chisq")
DMSmyDiff5X.all <- calculateDiffMeth(DMSmeth5X.all, mc.cores=2, covariates=covariates.pop, test="Chisq")

# Call significant methylation
DMSdiffMeth5X.dev <- getMethylDiff(DMSmyDiff5X.dev, difference = 20, qvalue = 0.0125)
DMSdiffMeth5X.st <- getMethylDiff(DMSmyDiff5X.st, difference = 20, qvalue = 0.0125)
DMSdiffMeth5X.all <- getMethylDiff(DMSmyDiff5X.all, difference = 20, qvalue = 0.0125)

# Check number of significant DMS
DMSdiffMeth5X.dev
DMSdiffMeth5X.st
DMSdiffMeth5X.all

# Get meth per chromosome
DMSdiffMethChr5X.dev <- diffMethPerChr(DMSmyDiff5X.dev, plot=FALSE, qvalue.cutoff=0.0125, meth.cutoff=20)
DMSdiffMethChr5X.st <- diffMethPerChr(DMSmyDiff5X.st, plot=FALSE, qvalue.cutoff=0.0125, meth.cutoff=20)
DMSdiffMethChr5X.all <- diffMethPerChr(DMSmyDiff5X.all, plot=FALSE, qvalue.cutoff=0.0125, meth.cutoff=20)

DMSdiffMethChr5X.dev
DMSdiffMethChr5X.st
DMSdiffMethChr5X.all

## Save R objects ##
saveRDS(myobj5X.dev.subset, file = "./DMS_res/myObj5X_dev.RDS")
saveRDS(DMSmeth5X.dev, file = "./DMS_res/DMSmeth5x_dev.RDS")
saveRDS(DMSmyDiff5X.dev, file = "./DMS_res/DMSmydiff5X_dev.RDS")
saveRDS(DMSdiffMeth5X.dev, file = "./DMS_res/DMSdiffmeth5X_dev.RDS")
saveRDS(DMSdiffMethChr5X.dev, file = "./DMS_res/DMSdiffmethchr5x_dev.RDS")
saveRDS(getData(DMSmeth5X.dev), file = "./DMS_res/DMSmeth5X_dev_data.RDS")
saveRDS(getData(DMSdiffMeth5X.dev), file = "./DMS_res/DMSdiffmeth5X_dev_data.RDS")
saveRDS(getData(DMSmyDiff5X.dev), file = "./DMS_res/DMSmydiff5X_dev_data.RDS")

saveRDS(myobj5X.st.subset, file = "./DMS_res/myObj5X_st.RDS")
saveRDS(DMSmeth5X.st, file = "./DMS_res/DMSmeth5x_st.RDS")
saveRDS(DMSmyDiff5X.st, file = "./DMS_res/DMSmydiff5X_st.RDS")
saveRDS(DMSdiffMeth5X.st, file = "./DMS_res/DMSdiffmeth5X_st.RDS")
saveRDS(DMSdiffMethChr5X.st, file = "./DMS_res/DMSdiffmethchr5x_st.RDS")
saveRDS(getData(DMSmeth5X.st), file = "./DMS_res/DMSmeth5X_st_data.RDS")
saveRDS(getData(DMSdiffMeth5X.st), file = "./DMS_res/DMSdiffmeth5X_st_data.RDS")
saveRDS(getData(DMSmyDiff5X.st), file = "./DMS_res/DMSmydiff5X_st_data.RDS")

saveRDS(myobj5X.all.subset, file = "./DMS_res/myObj5X_all.RDS")
saveRDS(DMSmeth5X.all, file = "./DMS_res/DMSmeth5x_all.RDS")
saveRDS(DMSmyDiff5X.all, file = "./DMS_res/DMSmydiff5X_all.RDS")
saveRDS(DMSdiffMeth5X.all, file = "./DMS_res/DMSdiffmeth5X_all.RDS")
saveRDS(DMSdiffMethChr5X.all, file = "./DMS_res/DMSdiffmethchr5x_all.RDS")
saveRDS(getData(DMSmeth5X.all), file = "./DMS_res/DMSmeth5X_all_data.RDS")
saveRDS(getData(DMSdiffMeth5X.all), file = "./DMS_res/DMSdiffmeth5X_all_data.RDS")
saveRDS(getData(DMSmyDiff5X.all), file = "./DMS_res/DMSmydiff5X_all_data.RDS")

## Find DMRs ##
#tile into 100 bp windows with min coverage 10X and 5X
tiles.5X.dev <- tileMethylCounts(myobj3X.dev.subset, win.size = 100, step.size = 100, cov.bases = 5)
tiles.5X.st <- tileMethylCounts(myobj3X.st.subset, win.size = 100, step.size = 100, cov.bases = 5)
tiles.5X.all <- tileMethylCounts(myobj3X.all.subset, win.size = 100, step.size = 100, cov.bases = 5)

#check number of tiles 
tiles.5X.dev
tiles.5X.st
tiles.5X.all

#unite calls for 60% of samples
DMRmeth5X.dev <- unite(tiles.5X.dev, save.db = FALSE)
DMRmeth5X.st <- unite(tiles.5X.st, save.db = FALSE)
DMRmeth5X.all <- unite(tiles.5X.all, save.db = FALSE)

#check number of regions retained 
DMRmeth5X.dev
DMRmeth5X.st
DMRmeth5X.all

#calculate differential methylation 
DMRmyDiff5X.dev <- calculateDiffMeth(DMRmeth5X.dev, mc.cores=2, test="Chisq")
DMRmyDiff5X.st <- calculateDiffMeth(DMRmeth5X.st, mc.cores=2, test="Chisq")
DMRmyDiff5X.all <- calculateDiffMeth(DMRmeth5X.all, mc.cores=2, test="Chisq", covariates=covariates.pop)

#call significant methylation
DMRdiffMeth5X.dev <- getMethylDiff(DMRmyDiff5X.dev, difference = 20, qvalue = 0.0125)
DMRdiffMeth5X.st <- getMethylDiff(DMRmyDiff5X.st, difference = 20, qvalue = 0.0125)
DMRdiffMeth5X.all <- getMethylDiff(DMRmyDiff5X.all, difference = 20, qvalue = 0.0125)

#check number of DMRs
DMRdiffMeth5X.dev
DMRdiffMeth5X.st
DMRdiffMeth5X.all

#get meth per chromosome
DMRdiffMethChr5X.dev <- diffMethPerChr(DMRmyDiff5X.dev, plot=FALSE,qvalue.cutoff=0.0125, meth.cutoff=20)
DMRdiffMethChr5X.st <- diffMethPerChr(DMRmyDiff5X.st, plot=FALSE,qvalue.cutoff=0.0125, meth.cutoff=20)
DMRdiffMethChr5X.all <- diffMethPerChr(DMRmyDiff5X.all, plot=FALSE,qvalue.cutoff=0.0125, meth.cutoff=20)

DMRdiffMethChr5X.dev
DMRdiffMethChr5X.st
DMRdiffMethChr5X.all

## Save R objects ##
saveRDS(tiles.5X.dev, file = "./DMR_res/DMRtiles5X_dev.RDS")
saveRDS(DMRmeth5X.dev, file = "./DMR_res/DMRmeth5X_dev.RDS")
saveRDS(DMRmyDiff5X.dev, file = "./DMR_res/DMRmydiff5X_dev.RDS")
saveRDS(DMRdiffMeth5X.dev, file = "./DMR_res/DMRdiffmeth5X_dev.RDS")
saveRDS(DMRdiffMethChr5X.dev, file = "./DMR_res/DMRdiffmethchr5X_dev.RDS")
saveRDS(getData(DMRmeth5X.dev), file = "./DMR_res/DMRmeth5X_data_dev.RDS")
saveRDS(getData(DMRmyDiff5X.dev), file = "./DMR_res/DMRmydiff5X_data_dev.RDS")
saveRDS(getData(DMRdiffMeth5X.dev), file = "./DMR_res/DMRdiffmeth5X_data_dev.RDS")

saveRDS(tiles.5X.st, file = "./DMR_res/DMRtiles5X_st.RDS")
saveRDS(DMRmeth5X.st, file = "./DMR_res/DMRmeth5X_st.RDS")
saveRDS(DMRmyDiff5X.st, file = "./DMR_res/DMRmydiff5X_st.RDS")
saveRDS(DMRdiffMeth5X.st, file = "./DMR_res/DMRdiffmeth5X_st.RDS")
saveRDS(DMRdiffMethChr5X.st, file = "./DMR_res/DMRdiffmethchr5X_st.RDS")
saveRDS(getData(DMRmeth5X.st), file = "./DMR_res/DMRmeth5X_data_st.RDS")
saveRDS(getData(DMRmyDiff5X.st), file = "./DMR_res/DMRmydiff5X_data_st.RDS")
saveRDS(getData(DMRdiffMeth5X.st), file = "./DMR_res/DMRdiffmeth5X_data_st.RDS")

saveRDS(tiles.5X.all, file = "./DMR_res/DMRtiles5X_all.RDS")
saveRDS(DMRmeth5X.all, file = "./DMR_res/DMRmeth5X_all.RDS")
saveRDS(DMRmyDiff5X.all, file = "./DMR_res/DMRmydiff5X_all.RDS")
saveRDS(DMRdiffMeth5X.all, file = "./DMR_res/DMRdiffmeth5X_all.RDS")
saveRDS(DMRdiffMethChr5X.all, file = "./DMR_res/DMRdiffmethchr5X_all.RDS")
saveRDS(getData(DMRmeth5X.all), file = "./DMR_res/DMRmeth5X_data_all.RDS")
saveRDS(getData(DMRmyDiff5X.all), file = "./DMR_res/DMRmydiff5X_data_all.RDS")
saveRDS(getData(DMRdiffMeth5X.all), file = "./DMR_res/DMRdiffmeth5X_data_all.RDS")

## Annotation ##
#load in data 
ref.anno <- readTranscriptFeatures("/scratch/janayfox/guppyWGBS/annotation/guppy.sorted.bed", unique.prom = FALSE, up.flank = 1500, down.flank = 500)

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

DMRdiffMeth5X.dev.gr.rename <- renameChr.noXY(DMRdiffMeth5X.dev)
DMRdiffMeth5X.st.gr.rename <- renameChr.noXY(DMRdiffMeth5X.st)
DMRdiffMeth5X.all.gr.rename <- renameChr.noXY(DMRdiffMeth5X.all)

DMRdiffMeth5X.dev.gr.rename <- renameChr.noXY(DMSdiffMeth5X.dev)
DMRdiffMeth5X.st.gr.rename <- renameChr.noXY(DMSdiffMeth5X.st)
DMRdiffMeth5X.all.gr.rename <- renameChr.noXY(DMSdiffMeth5X.all)

CpG.dev.gr.rename <- renameChr.noXY(DMSmyDiff5X.dev)
CpG.st.gr.rename <- renameChr.noXY(DMSmyDiff5X.st)
CpG.all.gr.rename <- renameChr.noXY(DMSmyDiff5X.all)

#make function for annotations
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

anno.func(DMRdiffMeth5X.dev.gr.rename, DMS.diffmeth.dev.gr.rename, CpG.dev.gr.rename,
          "./anno_res/DMR_anno_dev.RDS", "./anno_res/DMS_anno_dev.RDS", "./anno_res/CpG_anno_dev.RDS",
          "./anno_res/DMR_annStats_perc_dev.RDS", "./anno_res/DMS_annStats_perc_dev.RDS", "./anno_res/CpG_annStats_perc_dev.RDS",
          "./anno_res/DMR_annStats_num_dev.RDS", "./anno_res/DMS_annStats_num_dev.RDS", "./anno_res/CpG_annStats_num_dev.RDS",
          "./anno_res/DMR_TSS_dev.RDS", "./0anno_res/DMS_TSS_dev.RDS")

anno.func(DMRdiffMeth5X.st.gr.rename, DMS.diffmeth.st.gr.rename, CpG.st.gr.rename,
          "./anno_res/DMR_anno_st.RDS", "./anno_res/DMS_anno_st.RDS", "./anno_res/CpG_anno_st.RDS",
          "./anno_res/DMR_annStats_perc_st.RDS", "./anno_res/DMS_annStats_perc_st.RDS", "./anno_res/CpG_annStats_perc_st.RDS",
          "./anno_res/DMR_annStats_num_st.RDS", "./anno_res/DMS_annStats_num_st.RDS", "./anno_res/CpG_annStats_num_st.RDS",
          "./anno_res/DMR_TSS_st.RDS", "./0anno_res/DMS_TSS_st.RDS")

anno.func(DMRdiffMeth5X.all.gr.rename, DMS.diffmeth.all.gr.rename, CpG.all.gr.rename,
          "./anno_res/DMR_anno_all.RDS", "./anno_res/DMS_anno_all.RDS", "./anno_res/CpG_anno_all.RDS",
          "./anno_res/DMR_annStats_perc_all.RDS", "./anno_res/DMS_annStats_perc_all.RDS", "./anno_res/CpG_annStats_perc_all.RDS",
          "./anno_res/DMR_annStats_num_all.RDS", "./anno_res/DMS_annStats_num_all.RDS", "./anno_res/CpG_annStats_num_all.RDS",
          "./anno_res/DMR_TSS_all.RDS", "./0anno_res/DMS_TSS_all.RDS")