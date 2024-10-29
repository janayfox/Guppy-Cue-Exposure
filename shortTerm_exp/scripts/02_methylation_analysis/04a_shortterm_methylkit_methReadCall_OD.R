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

setwd("/scratch/janayfox/guppyWGBS_final/methylKit/st/perc20_od")

#Prepare function for getting DMS and DMRs 
run.DMS.DMR <- function(myobj,cov,chr.to.keep,unite.val,filename.myobj3x,filename.myobj5x,
                        file.name.meth.DMS,file.name.meth.DMS.data,
                        file.name.myDiff.DMS,file.name.myDiff.data,
                        file.name.diffMeth.DMS,filename.diffMeth.DMS.data,
                        filename.chr.DMS,
                        file.name.meth.DMR,file.name.meth.DMR.data,file.name.tile.meth,
                        file.name.myDiff.DMR,file.name.myDiff.DMR.data,
                        file.name.diffMeth.DMR,file.name.diffMeth.DMR.data,
                        filename.chr.DMR){

      #filter out sites in the 99.9th percentile of coverage (PCR bias) for 3x and 5X
      myobj.3X=filterByCoverage(myobj,lo.count=3,lo.perc=NULL,
                                    hi.count=NULL, hi.perc=99.9)


      myobj.5X=filterByCoverage(myobj,lo.count=5,lo.perc=NULL,
                                    hi.count=NULL, hi.perc=99.9)

      #normalize by median coverage
      norm.myobj.3X=normalizeCoverage(myobj.3X, method="median")
      norm.myobj.5X=normalizeCoverage(myobj.5X, method="median")

      #remove sex chr (LG12) and unplaced scaffolds
      myobj.subset.3X <- selectByOverlap(norm.myobj.3X, chr.to.keep)
      myobj.subset.5X <- selectByOverlap(norm.myobj.5X, chr.to.keep)

      ##Find DMS##
      #unite sites 
      DMS.meth.5X=unite(myobj.subset.5X, min.per.group = unite.val, destrand=FALSE, save.db = FALSE)

      #convert to non DB object 
      DMS.meth.5X <- as(DMS.meth.5X, "methylBase")

      #filter out low variation sites 
      pm.5X <- percMethylation(DMS.meth.5X) #get percent methylation matrix
      sds.5X <- matrixStats::rowSds(pm.5X, na.rm = TRUE) #calculate standard deviation of CpGs 
      DMS.meth.5X <- DMS.meth.5X[sds.5X > 2]

      #filter out SNPs
      snp <- read.csv("../../../BS-SNPer/shortterm_CT_SNP_edit.csv") #read in snps
      snp.granges <- makeGRangesFromDataFrame(snp, ignore.strand = TRUE) #convert to granges 
      DMS.meth.5X <- DMS.meth.5X[!as(DMS.meth.5X, "GRanges") %over% snp.granges, ] #select CpGs that do not overlap

      # Check number of CpGs 
      head(DMS.meth.5X)

      #calculate differential methylation
      DMS.myDiff.5X <- calculateDiffMeth(DMS.meth.5X, overdispersion="MN", 
                                         covariates=cov, mc.cores=2, test="Chisq", save.db = FALSE)

      #call significant methylation
      DMS.diffMeth.5X <- getMethylDiff(DMS.myDiff.5X, difference = 20, qvalue = 0.0125, save.db = FALSE)

      #check number of significant DMS
      head(DMS.diffMeth.5X)

      # Get meth per chromosome
      DMS.diffMeth.5X.chr <- diffMethPerChr(DMS.myDiff.5X, plot=FALSE, qvalue.cutoff=0.0125, meth.cutoff=20, save.db = FALSE)

      ##Find DMRs ##
      #unite sites 
      DMR.meth.3X=unite(myobj.subset.3X, min.per.group = unite.val, destrand=FALSE, save.db = FALSE)

      #convert to non DB object 
      DMR.meth.3X <- as(DMR.meth.3X, "methylBase")

      #filter out low variation sites 
      pm.3x <- percMethylation(DMR.meth.3X) #get percent methylation matrix
      sds.3x <- matrixStats::rowSds(pm.3x, na.rm = TRUE) #calculate standard deviation of CpGs 
      DMR.meth.3X <- DMR.meth.3X[sds.3x > 2]

      #filter out SNPs
      DMR.meth.3X <- DMR.meth.3X[!as(DMR.meth.3X, "GRanges") %over% snp.granges, ] #select CpGs that do not overlap

      # Check number of CpGs 
      head(DMR.meth.3X)

      #tile into 100 bp windows with min coverage 5X
      tile.meth.5X <- tileMethylCounts(DMR.meth.3X, win.size = 100, step.size = 100, cov.bases = 5)

      #check number of regions retained 
      head(tile.meth.5X)

      #calculate differential methylation 
      DMR.myDiff.5X <- calculateDiffMeth(tile.meth.5X, mc.cores=2, overdispersion="MN", 
                                         test="Chisq", covariates=cov, save.db = FALSE)
       
      #call significant methylation
      DMR.diffMeth.5X <- getMethylDiff(DMR.myDiff.5X, difference = 20, qvalue = 0.0125, save.db = FALSE)

      #check number of DMRs 
      head(DMR.diffMeth.5X)

      #get meth per chromosome 
      DMR.diffMeth.5X.chr <- diffMethPerChr(DMR.diffMeth.5X, plot = FALSE,qvalue.cutoff=0.0125, meth.cutoff=20, save.db = FALSE)

      # Save R objects ##
      saveRDS(myobj.subset.3X, file = filename.myobj3x)
      saveRDS(myobj.subset.5X, file = filename.myobj5x)

      saveRDS(DMS.meth.5X, file = file.name.meth.DMS)
      saveRDS(getData(DMS.meth.5X), file = file.name.meth.DMS.data)
      saveRDS(DMS.myDiff.5X, file = file.name.myDiff.DMS)
      saveRDS(getData(DMS.myDiff.5X), file = file.name.myDiff.data)
      saveRDS(DMS.diffMeth.5X, file = file.name.diffMeth.DMS)
      saveRDS(getData(DMS.diffMeth.5X), file = filename.diffMeth.DMS.data)
      saveRDS(DMS.diffMeth.5X.chr, file = filename.chr.DMS)

      
      saveRDS(DMR.meth.3X, file = file.name.meth.DMR)
      saveRDS(getData(DMR.meth.3X), file = file.name.meth.DMR.data)
      saveRDS(tile.meth.5X, file = file.name.tile.meth)
      saveRDS(DMR.myDiff.5X, file = file.name.myDiff.DMR)
      saveRDS(getData(DMR.myDiff.5X), file = file.name.myDiff.DMR.data)
      saveRDS(DMR.diffMeth.5X, file = file.name.diffMeth.DMR)
      saveRDS(getData(DMR.diffMeth.5X), file = file.name.diffMeth.DMR.data)
      saveRDS(DMR.diffMeth.5X.chr, file = filename.chr.DMR)

}

#prepare GRanges object for chromosomes to keep 
#to remove unplaced scaffolds and sex chromosomes
keep.chr.noXY <- GRanges(seqnames = c("NC_024331.1", "NC_024332.1", "NC_024333.1", "NC_024334.1", "NC_024335.1",
                                    "NC_024336.1", "NC_024337.1", "NC_024338.1", "NC_024339.1", "NC_024340.1",
                                    "NC_024341.1", "NC_024343.1", "NC_024344.1", "NC_024345.1", "NC_024346.1",
                                    "NC_024347.1", "NC_024348.1", "NC_024349.1", "NC_024350.1", "NC_024351.1",
                                    "NC_024352.1", "NC_024353.1"),
                        ranges=IRanges(start = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
                                    end = c(34115677,46286544,35265442,31497199,33908744,
                                                31529174,31413364,27946405,34117797,32819797,
                                                28875558,33524197,28338960,30644713,33199742,
                                                30788009,22026651,28470737,26385442,25773841,
                                                25248790,18084596)),
                        strand="*")

seqlengths(keep.chr.noXY)=c(34115677,46286544,35265442,31497199,33908744,
                        31529174,31413364,27946405,34117797,32819797,
                        28875558,33524197,28338960,30644713,33199742,
                        30788009,22026651,28470737,26385442,25773841,
                        25248790,18084596)

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


#0.5h
#create lists of file locations
file.list.05h = list("../../../merged_sequences/st/ST2AC10F.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2AC10M.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2AC16F.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2AC16M.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2AC3F.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2AC3M.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2C10F.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2C10M.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2C13F.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2C13M.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2C2F.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2C2M.CpG_merged.cov")

file.list.05h.fem = list("../../../merged_sequences/st/ST2AC10F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC16F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC3F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C10F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C13F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C2F.CpG_merged.cov")

file.list.05h.mal = list("../../../merged_sequences/st/ST2AC10M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC16M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC3M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C10M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C13M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C2M.CpG_merged.cov")

#create tabix file
myobj.05h=methRead(file.list.05h,
                   sample.id=list("ST2AC10F","ST2AC10M","ST2AC16F","ST2AC16M","ST2AC3F","ST2AC3M",
                                  "ST2C10F","ST2C10M","ST2C13F","ST2C13M","ST2C2F","ST2C2M"),
                   assembly="guppyWGBS_shortterm",
                   pipeline="bismarkCoverage",
                   treatment=c(1,1,1,1,1,1,
                               0,0,0,0,0,0),
                   context="CpG",
                   mincov = 1,
                   dbtype = "tabix",
                   dbdir = "shortterm_05h_DB_merged"
)

myobj.05h.fem=methRead(file.list.05h.fem,
                       sample.id=list("ST2AC10F","ST2AC16F","ST2AC3F",
                                      "ST2C10F","ST2C13F","ST2C2F"),
                       assembly="guppyWGBS_shortterm",
                       pipeline="bismarkCoverage",
                       treatment=c(1,1,1,
                                   0,0,0),
                       context="CpG",
                       mincov = 1,
                       dbtype = "tabix",
                       dbdir = "shortterm_05hF_DB_merged"
)

myobj.05h.mal=methRead(file.list.05h.mal,
                       sample.id=list("ST2AC10M","ST2AC16M","ST2AC3M",
                                      "ST2C10M","ST2C13M","ST2C2M"),
                       assembly="guppyWGBS_shortterm",
                       pipeline="bismarkCoverage",
                       treatment=c(1,1,1,
                                   0,0,0),
                       context="CpG",
                       mincov = 1,
                       dbtype = "tabix",
                       dbdir = "shortterm_05hM_DB_merged"
)

#enter covariates
covariates.05h <- data.frame(tank=c("AC10","AC10","AC16","AC16","AC3","AC3",
                                    "C10","C10","C13","C13","C2","C2"),
                             sex=c("F","M","F","M","F","M",
                                   "F","M","F","M","F","M"),
                             stringsAsFactors = TRUE)

run.DMS.DMR(myobj.05h,covariates.05h,keep.chr.noXY,6L,
                        "./05h/DMS_res_05h/myobj_05h_3X.RDS","./05h/DMS_res_05h/myobj_05h_5X.RDS",
                        "./05h/DMS_res_05h/DMSmeth_05h_5X.RDS","./05h/DMS_res_05h/meth_05h_5X_data.RDS",
                        "./05h/DMS_res_05h/DMSmydiff_05h_5X.RDS","./05h/DMS_res_05h/DMSmydiff_05h_5X_data.RDS",
                        "./05h/DMS_res_05h/DMSdiffmeth_05h_5X.RDS","./05h/DMS_res_05h/DMSdiffmeth_05h_5X_data.RDS",
                        "./05h/DMS_res_05h/DMSdiffmethchr_05h_5X.RDS",
                        "./05h/DMR_res_05h/DMRmeth_05h_5X.RDS","./05h/DMR_res_05h/DMRmeth_05h_5X_data.RDS",
                        "./05h/DMR_res_05h/DMR_tile_meth_05h_5X.RDS",
                        "./05h/DMR_res_05h/DMRmydiff_05h_5X.RDS","./05h/DMR_res_05h/DMRmydiff_05h_5X_data.RDS",
                        "./05h/DMR_res_05h/DMRdiffmeth_05h_5X.RDS","./05h/DMR_res_05h/DMRdiffmeth_05h_5X_data.RDS",
                        "./05h/DMR_res_05h/DMRdiffmethchr_05h_5X.RDS")

run.DMS.DMR(myobj.05h.fem,NULL,keep.chr.allchr,3L,
                        "./05h/DMS_res_05h/myobj_05h_fem_3X.RDS","./05h/DMS_res_05h/myobj_05h_fem_5X.RDS",
                        "./05h/DMS_res_05h/DMSmeth_05h_fem_5X.RDS","./05h/DMS_res_05h/meth_05h_fem_5X_data.RDS",
                        "./05h/DMS_res_05h/DMSmydiff_05h_fem_5X.RDS","./05h/DMS_res_05h/DMSmydiff_05h_fem_5X_data.RDS",
                        "./05h/DMS_res_05h/DMSdiffmeth_05h_fem_5X.RDS","./05h/DMS_res_05h/DMSdiffmeth_05h_fem_5X_data.RDS",
                        "./05h/DMS_res_05h/DMSdiffmethchr_05h_fem_5X.RDS",
                        "./05h/DMR_res_05h/DMRmeth_05h_fem_5X.RDS","./05h/DMR_res_05h/DMRmeth_05h_fem_5X_data.RDS",
                        "./05h/DMR_res_05h/DMR_tile_meth_05h_fem_5X.RDS",
                        "./05h/DMR_res_05h/DMRmydiff_05h_fem_5X.RDS","./05h/DMR_res_05h/DMRmydiff_05h_fem_5X_data.RDS",
                        "./05h/DMR_res_05h/DMRdiffmeth_05h_fem_5X.RDS","./05h/DMR_res_05h/DMRdiffmeth_05h_fem_5X_data.RDS",
                        "./05h/DMR_res_05h/DMRdiffmethchr_05h_fem_5X.RDS")

run.DMS.DMR(myobj.05h.mal,NULL,keep.chr.allchr,3L,
                        "./05h/DMS_res_05h/myobj_05h_mal_3X.RDS","./05h/DMS_res_05h/myobj_05h_mal_5X.RDS",
                        "./05h/DMS_res_05h/DMSmeth_05h_mal_5X.RDS","./05h/DMS_res_05h/meth_05h_mal_5X_data.RDS",
                        "./05h/DMS_res_05h/DMSmydiff_05h_mal_5X.RDS","./05h/DMS_res_05h/DMSmydiff_05h_mal_5X_data.RDS",
                        "./05h/DMS_res_05h/DMSdiffmeth_05h_mal_5X.RDS","./05h/DMS_res_05h/DMSdiffmeth_05h_mal_5X_data.RDS",
                        "./05h/DMS_res_05h/DMSdiffmethchr_05h_mal_5X.RDS",
                        "./05h/DMR_res_05h/DMRmeth_05h_mal_5X.RDS","./05h/DMR_res_05h/DMRmeth_05h_mal_5X_data.RDS",
                        "./05h/DMR_res_05h/DMR_tile_meth_05h_mal_5X.RDS",
                        "./05h/DMR_res_05h/DMRmydiff_05h_mal_5X.RDS","./05h/DMR_res_05h/DMRmydiff_05h_mal_5X_data.RDS",
                        "./05h/DMR_res_05h/DMRdiffmeth_05h_mal_5X.RDS","./05h/DMR_res_05h/DMRdiffmeth_05h_mal_5X_data.RDS",
                        "./05h/DMR_res_05h/DMRdiffmethchr_05h_mal_5X.RDS")

rm(myobj.05h,myobj.05h.fem,myobj.05h.mal)

## Prepare tabix files
#run 1h
#create lists of file locations
file.list.1h = list("../../../merged_sequences/st/ST2AC14F.CpG_merged.cov",
                    "../../../merged_sequences/st/ST2AC14M.CpG_merged.cov",
                    "../../../merged_sequences/st/ST2AC2F.CpG_merged.cov",
                    "../../../merged_sequences/st/ST2AC2M.CpG_merged.cov",
                    "../../../merged_sequences/st/ST2AC9F.CpG_merged.cov",
                    "../../../merged_sequences/st/ST2AC9M.CpG_merged.cov",
                    "../../../merged_sequences/st/ST2C14F.CpG_merged.cov",
                    "../../../merged_sequences/st/ST2C14M.CpG_merged.cov",
                    "../../../merged_sequences/st/ST2C3F.CpG_merged.cov",
                    "../../../merged_sequences/st/ST2C3M.CpG_merged.cov",
                    "../../../merged_sequences/st/ST2C9F.CpG_merged.cov",
                    "../../../merged_sequences/st/ST2C9M.CpG_merged.cov")

file.list.1h.fem = list("../../../merged_sequences/st/ST2AC14F.CpG_merged.cov",
                        "../../../merged_sequences/st/ST2AC2F.CpG_merged.cov",
                        "../../../merged_sequences/st/ST2AC9F.CpG_merged.cov",
                        "../../../merged_sequences/st/ST2C14F.CpG_merged.cov",
                        "../../../merged_sequences/st/ST2C3F.CpG_merged.cov",
                        "../../../merged_sequences/st/ST2C9F.CpG_merged.cov")

file.list.1h.mal = list("../../../merged_sequences/st/ST2AC14M.CpG_merged.cov",
                        "../../../merged_sequences/st/ST2AC2M.CpG_merged.cov",
                        "../../../merged_sequences/st/ST2AC9M.CpG_merged.cov",
                        "../../../merged_sequences/st/ST2C14M.CpG_merged.cov",
                        "../../../merged_sequences/st/ST2C3M.CpG_merged.cov",
                        "../../../merged_sequences/st/ST2C9M.CpG_merged.cov")

#create tabix file
myobj.1h=methRead(file.list.1h,
                  sample.id=list("ST2AC14F","ST2AC14M","ST2AC2F","ST2AC2M","ST2AC9F","ST2AC9M",
                                 "ST2C14F","ST2C14M","ST2C3F","ST2C3M","ST2C9F", "ST2C9M"),
                  assembly="guppyWGBS_shortterm",
                  pipeline="bismarkCoverage",
                  treatment=c(1,1,1,1,1,1,
                              0,0,0,0,0,0),
                  context="CpG",
                  mincov = 1,
                  dbtype = "tabix",
                  dbdir = "shortterm_1h_DB_merged"
)

myobj.1h.fem=methRead(file.list.1h.fem,
                      sample.id=list("ST2AC14F","ST2AC2F","ST2AC9F",
                                     "ST2C14F","ST2C3F","ST2C9F"),
                      assembly="guppyWGBS_shortterm",
                      pipeline="bismarkCoverage",
                      treatment=c(1,1,1,
                                  0,0,0),
                      context="CpG",
                      mincov = 1,
                      dbtype = "tabix",
                      dbdir = "shortterm_1hF_DB_merged"
)

myobj.1h.mal=methRead(file.list.1h.mal,
                      sample.id=list("ST2AC14M","ST2AC2M","ST2AC9M",
                                     "ST2C14M","ST2C3M","ST2C9M"),
                      assembly="guppyWGBS_shortterm",
                      pipeline="bismarkCoverage",
                      treatment=c(1,1,1,
                                  0,0,0),
                      context="CpG",
                      mincov = 1,
                      dbtype = "tabix",
                      dbdir = "shortterm_1hM_DB_merged"
)

#enter covariates
covariates.1h <- data.frame(tank=c("AC14","AC14","AC2","AC2","AC9","AC9",
                                   "C14","C14","C3","C3","C9","C9"),
                            sex=c("F","M","F","M","F","M",
                                  "F","M","F","M","F","M"),
                            stringsAsFactors = TRUE)

run.DMS.DMR(myobj.1h,covariates.1h,keep.chr.noXY,6L,
                        "./1h/DMS_res_1h/myobj_1h_3X.RDS","./1h/DMS_res_1h/myobj_1h_5X.RDS",
                        "./1h/DMS_res_1h/DMSmeth_1h_5X.RDS","./1h/DMS_res_1h/meth_1h_5X_data.RDS",
                        "./1h/DMS_res_1h/DMSmydiff_1h_5X.RDS","./1h/DMS_res_1h/DMSmydiff_1h_5X_data.RDS",
                        "./1h/DMS_res_1h/DMSdiffmeth_1h_5X.RDS","./1h/DMS_res_1h/DMSdiffmeth_1h_5X_data.RDS",
                        "./1h/DMS_res_1h/DMSdiffmethchr_1h_5X.RDS",
                        "./1h/DMR_res_1h/DMRmeth_1h_5X.RDS","./1h/DMR_res_1h/DMRmeth_1h_5X_data.RDS",
                        "./1h/DMR_res_1h/DMR_tile_meth_1h_5X.RDS",
                        "./1h/DMR_res_1h/DMRmydiff_1h_5X.RDS","./1h/DMR_res_1h/DMRmydiff_1h_5X_data.RDS",
                        "./1h/DMR_res_1h/DMRdiffmeth_1h_5X.RDS","./1h/DMR_res_1h/DMRdiffmeth_1h_5X_data.RDS",
                        "./1h/DMR_res_1h/DMRdiffmethchr_1h_5X.RDS")

run.DMS.DMR(myobj.1h.fem,NULL,keep.chr.allchr,3L,
                        "./1h/DMS_res_1h/myobj_1h_fem_3X.RDS","./1h/DMS_res_1h/myobj_1h_fem_5X.RDS",
                        "./1h/DMS_res_1h/DMSmeth_1h_fem_5X.RDS","./1h/DMS_res_1h/meth_1h_fem_5X_data.RDS",
                        "./1h/DMS_res_1h/DMSmydiff_1h_fem_5X.RDS","./1h/DMS_res_1h/DMSmydiff_1h_fem_5X_data.RDS",
                        "./1h/DMS_res_1h/DMSdiffmeth_1h_fem_5X.RDS","./1h/DMS_res_1h/DMSdiffmeth_1h_fem_5X_data.RDS",
                        "./1h/DMS_res_1h/DMSdiffmethchr_1h_fem_5X.RDS",
                        "./1h/DMR_res_1h/DMRmeth_1h_fem_5X.RDS","./1h/DMR_res_1h/DMRmeth_1h_fem_5X_data.RDS",
                        "./1h/DMR_res_1h/DMR_tile_meth_1h_fem_5X.RDS",
                        "./1h/DMR_res_1h/DMRmydiff_1h_fem_5X.RDS","./1h/DMR_res_1h/DMRmydiff_1h_fem_5X_data.RDS",
                        "./1h/DMR_res_1h/DMRdiffmeth_1h_fem_5X.RDS","./1h/DMR_res_1h/DMRdiffmeth_1h_fem_5X_data.RDS",
                        "./1h/DMR_res_1h/DMRdiffmethchr_1h_fem_5X.RDS")

run.DMS.DMR(myobj.1h.mal,NULL,keep.chr.allchr,3L,
                        "./1h/DMS_res_1h/myobj_1h_mal_3X.RDS","./1h/DMS_res_1h/myobj_1h_mal_5X.RDS",
                        "./1h/DMS_res_1h/DMSmeth_1h_mal_5X.RDS","./1h/DMS_res_1h/meth_1h_mal_5X_data.RDS",
                        "./1h/DMS_res_1h/DMSmydiff_1h_mal_5X.RDS","./1h/DMS_res_1h/DMSmydiff_1h_mal_5X_data.RDS",
                        "./1h/DMS_res_1h/DMSdiffmeth_1h_mal_5X.RDS","./1h/DMS_res_1h/DMSdiffmeth_1h_mal_5X_data.RDS",
                        "./1h/DMS_res_1h/DMSdiffmethchr_1h_mal_5X.RDS",
                        "./1h/DMR_res_1h/DMRmeth_1h_mal_5X.RDS","./1h/DMR_res_1h/DMRmeth_1h_mal_5X_data.RDS",
                        "./1h/DMR_res_1h/DMR_tile_meth_1h_mal_5X.RDS",
                        "./1h/DMR_res_1h/DMRmydiff_1h_mal_5X.RDS","./1h/DMR_res_1h/DMRmydiff_1h_mal_5X_data.RDS",
                        "./1h/DMR_res_1h/DMRdiffmeth_1h_mal_5X.RDS","./1h/DMR_res_1h/DMRdiffmeth_1h_mal_5X_data.RDS",
                        "./1h/DMR_res_1h/DMRdiffmethchr_1h_mal_5X.RDS")

rm(myobj.1h,myobj.1h.fem,myobj.1h.mal)

#run 4h
#create lists of file locations
file.list.4h = list("../../../merged_sequences/st/ST2AC13F.CpG_merged.cov",
                    "../../../merged_sequences/st/ST2AC13M.CpG_merged.cov",
                    "../../../merged_sequences/st/ST2AC1F.CpG_merged.cov",
                    "../../../merged_sequences/st/ST2AC1M.CpG_merged.cov",
                    "../../../merged_sequences/st/ST2AC7F.CpG_merged.cov",
                    "../../../merged_sequences/st/ST2AC7M.CpG_merged.cov",
                    "../../../merged_sequences/st/ST2C12F.CpG_merged.cov",
                    "../../../merged_sequences/st/ST2C12M.CpG_merged.cov",
                    "../../../merged_sequences/st/ST2C1F.CpG_merged.cov",
                    "../../../merged_sequences/st/ST2C1M.CpG_merged.cov",
                    "../../../merged_sequences/st/ST2C7F.CpG_merged.cov",
                    "../../../merged_sequences/st/ST2C7M.CpG_merged.cov")

file.list.4h.fem = list("../../../merged_sequences/st/ST2AC13F.CpG_merged.cov",
                        "../../../merged_sequences/st/ST2AC1F.CpG_merged.cov",
                        "../../../merged_sequences/st/ST2AC7F.CpG_merged.cov",
                        "../../../merged_sequences/st/ST2C12F.CpG_merged.cov",
                        "../../../merged_sequences/st/ST2C1F.CpG_merged.cov",
                        "../../../merged_sequences/st/ST2C7F.CpG_merged.cov")

file.list.4h.mal = list("../../../merged_sequences/st/ST2AC13M.CpG_merged.cov",
                        "../../../merged_sequences/st/ST2AC1M.CpG_merged.cov",
                        "../../../merged_sequences/st/ST2AC7M.CpG_merged.cov",
                        "../../../merged_sequences/st/ST2C12M.CpG_merged.cov",
                        "../../../merged_sequences/st/ST2C1M.CpG_merged.cov",
                        "../../../merged_sequences/st/ST2C7M.CpG_merged.cov")

#create tabix file
myobj.4h=methRead(file.list.4h,
                  sample.id=list("ST2AC13F","ST2AC13M","ST2AC1F","ST2AC1M","ST2AC7F","ST2AC7M",
                                 "ST2C12F","ST2C12M","ST2C1F","ST2C1M","ST2C7F","ST2C7M"),
                  assembly="guppyWGBS_shortterm",
                  pipeline="bismarkCoverage",
                  treatment=c(1,1,1,1,1,1,
                              0,0,0,0,0,0),
                  context="CpG",
                  mincov = 1,
                  dbtype = "tabix",
                  dbdir = "shortterm_4h_DB"
)

myobj.4h.fem=methRead(file.list.4h.fem,
                      sample.id=list("ST2AC13F","ST2AC1F","ST2AC7F",
                                     "ST2C12F","ST2C1F","ST2C7F"),
                      assembly="guppyWGBS_shortterm",
                      pipeline="bismarkCoverage",
                      treatment=c(1,1,1,
                                  0,0,0),
                      context="CpG",
                      mincov = 1,
                      dbtype = "tabix",
                      dbdir = "shortterm_4hF_DB"
)

myobj.4h.mal=methRead(file.list.4h.mal,
                      sample.id=list("ST2AC13M","ST2AC1M","ST2AC7M",
                                     "ST2C12M","ST2C1M","ST2C7M"),
                      assembly="guppyWGBS_shortterm",
                      pipeline="bismarkCoverage",
                      treatment=c(1,1,1,
                                  0,0,0),
                      context="CpG",
                      mincov = 1,
                      dbtype = "tabix",
                      dbdir = "shortterm_4hmal_DB"
)

covariates.4h <- data.frame(tank=c("AC13","AC13","AC1","AC1","AC7","AC7",
                                   "C12","C12","C1","C1","C7","C7"),
                            sex=c("F","M","F","M","F","M",
                                  "F","M","F","M","F","M"),
                            stringsAsFactors = TRUE)

run.DMS.DMR(myobj.4h,covariates.4h,keep.chr.noXY,6L,
                        "./4h/DMS_res_4h/myobj_4h_3X.RDS","./4h/DMS_res_4h/myobj_4h_5X.RDS",
                        "./4h/DMS_res_4h/DMSmeth_4h_5X.RDS","./4h/DMS_res_4h/meth_4h_5X_data.RDS",
                        "./4h/DMS_res_4h/DMSmydiff_4h_5X.RDS","./4h/DMS_res_4h/DMSmydiff_4h_5X_data.RDS",
                        "./4h/DMS_res_4h/DMSdiffmeth_4h_5X.RDS","./4h/DMS_res_4h/DMSdiffmeth_4h_5X_data.RDS",
                        "./4h/DMS_res_4h/DMSdiffmethchr_4h_5X.RDS",
                        "./4h/DMR_res_4h/DMRmeth_4h_5X.RDS","./4h/DMR_res_4h/DMRmeth_4h_5X_data.RDS",
                        "./4h/DMR_res_4h/DMR_tile_meth_4h_5X.RDS",
                        "./4h/DMR_res_4h/DMRmydiff_4h_5X.RDS","./4h/DMR_res_4h/DMRmydiff_4h_5X_data.RDS",
                        "./4h/DMR_res_4h/DMRdiffmeth_4h_5X.RDS","./4h/DMR_res_4h/DMRdiffmeth_4h_5X_data.RDS",
                        "./4h/DMR_res_4h/DMRdiffmethchr_4h_5X.RDS")

run.DMS.DMR(myobj.4h.fem,NULL,keep.chr.allchr,3L,
                        "./4h/DMS_res_4h/myobj_4h_fem_3X.RDS","./4h/DMS_res_4h/myobj_4h_fem_5X.RDS",
                        "./4h/DMS_res_4h/DMSmeth_4h_fem_5X.RDS","./4h/DMS_res_4h/meth_4h_fem_5X_data.RDS",
                        "./4h/DMS_res_4h/DMSmydiff_4h_fem_5X.RDS","./4h/DMS_res_4h/DMSmydiff_4h_fem_5X_data.RDS",
                        "./4h/DMS_res_4h/DMSdiffmeth_4h_fem_5X.RDS","./4h/DMS_res_4h/DMSdiffmeth_4h_fem_5X_data.RDS",
                        "./4h/DMS_res_4h/DMSdiffmethchr_4h_fem_5X.RDS",
                        "./4h/DMR_res_4h/DMRmeth_4h_fem_5X.RDS","./4h/DMR_res_4h/DMRmeth_4h_fem_5X_data.RDS",
                        "./4h/DMR_res_4h/DMR_tile_meth_4h_fem_5X.RDS",
                        "./4h/DMR_res_4h/DMRmydiff_4h_fem_5X.RDS","./4h/DMR_res_4h/DMRmydiff_4h_fem_5X_data.RDS",
                        "./4h/DMR_res_4h/DMRdiffmeth_4h_fem_5X.RDS","./4h/DMR_res_4h/DMRdiffmeth_4h_fem_5X_data.RDS",
                        "./4h/DMR_res_4h/DMRdiffmethchr_4h_fem_5X.RDS")

run.DMS.DMR(myobj.4h.mal,NULL,keep.chr.allchr,3L,
                        "./4h/DMS_res_4h/myobj_4h_mal_3X.RDS","./4h/DMS_res_4h/myobj_4h_mal_5X.RDS",
                        "./4h/DMS_res_4h/DMSmeth_4h_mal_5X.RDS","./4h/DMS_res_4h/meth_4h_mal_5X_data.RDS",
                        "./4h/DMS_res_4h/DMSmydiff_4h_mal_5X.RDS","./4h/DMS_res_4h/DMSmydiff_4h_mal_5X_data.RDS",
                        "./4h/DMS_res_4h/DMSdiffmeth_4h_mal_5X.RDS","./4h/DMS_res_4h/DMSdiffmeth_4h_mal_5X_data.RDS",
                        "./4h/DMS_res_4h/DMSdiffmethchr_4h_mal_5X.RDS",
                        "./4h/DMR_res_4h/DMRmeth_4h_mal_5X.RDS","./4h/DMR_res_4h/DMRmeth_4h_mal_5X_data.RDS",
                        "./4h/DMR_res_4h/DMR_tile_meth_4h_mal_5X.RDS",
                        "./4h/DMR_res_4h/DMRmydiff_4h_mal_5X.RDS","./4h/DMR_res_4h/DMRmydiff_4h_mal_5X_data.RDS",
                        "./4h/DMR_res_4h/DMRdiffmeth_4h_mal_5X.RDS","./4h/DMR_res_4h/DMRdiffmeth_4h_mal_5X_data.RDS",
                        "./4h/DMR_res_4h/DMRdiffmethchr_4h_mal_5X.RDS")

rm(myobj.4h,myobj.4h.fem,myobj.4h.mal)

#24h
#create lists of file locations
file.list.24h = list("../../../merged_sequences/st/ST2AC15F.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2AC15M.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2AC4F.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2AC4M.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2AC8F.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2AC8M.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2C15F.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2C15M.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2C4F.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2C4M.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2C8F.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2C8M.CpG_merged.cov")

file.list.24h.fem = list("../../../merged_sequences/st/ST2AC15F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC4F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC8F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C15F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C4F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C8F.CpG_merged.cov")

file.list.24h.mal = list("../../../merged_sequences/st/ST2AC15M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC4M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC8M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C15M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C4M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C8M.CpG_merged.cov")

#create tabix file
myobj.24h=methRead(file.list.24h,
                   sample.id=list("ST2AC15F","ST2AC15M","ST2AC4F","ST2AC4M","ST2AC8F","ST2AC8M",
                                  "ST2C15F","ST2C15M","ST2C4F","ST2C4M","ST2C8F","ST2C8M"),
                   assembly="guppyWGBS_shortterm",
                   pipeline="bismarkCoverage",
                   treatment=c(1,1,1,1,1,1,
                               0,0,0,0,0,0),
                   context="CpG",
                   mincov = 1,
                   dbtype = "tabix",
                   dbdir = "shortterm_24h_DB_merged"
)

myobj.24h.fem=methRead(file.list.24h.fem,
                       sample.id=list("ST2AC15F","ST2AC4F","ST2AC8F",
                                      "ST2C15F","ST2C4F","ST2C8F"),
                       assembly="guppyWGBS_shortterm",
                       pipeline="bismarkCoverage",
                       treatment=c(1,1,1,
                                   0,0,0),
                       context="CpG",
                       mincov = 1,
                       dbtype = "tabix",
                       dbdir = "shortterm_24hF_DB_merged"
)

myobj.24h.mal=methRead(file.list.24h.mal,
                       sample.id=list("ST2AC15M","ST2AC4M","ST2AC8M",
                                      "ST2C15M","ST2C4M","ST2C8M"),
                       assembly="guppyWGBS_shortterm",
                       pipeline="bismarkCoverage",
                       treatment=c(1,1,1,
                                   0,0,0),
                       context="CpG",
                       mincov = 1,
                       dbtype = "tabix",
                       dbdir = "shortterm_24hM_DB_merged"
)

# #enter covariates
covariates.24h <- data.frame(tank=c("AC15","AC15","AC4","AC4","AC8","AC8",
                                    "C15","C15","C4","C4","C8","C8"),
                             sex=c("F","M","F","M","F","M",
                                   "F","M","F","M","F","M"),
                             stringsAsFactors = TRUE)

run.DMS.DMR(myobj.24h,covariates.24h,keep.chr.noXY,6L,
                        "./24h/DMS_res_24h/myobj_24h_3X.RDS","./24h/DMS_res_24h/myobj_24h_5X.RDS",
                        "./24h/DMS_res_24h/DMSmeth_24h_5X.RDS","./24h/DMS_res_24h/meth_24h_5X_data.RDS",
                        "./24h/DMS_res_24h/DMSmydiff_24h_5X.RDS","./24h/DMS_res_24h/DMSmydiff_24h_5X_data.RDS",
                        "./24h/DMS_res_24h/DMSdiffmeth_24h_5X.RDS","./24h/DMS_res_24h/DMSdiffmeth_24h_5X_data.RDS",
                        "./24h/DMS_res_24h/DMSdiffmethchr_24h_5X.RDS",
                        "./24h/DMR_res_24h/DMRmeth_24h_5X.RDS","./24h/DMR_res_24h/DMRmeth_24h_5X_data.RDS",
                        "./24h/DMR_res_24h/DMR_tile_meth_24h_5X.RDS",
                        "./24h/DMR_res_24h/DMRmydiff_24h_5X.RDS","./24h/DMR_res_24h/DMRmydiff_24h_5X_data.RDS",
                        "./24h/DMR_res_24h/DMRdiffmeth_24h_5X.RDS","./24h/DMR_res_24h/DMRdiffmeth_24h_5X_data.RDS",
                        "./24h/DMR_res_24h/DMRdiffmethchr_24h_5X.RDS")

run.DMS.DMR(myobj.24h.fem,NULL,keep.chr.allchr,3L,
                        "./24h/DMS_res_24h/myobj_24h_fem_3X.RDS","./24h/DMS_res_24h/myobj_24h_fem_5X.RDS",
                        "./24h/DMS_res_24h/DMSmeth_24h_fem_5X.RDS","./24h/DMS_res_24h/meth_24h_fem_5X_data.RDS",
                        "./24h/DMS_res_24h/DMSmydiff_24h_fem_5X.RDS","./24h/DMS_res_24h/DMSmydiff_24h_fem_5X_data.RDS",
                        "./24h/DMS_res_24h/DMSdiffmeth_24h_fem_5X.RDS","./24h/DMS_res_24h/DMSdiffmeth_24h_fem_5X_data.RDS",
                        "./24h/DMS_res_24h/DMSdiffmethchr_24h_fem_5X.RDS",
                        "./24h/DMR_res_24h/DMRmeth_24h_fem_5X.RDS","./24h/DMR_res_24h/DMRmeth_24h_fem_5X_data.RDS",
                        "./24h/DMR_res_24h/DMR_tile_meth_24h_fem_5X.RDS",
                        "./24h/DMR_res_24h/DMRmydiff_24h_fem_5X.RDS","./24h/DMR_res_24h/DMRmydiff_24h_fem_5X_data.RDS",
                        "./24h/DMR_res_24h/DMRdiffmeth_24h_fem_5X.RDS","./24h/DMR_res_24h/DMRdiffmeth_24h_fem_5X_data.RDS",
                        "./24h/DMR_res_24h/DMRdiffmethchr_24h_fem_5X.RDS")

run.DMS.DMR(myobj.24h.mal,NULL,keep.chr.allchr,3L,
                        "./24h/DMS_res_24h/myobj_24h_mal_3X.RDS","./24h/DMS_res_24h/myobj_24h_mal_5X.RDS",
                        "./24h/DMS_res_24h/DMSmeth_24h_mal_5X.RDS","./24h/DMS_res_24h/meth_24h_mal_5X_data.RDS",
                        "./24h/DMS_res_24h/DMSmydiff_24h_mal_5X.RDS","./24h/DMS_res_24h/DMSmydiff_24h_mal_5X_data.RDS",
                        "./24h/DMS_res_24h/DMSdiffmeth_24h_mal_5X.RDS","./24h/DMS_res_24h/DMSdiffmeth_24h_mal_5X_data.RDS",
                        "./24h/DMS_res_24h/DMSdiffmethchr_24h_mal_5X.RDS",
                        "./24h/DMR_res_24h/DMRmeth_24h_mal_5X.RDS","./24h/DMR_res_24h/DMRmeth_24h_mal_5X_data.RDS",
                        "./24h/DMR_res_24h/DMR_tile_meth_24h_mal_5X.RDS",
                        "./24h/DMR_res_24h/DMRmydiff_24h_mal_5X.RDS","./24h/DMR_res_24h/DMRmydiff_24h_mal_5X_data.RDS",
                        "./24h/DMR_res_24h/DMRdiffmeth_24h_mal_5X.RDS","./24h/DMR_res_24h/DMRdiffmeth_24h_mal_5X_data.RDS",
                        "./24h/DMR_res_24h/DMRdiffmethchr_24h_mal_5X.RDS")

rm(myobj.24h,myobj.24h.fem,myobj.24h.mal)

#72h
#create lists of file locations
file.list.72h = list("../../../merged_sequences/st/ST2AC11F.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2AC11M.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2AC5F.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2AC5M.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2AC6F.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2AC6M.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2C11F.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2C11M.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2C5F.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2C5M.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2C6F.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2C6M.CpG_merged.cov")

file.list.72h.fem = list("../../../merged_sequences/st/ST2AC11F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC5F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC6F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C11F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C5F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C6F.CpG_merged.cov")

file.list.72h.mal = list("../../../merged_sequences/st/ST2AC11M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC5M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC6M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C11M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C5M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C6M.CpG_merged.cov")

#create tabix file
myobj.72h=methRead(file.list.72h,
                   sample.id=list("ST2AC11F","ST2AC11M","ST2AC5F","ST2AC5M","ST2AC6F","ST2AC6M",
                                  "ST2C11F","ST2C11M","ST2C5F","ST2C5M","ST2C6F","ST2C6M"),
                   assembly="guppyWGBS_shortterm",
                   pipeline="bismarkCoverage",
                   treatment=c(1,1,1,1,1,1,
                               0,0,0,0,0,0),
                   context="CpG",
                   mincov = 1,
                   dbtype = "tabix",
                   dbdir = "shortterm_72h_DB_merged"
)

myobj.72h.fem=methRead(file.list.72h.fem,
                       sample.id=list("ST2AC11F","ST2AC5F","ST2AC6F",
                                      "ST2C11F","ST2C5F","ST2C6F"),
                       assembly="guppyWGBS_shortterm",
                       pipeline="bismarkCoverage",
                       treatment=c(1,1,1,
                                   0,0,0),
                       context="CpG",
                       mincov = 1,
                       dbtype = "tabix",
                       dbdir = "shortterm_72hF_DB_merged"
)

myobj.72h.mal=methRead(file.list.72h.mal,
                       sample.id=list("ST2AC11M","ST2AC5M","ST2AC6M",
                                      "ST2C11M","ST2C5M","ST2C6M"),
                       assembly="guppyWGBS_shortterm",
                       pipeline="bismarkCoverage",
                       treatment=c(1,1,1,
                                   0,0,0),
                       context="CpG",
                       mincov = 1,
                       dbtype = "tabix",
                       dbdir = "shortterm_72hM_DB_merged"
)

covariates.72h <- data.frame(tank=c("AC11","AC11","AC5","AC5","AC6","AC6",
                                    "C11","C11","C5","C5","C6","C6"), 
                             sex=c("F","M","F","M","F","M",
                                   "F","M","F","M","F","M"), 
                             stringsAsFactors = TRUE)

run.DMS.DMR(myobj.72h,covariates.72h,keep.chr.noXY,6L,
                        "./72h/DMS_res_72h/myobj_72h_3X.RDS","./72h/DMS_res_72h/myobj_72h_5X.RDS",
                        "./72h/DMS_res_72h/DMSmeth_72h_5X.RDS","./72h/DMS_res_72h/meth_72h_5X_data.RDS",
                        "./72h/DMS_res_72h/DMSmydiff_72h_5X.RDS","./72h/DMS_res_72h/DMSmydiff_72h_5X_data.RDS",
                        "./72h/DMS_res_72h/DMSdiffmeth_72h_5X.RDS","./72h/DMS_res_72h/DMSdiffmeth_72h_5X_data.RDS",
                        "./72h/DMS_res_72h/DMSdiffmethchr_72h_5X.RDS",
                        "./72h/DMR_res_72h/DMRmeth_72h_5X.RDS","./72h/DMR_res_72h/DMRmeth_72h_5X_data.RDS",
                        "./72h/DMR_res_72h/DMR_tile_meth_72h_5X.RDS",
                        "./72h/DMR_res_72h/DMRmydiff_72h_5X.RDS","./72h/DMR_res_72h/DMRmydiff_72h_5X_data.RDS",
                        "./72h/DMR_res_72h/DMRdiffmeth_72h_5X.RDS","./72h/DMR_res_72h/DMRdiffmeth_72h_5X_data.RDS",
                        "./72h/DMR_res_72h/DMRdiffmethchr_72h_5X.RDS")

run.DMS.DMR(myobj.72h.fem,NULL,keep.chr.allchr,3L,
                        "./72h/DMS_res_72h/myobj_72h_fem_3X.RDS","./72h/DMS_res_72h/myobj_72h_fem_5X.RDS",
                        "./72h/DMS_res_72h/DMSmeth_72h_fem_5X.RDS","./72h/DMS_res_72h/meth_72h_fem_5X_data.RDS",
                        "./72h/DMS_res_72h/DMSmydiff_72h_fem_5X.RDS","./72h/DMS_res_72h/DMSmydiff_72h_fem_5X_data.RDS",
                        "./72h/DMS_res_72h/DMSdiffmeth_72h_fem_5X.RDS","./72h/DMS_res_72h/DMSdiffmeth_72h_fem_5X_data.RDS",
                        "./72h/DMS_res_72h/DMSdiffmethchr_72h_fem_5X.RDS",
                        "./72h/DMR_res_72h/DMRmeth_72h_fem_5X.RDS","./72h/DMR_res_72h/DMRmeth_72h_fem_5X_data.RDS",
                        "./72h/DMR_res_72h/DMR_tile_meth_72h_fem_5X.RDS",
                        "./72h/DMR_res_72h/DMRmydiff_72h_fem_5X.RDS","./72h/DMR_res_72h/DMRmydiff_72h_fem_5X_data.RDS",
                        "./72h/DMR_res_72h/DMRdiffmeth_72h_fem_5X.RDS","./72h/DMR_res_72h/DMRdiffmeth_72h_fem_5X_data.RDS",
                        "./72h/DMR_res_72h/DMRdiffmethchr_72h_fem_5X.RDS")

run.DMS.DMR(myobj.72h.mal,NULL,keep.chr.allchr,3L,
                        "./72h/DMS_res_72h/myobj_72h_mal_3X.RDS","./72h/DMS_res_72h/myobj_72h_mal_5X.RDS",
                        "./72h/DMS_res_72h/DMSmeth_72h_mal_5X.RDS","./72h/DMS_res_72h/meth_72h_mal_5X_data.RDS",
                        "./72h/DMS_res_72h/DMSmydiff_72h_mal_5X.RDS","./72h/DMS_res_72h/DMSmydiff_72h_mal_5X_data.RDS",
                        "./72h/DMS_res_72h/DMSdiffmeth_72h_mal_5X.RDS","./72h/DMS_res_72h/DMSdiffmeth_72h_mal_5X_data.RDS",
                        "./72h/DMS_res_72h/DMSdiffmethchr_72h_mal_5X.RDS",
                        "./72h/DMR_res_72h/DMRmeth_72h_mal_5X.RDS","./72h/DMR_res_72h/DMRmeth_72h_mal_5X_data.RDS",
                        "./72h/DMR_res_72h/DMR_tile_meth_72h_mal_5X.RDS",
                        "./72h/DMR_res_72h/DMRmydiff_72h_mal_5X.RDS","./72h/DMR_res_72h/DMRmydiff_72h_mal_5X_data.RDS",
                        "./72h/DMR_res_72h/DMRdiffmeth_72h_mal_5X.RDS","./72h/DMR_res_72h/DMRdiffmeth_72h_mal_5X_data.RDS",
                        "./72h/DMR_res_72h/DMRdiffmethchr_72h_mal_5X.RDS")

rm(myobj.72h,myobj.72h.fem,myobj.72h.mal)

# #create lists of file locations
file.list.all = list("../../../merged_sequences/st/ST2AC10F.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2AC10M.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2AC11F.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2AC11M.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2AC13F.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2AC13M.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2AC14F.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2AC14M.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2AC15F.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2AC15M.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2AC16F.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2AC16M.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2AC1F.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2AC1M.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2AC2F.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2AC2M.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2AC3F.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2AC3M.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2AC4F.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2AC4M.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2AC5F.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2AC5M.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2AC6F.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2AC6M.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2AC7F.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2AC7M.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2AC8F.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2AC8M.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2AC9F.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2AC9M.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2C10F.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2C10M.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2C11F.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2C11M.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2C12F.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2C12M.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2C13F.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2C13M.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2C14F.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2C14M.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2C15F.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2C15M.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2C1F.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2C1M.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2C2F.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2C2M.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2C3F.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2C3M.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2C4F.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2C4M.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2C5F.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2C5M.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2C6F.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2C6M.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2C7F.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2C7M.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2C8F.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2C8M.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2C9F.CpG_merged.cov",
                     "../../../merged_sequences/st/ST2C9M.CpG_merged.cov")

file.list.all.fem = list("../../../merged_sequences/st/ST2AC10F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC11F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC13F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC14F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC15F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC16F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC1F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC2F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC3F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC4F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC5F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC6F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC7F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC8F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC9F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C10F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C11F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C12F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C13F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C14F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C15F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C1F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C2F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C3F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C4F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C5F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C6F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C7F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C8F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C9F.CpG_merged.cov")

file.list.all.mal = list("../../../merged_sequences/st/ST2AC10M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC11M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC13M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC14M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC15M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC16M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC1M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC2M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC3M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC4M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC5M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC6M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC7M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC8M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC9M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C10M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C11M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C12M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C13M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C14M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C15M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C1M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C2M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C3M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C4M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C5M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C6M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C7M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C8M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C9M.CpG_merged.cov")

#create tabix file
myobj.all=methRead(file.list.all,
                   sample.id=list("ST2AC10F","ST2AC10M","ST2AC11F","ST2AC11M","ST2AC13F",
                                  "ST2AC13M","ST2AC14F","ST2AC14M","ST2AC15F","ST2AC15M",
                                  "ST2AC16F","ST2AC16M","ST2AC1F","ST2AC1M","ST2AC2F",
                                  "ST2AC2M","ST2AC3F","ST2AC3M","ST2AC4F","ST2AC4M",
                                  "ST2AC5F","ST2AC5M","ST2AC6F","ST2AC6M","ST2AC7F",
                                  "ST2AC7M","ST2AC8F","ST2AC8M","ST2AC9F","ST2AC9M",
                                  "ST2C10F","ST2C10M","ST2C11F","ST2C11M","ST2C12F",
                                  "ST2C12M","ST2C13F","ST2C13M","ST2C14F","ST2C14M",
                                  "ST2C15F","ST2C15M","ST2C1F","ST2C1M","ST2C2F",
                                  "ST2C2M","ST2C3F","ST2C3M","ST2C4F","ST2C4M",
                                  "ST2C5F","ST2C5M","ST2C6F","ST2C6M","ST2C7F",
                                  "ST2C7M","ST2C8F","ST2C8M","ST2C9F", "ST2C9M"),
                   assembly="guppyWGBS_shortterm",
                   pipeline="bismarkCoverage",
                   treatment=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                   context="CpG",
                   mincov = 1,
                   dbtype = "tabix",
                   dbdir = "shortterm_all_DB"
)

myobj.all.fem=methRead(file.list.all.fem,
                       sample.id=list("ST2AC10F","ST2AC11F","ST2AC13F","ST2AC14F","ST2AC15F",
                                      "ST2AC16F","ST2AC1F","ST2AC2F","ST2AC3F","ST2AC4F",
                                      "ST2AC5F","ST2AC6F","ST2AC7F","ST2AC8F","ST2AC9F",
                                      "ST2C10F","ST2C11F","ST2C12F","ST2C13F","ST2C14F",
                                      "ST2C15F","ST2C1F","ST2C2F","ST2C3F","ST2C4F",
                                      "ST2C5F","ST2C6F","ST2C7F","ST2C8F","ST2C9F"),
                       assembly="guppyWGBS_shortterm",
                       pipeline="bismarkCoverage",
                       treatment=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                       context="CpG",
                       mincov = 1,
                       dbtype = "tabix",
                       dbdir = "shortterm_allfem_DB"
)

myobj.all.mal=methRead(file.list.all.mal,
                       sample.id=list("ST2AC10M","ST2AC11M","ST2AC13M","ST2AC14M","ST2AC15M",
                                      "ST2AC16M","ST2AC1M","ST2AC2M","ST2AC3M","ST2AC4M",
                                      "ST2AC5M","ST2AC6M","ST2AC7M","ST2AC8M","ST2AC9M",
                                      "ST2C10M","ST2C11M","ST2C12M","ST2C13M","ST2C14M",
                                      "ST2C15M","ST2C1M","ST2C2M","ST2C3M","ST2C4M",
                                      "ST2C5M","ST2C6M","ST2C7M","ST2C8M", "ST2C9M"),
                       assembly="guppyWGBS_shortterm",
                       pipeline="bismarkCoverage",
                       treatment=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                       context="CpG",
                       mincov = 1,
                       dbtype = "tabix",
                       dbdir = "shortterm_allmal_DB"
)

# #enter covariates 
covariates.all <- data.frame(tank=c("AC10","AC10","AC11","AC11","AC13",
                                  "AC13","AC14","AC14","AC15","AC15",
                                  "AC16","AC16","AC1","AC1","AC2",
                                  "AC2","AC3","AC3","AC4","AC4",
                                  "AC5","AC5","AC6","AC6","AC7",
                                  "AC7","AC8","AC8","AC9","AC9",
                                  "C10","C10","C11","C11","C12",
                                  "C12","C13","C13","C14","C14",
                                  "C15","C15","C1","C1","C2",
                                  "C2","C3","C3","C4","C4",
                                  "C5","C5","C6","C6","C7",
                                  "C7","C8","C8","C9", "C9"), 
                            sex=c("F","M","F","M","F",
                                  "M","F","M","F","M",
                                  "F","M","F","M","F",
                                  "M","F","M","F","M",
                                  "F","M","F","M","F",
                                  "M","F","M","F","M",
                                  "F","M","F","M","F",
                                  "M","F","M","F","M",
                                  "F","M","F","M","F",
                                  "M","F","M","F","M",
                                  "F","M","F","M","F",
                                  "M","F","M","F", "M"), 
                            stringsAsFactors = TRUE)

run.DMS.DMR(myobj.all,covariates.all,keep.chr.noXY,30L,
                        "./all/DMS_res_all/myobj_all_3X.RDS","./all/DMS_res_all/myobj_all_5X.RDS",
                        "./all/DMS_res_all/DMSmeth_all_5X.RDS","./all/DMS_res_all/meth_all_5X_data.RDS",
                        "./all/DMS_res_all/DMSmydiff_all_5X.RDS","./all/DMS_res_all/DMSmydiff_all_5X_data.RDS",
                        "./all/DMS_res_all/DMSdiffmeth_all_5X.RDS","./all/DMS_res_all/DMSdiffmeth_all_5X_data.RDS",
                        "./all/DMS_res_all/DMSdiffmethchr_all_5X.RDS",
                        "./all/DMR_res_all/DMRmeth_all_5X.RDS","./all/DMR_res_all/DMRmeth_all_5X_data.RDS",
                        "./all/DMR_res_all/DMR_tile_meth_all_5X.RDS",
                        "./all/DMR_res_all/DMRmydiff_all_5X.RDS","./all/DMR_res_all/DMRmydiff_all_5X_data.RDS",
                        "./all/DMR_res_all/DMRdiffmeth_all_5X.RDS","./all/DMR_res_all/DMRdiffmeth_all_5X_data.RDS",
                        "./all/DMR_res_all/DMRdiffmethchr_all_5X.RDS")

run.DMS.DMR(myobj.all.fem,NULL,keep.chr.allchr,15L,
                        "./all/DMS_res_all/myobj_all_fem_3X.RDS","./all/DMS_res_all/myobj_all_fem_5X.RDS",
                        "./all/DMS_res_all/DMSmeth_all_fem_5X.RDS","./all/DMS_res_all/meth_all_fem_5X_data.RDS",
                        "./all/DMS_res_all/DMSmydiff_all_fem_5X.RDS","./all/DMS_res_all/DMSmydiff_all_fem_5X_data.RDS",
                        "./all/DMS_res_all/DMSdiffmeth_all_fem_5X.RDS","./all/DMS_res_all/DMSdiffmeth_all_fem_5X_data.RDS",
                        "./all/DMS_res_all/DMSdiffmethchr_all_fem_5X.RDS",
                        "./all/DMR_res_all/DMRmeth_all_fem_5X.RDS","./all/DMR_res_all/DMRmeth_all_fem_5X_data.RDS",
                        "./all/DMR_res_all/DMR_tile_meth_all_fem_5X.RDS",
                        "./all/DMR_res_all/DMRmydiff_all_fem_5X.RDS","./all/DMR_res_all/DMRmydiff_all_fem_5X_data.RDS",
                        "./all/DMR_res_all/DMRdiffmeth_all_fem_5X.RDS","./all/DMR_res_all/DMRdiffmeth_all_fem_5X_data.RDS",
                        "./all/DMR_res_all/DMRdiffmethchr_all_fem_5X.RDS")

run.DMS.DMR(myobj.all.mal,NULL,keep.chr.allchr,15L,
                        "./all/DMS_res_all/myobj_all_mal_3X.RDS","./all/DMS_res_all/myobj_all_mal_5X.RDS",
                        "./all/DMS_res_all/DMSmeth_all_mal_5X.RDS","./all/DMS_res_all/meth_all_mal_5X_data.RDS",
                        "./all/DMS_res_all/DMSmydiff_all_mal_5X.RDS","./all/DMS_res_all/DMSmydiff_all_mal_5X_data.RDS",
                        "./all/DMS_res_all/DMSdiffmeth_all_mal_5X.RDS","./all/DMS_res_all/DMSdiffmeth_all_mal_5X_data.RDS",
                        "./all/DMS_res_all/DMSdiffmethchr_all_mal_5X.RDS",
                        "./all/DMR_res_all/DMRmeth_all_mal_5X.RDS","./all/DMR_res_all/DMRmeth_all_mal_5X_data.RDS",
                        "./all/DMR_res_all/DMR_tile_meth_all_mal_5X.RDS",
                        "./all/DMR_res_all/DMRmydiff_all_mal_5X.RDS","./all/DMR_res_all/DMRmydiff_all_mal_5X_data.RDS",
                        "./all/DMR_res_all/DMRdiffmeth_all_mal_5X.RDS","./all/DMR_res_all/DMRdiffmeth_all_mal_5X_data.RDS",
                        "./all/DMR_res_all/DMRdiffmethchr_all_mal_5X.RDS")

rm(myobj.all,myobj.all.fem,myobj.all.mal)