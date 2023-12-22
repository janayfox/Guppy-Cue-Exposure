#######################################################
### Goal: Run and plot GO-seq analysis
### Author: Janay Fox
### R script
#######################################################

## Set up ##
#install packages
#install.packages("ggplot2")
#install.packages("dplyr")

#load packaes
library(ggplot2)
library(dplyr)
library(genomation)
library(DescTools)
library(tibble)
library(GenomicRanges)
library(methylKit)
library(GOstats)
library(biomaRt)
library(GSEABase)
library(ggpubr)

#load in data
DMS_TSS_all <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_all/DMS_TSS_all.RDS")
DMS_TSS_05h <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_05h/DMS_TSS_05h.RDS")
DMS_TSS_1h <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_1h/DMS_TSS_1h.RDS")
DMS_TSS_4h <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_4h/DMS_TSS_4h.RDS")
DMS_TSS_24h <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_24h/DMS_TSS_24h.RDS")
DMS_TSS_72h <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_72h/DMS_TSS_72h.RDS")

DMS_TSS_all_fem <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_all/DMS_TSS_all_fem.RDS")
DMS_TSS_05h_fem <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_05h/DMS_TSS_05h_fem.RDS")
DMS_TSS_1h_fem <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_1h/DMS_TSS_1h_fem.RDS")
DMS_TSS_4h_fem <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_4h/DMS_TSS_4h_fem.RDS")
DMS_TSS_24h_fem <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_24h/DMS_TSS_24h_fem.RDS")
DMS_TSS_72h_fem <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_72h/DMS_TSS_72h_fem.RDS")

DMS_TSS_all_mal <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_all/DMS_TSS_all_mal.RDS")
DMS_TSS_05h_mal <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_05h/DMS_TSS_05h_mal.RDS")
DMS_TSS_1h_mal <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_1h/DMS_TSS_1h_mal.RDS")
DMS_TSS_4h_mal <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_4h/DMS_TSS_4h_mal.RDS")
DMS_TSS_24h_mal <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_24h/DMS_TSS_24h_mal.RDS")
DMS_TSS_72h_mal <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_72h/DMS_TSS_72h_mal.RDS")

DMR_TSS_all <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_all/DMR_TSS_all.RDS")
DMR_TSS_05h <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_05h/DMR_TSS_05h.RDS")
DMR_TSS_1h <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_1h/DMR_TSS_1h.RDS")
DMR_TSS_4h <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_4h/DMR_TSS_4h.RDS")
DMR_TSS_24h <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_24h/DMR_TSS_24h.RDS")
DMR_TSS_72h <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_72h/DMR_TSS_72h.RDS")

DMR_TSS_all_fem <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_all/DMR_TSS_all_fem.RDS")
DMR_TSS_05h_fem <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_05h/DMR_TSS_05h_fem.RDS")
DMR_TSS_1h_fem <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_1h/DMR_TSS_1h_fem.RDS")
DMR_TSS_4h_fem <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_4h/DMR_TSS_4h_fem.RDS")
DMR_TSS_24h_fem <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_24h/DMR_TSS_24h_fem.RDS")
DMR_TSS_72h_fem <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_72h/DMR_TSS_72h_fem.RDS")

DMR_TSS_all_mal <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_all/DMR_TSS_all_mal.RDS")
DMR_TSS_05h_mal <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_05h/DMR_TSS_05h_mal.RDS")
DMR_TSS_1h_mal <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_1h/DMR_TSS_1h_mal.RDS")
DMR_TSS_4h_mal <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_4h/DMR_TSS_4h_mal.RDS")
DMR_TSS_24h_mal <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_24h/DMR_TSS_24h_mal.RDS")
DMR_TSS_72h_mal <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_72h/DMR_TSS_72h_mal.RDS")

CpG_anno_all <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_all/CpG_anno_all.RDS")
CpG_anno_05h <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_05h/CpG_anno_05h.RDS")
CpG_anno_1h <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_1h/CpG_anno_1h.RDS")
CpG_anno_4h <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_4h/CpG_anno_4h.RDS")
CpG_anno_24h <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_24h/CpG_anno_24h.RDS")
CpG_anno_72h <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_72h/CpG_anno_72h.RDS")

CpG_anno_all_fem <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_all/CpG_anno_all_fem.RDS")
CpG_anno_05h_fem <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_05h/CpG_anno_05h_fem.RDS")
CpG_anno_1h_fem <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_1h/CpG_anno_1h_fem.RDS")
CpG_anno_4h_fem <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_4h/CpG_anno_4h_fem.RDS")
CpG_anno_24h_fem <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_24h/CpG_anno_24h_fem.RDS")
CpG_anno_72h_fem <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_72h/CpG_anno_72h_fem.RDS")

CpG_anno_all_mal <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_all/CpG_anno_all_mal.RDS")
CpG_anno_05h_mal <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_05h/CpG_anno_05h_mal.RDS")
CpG_anno_1h_mal <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_1h/CpG_anno_1h_mal.RDS")
CpG_anno_4h_mal <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_4h/CpG_anno_4h_mal.RDS")
CpG_anno_24h_mal <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_24h/CpG_anno_24h_mal.RDS")
CpG_anno_72h_mal <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_72h/CpG_anno_72h_mal.RDS")

## Go seq analysis ##
#get proper biomart
bm <- useEnsembl(biomart = "ensembl")
bm <- useDataset("preticulata_gene_ensembl", mart = bm)

#grab go terms
go_list <- getBM(mart = bm, attributes = c('ensembl_transcript_id','external_gene_name','go_id', "name_1006", 
                                                                              "namespace_1003", "go_linkage_type"))
#prep go terms for gene mappings
go_list_data <- go_list[,c(3,6,1)] #limit to needed data
go_list_data <- go_list_data[go_list_data$go_id != '',] #remove empty go terms

#set go frame
goFrame <- GOFrame(go_list_data, organism = "guppy")
goAllFrame <- GOAllFrame(goFrame)
#set gene set
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())

run.gostat <- function(test.data, cpg.data, ont){
  #get background gene sets from all CpGs
  cpg.data <- getAssociationWithTSS(cpg.data)
  
  #limit to genes within 10kb 
  test.data <- subset(test.data, dist.to.feature >= -10000 & dist.to.feature <= 10000)
  
  universe <- cpg.data$feature.name
  genes <- test.data$feature.name
  params <- GSEAGOHyperGParams(name = "shorterm GSEA based annot Params",
                               geneSetCollection = gsc,
                               geneIds = genes,
                               universeGeneIds = universe,
                               ontology = ont,
                               pvalueCutoff = 0.05,
                               conditional = FALSE,
                               testDirection = "over")
  
  over <- hyperGTest(params)
  print(over)
  sum.over <- summary(over)
  sum.over$FDR <-  p.adjust(sum.over$Pvalue, method = "fdr")
  return(sum.over)
}

#run on all gene sets
DMS_all_gostat_BP <- run.gostat(DMS_TSS_all, CpG_anno_all, "BP")
DMS_05h_gostat_BP <- run.gostat(DMS_TSS_05h, CpG_anno_05h, "BP")
DMS_1h_gostat_BP <- run.gostat(DMS_TSS_1h, CpG_anno_1h, "BP")
DMS_4h_gostat_BP <- run.gostat(DMS_TSS_4h, CpG_anno_4h, "BP")
DMS_24h_gostat_BP <- run.gostat(DMS_TSS_24h, CpG_anno_24h, "BP")
DMS_72h_gostat_BP<- run.gostat(DMS_TSS_72h, CpG_anno_72h, "BP")

DMS_all_fem_gostat_BP <- run.gostat(DMS_TSS_all_fem, CpG_anno_all_fem, "BP")
DMS_05h_fem_gostat_BP <- run.gostat(DMS_TSS_05h_fem, CpG_anno_05h_fem, "BP")
DMS_1h_fem_gostat_BP <- run.gostat(DMS_TSS_1h_fem, CpG_anno_1h_fem, "BP")
DMS_4h_fem_gostat_BP<- run.gostat(DMS_TSS_4h_fem, CpG_anno_4h_fem, "BP")
DMS_24h_fem_gostat_BP <- run.gostat(DMS_TSS_24h_fem, CpG_anno_24h_fem, "BP")
DMS_72h_fem_gostat_BP <- run.gostat(DMS_TSS_72h_fem, CpG_anno_72h_fem, "BP")

DMS_all_mal_gostat_BP <- run.gostat(DMS_TSS_all_mal, CpG_anno_all_mal, "BP")
DMS_05h_mal_gostat_BP <- run.gostat(DMS_TSS_05h_mal, CpG_anno_05h_mal, "BP")
DMS_1h_mal_gostat_BP <- run.gostat(DMS_TSS_1h_mal, CpG_anno_1h_mal, "BP")
DMS_4h_mal_gostat_BP <- run.gostat(DMS_TSS_4h_mal, CpG_anno_4h_mal, "BP")
DMS_24h_mal_gostat_BP <- run.gostat(DMS_TSS_24h_mal, CpG_anno_24h_mal, "BP")
DMS_72h_mal_gostat_BP <- run.gostat(DMS_TSS_72h_mal, CpG_anno_72h_mal, "BP")

DMS_all_gostat_MF <- run.gostat(DMS_TSS_all, CpG_anno_all, "MF")
DMS_05h_gostat_MF <- run.gostat(DMS_TSS_05h, CpG_anno_05h, "MF")
DMS_1h_gostat_MF <- run.gostat(DMS_TSS_1h, CpG_anno_1h, "MF")
DMS_4h_gostat_MF <- run.gostat(DMS_TSS_4h, CpG_anno_4h, "MF")
DMS_24h_gostat_MF <- run.gostat(DMS_TSS_24h, CpG_anno_24h, "MF")
DMS_72h_gostat_MF<- run.gostat(DMS_TSS_72h, CpG_anno_72h, "MF")

DMS_all_fem_gostat_MF <- run.gostat(DMS_TSS_all_fem, CpG_anno_all_fem, "MF")
DMS_05h_fem_gostat_MF <- run.gostat(DMS_TSS_05h_fem, CpG_anno_05h_fem, "MF")
DMS_1h_fem_gostat_MF <- run.gostat(DMS_TSS_1h_fem, CpG_anno_1h_fem, "MF")
DMS_4h_fem_gostat_MF<- run.gostat(DMS_TSS_4h_fem, CpG_anno_4h_fem, "MF")
DMS_24h_fem_gostat_MF <- run.gostat(DMS_TSS_24h_fem, CpG_anno_24h_fem, "MF")
DMS_72h_fem_gostat_MF <- run.gostat(DMS_TSS_72h_fem, CpG_anno_72h_fem, "MF")

DMS_all_mal_gostat_MF <- run.gostat(DMS_TSS_all_mal, CpG_anno_all_mal, "MF")
DMS_05h_mal_gostat_MF <- run.gostat(DMS_TSS_05h_mal, CpG_anno_05h_mal, "MF")
DMS_1h_mal_gostat_MF <- run.gostat(DMS_TSS_1h_mal, CpG_anno_1h_mal, "MF")
DMS_4h_mal_gostat_MF <- run.gostat(DMS_TSS_4h_mal, CpG_anno_4h_mal, "MF")
DMS_24h_mal_gostat_MF <- run.gostat(DMS_TSS_24h_mal, CpG_anno_24h_mal, "MF")
DMS_72h_mal_gostat_MF <- run.gostat(DMS_TSS_72h_mal, CpG_anno_72h_mal, "MF")

DMS_all_gostat_CC <- run.gostat(DMS_TSS_all, CpG_anno_all, "CC")
DMS_05h_gostat_CC <- run.gostat(DMS_TSS_05h, CpG_anno_05h, "CC")
DMS_1h_gostat_CC <- run.gostat(DMS_TSS_1h, CpG_anno_1h, "CC")
DMS_4h_gostat_CC <- run.gostat(DMS_TSS_4h, CpG_anno_4h, "CC")
DMS_24h_gostat_CC <- run.gostat(DMS_TSS_24h, CpG_anno_24h, "CC")
DMS_72h_gostat_CC<- run.gostat(DMS_TSS_72h, CpG_anno_72h, "CC")

DMS_all_fem_gostat_CC <- run.gostat(DMS_TSS_all_fem, CpG_anno_all_fem, "CC")
DMS_05h_fem_gostat_CC <- run.gostat(DMS_TSS_05h_fem, CpG_anno_05h_fem, "CC")
DMS_1h_fem_gostat_CC <- run.gostat(DMS_TSS_1h_fem, CpG_anno_1h_fem, "CC")
DMS_4h_fem_gostat_CC<- run.gostat(DMS_TSS_4h_fem, CpG_anno_4h_fem, "CC")
DMS_24h_fem_gostat_CC <- run.gostat(DMS_TSS_24h_fem, CpG_anno_24h_fem, "CC")
DMS_72h_fem_gostat_CC <- run.gostat(DMS_TSS_72h_fem, CpG_anno_72h_fem, "CC")

DMS_all_mal_gostat_CC <- run.gostat(DMS_TSS_all_mal, CpG_anno_all_mal, "CC")
DMS_05h_mal_gostat_CC <- run.gostat(DMS_TSS_05h_mal, CpG_anno_05h_mal, "CC")
DMS_1h_mal_gostat_CC <- run.gostat(DMS_TSS_1h_mal, CpG_anno_1h_mal, "CC")
DMS_4h_mal_gostat_CC <- run.gostat(DMS_TSS_4h_mal, CpG_anno_4h_mal, "CC")
DMS_24h_mal_gostat_CC <- run.gostat(DMS_TSS_24h_mal, CpG_anno_24h_mal, "CC")
DMS_72h_mal_gostat_CC <- run.gostat(DMS_TSS_72h_mal, CpG_anno_72h_mal, "CC")

DMR_all_gostat_BP <- run.gostat(DMR_TSS_all, CpG_anno_all, "BP")
DMR_05h_gostat_BP <- run.gostat(DMR_TSS_05h, CpG_anno_05h, "BP")
DMR_1h_gostat_BP <- run.gostat(DMR_TSS_1h, CpG_anno_1h, "BP")
DMR_4h_gostat_BP <- run.gostat(DMR_TSS_4h, CpG_anno_4h, "BP")
DMR_24h_gostat_BP <- run.gostat(DMR_TSS_24h, CpG_anno_24h, "BP")
DMR_72h_gostat_BP<- run.gostat(DMR_TSS_72h, CpG_anno_72h, "BP")

DMR_all_fem_gostat_BP <- run.gostat(DMR_TSS_all_fem, CpG_anno_all_fem, "BP")
DMR_05h_fem_gostat_BP <- run.gostat(DMR_TSS_05h_fem, CpG_anno_05h_fem, "BP")
DMR_1h_fem_gostat_BP <- run.gostat(DMR_TSS_1h_fem, CpG_anno_1h_fem, "BP")
DMR_4h_fem_gostat_BP<- run.gostat(DMR_TSS_4h_fem, CpG_anno_4h_fem, "BP")
DMR_24h_fem_gostat_BP <- run.gostat(DMR_TSS_24h_fem, CpG_anno_24h_fem, "BP")
DMR_72h_fem_gostat_BP <- run.gostat(DMR_TSS_72h_fem, CpG_anno_72h_fem, "BP")

DMR_all_mal_gostat_BP <- run.gostat(DMR_TSS_all_mal, CpG_anno_all_mal, "BP")
DMR_05h_mal_gostat_BP <- run.gostat(DMR_TSS_05h_mal, CpG_anno_05h_mal, "BP")
DMR_1h_mal_gostat_BP <- run.gostat(DMR_TSS_1h_mal, CpG_anno_1h_mal, "BP")
DMR_4h_mal_gostat_BP <- run.gostat(DMR_TSS_4h_mal, CpG_anno_4h_mal, "BP")
DMR_24h_mal_gostat_BP <- run.gostat(DMR_TSS_24h_mal, CpG_anno_24h_mal, "BP")
DMR_72h_mal_gostat_BP <- run.gostat(DMR_TSS_72h_mal, CpG_anno_72h_mal, "BP")

DMR_all_gostat_MF <- run.gostat(DMR_TSS_all, CpG_anno_all, "MF")
DMR_05h_gostat_MF <- run.gostat(DMR_TSS_05h, CpG_anno_05h, "MF")
DMR_1h_gostat_MF <- run.gostat(DMR_TSS_1h, CpG_anno_1h, "MF")
DMR_4h_gostat_MF <- run.gostat(DMR_TSS_4h, CpG_anno_4h, "MF")
DMR_24h_gostat_MF <- run.gostat(DMR_TSS_24h, CpG_anno_24h, "MF")
DMR_72h_gostat_MF<- run.gostat(DMR_TSS_72h, CpG_anno_72h, "MF")

DMR_all_fem_gostat_MF <- run.gostat(DMR_TSS_all_fem, CpG_anno_all_fem, "MF")
DMR_05h_fem_gostat_MF <- run.gostat(DMR_TSS_05h_fem, CpG_anno_05h_fem, "MF")
DMR_1h_fem_gostat_MF <- run.gostat(DMR_TSS_1h_fem, CpG_anno_1h_fem, "MF")
DMR_4h_fem_gostat_MF<- run.gostat(DMR_TSS_4h_fem, CpG_anno_4h_fem, "MF")
DMR_24h_fem_gostat_MF <- run.gostat(DMR_TSS_24h_fem, CpG_anno_24h_fem, "MF")
DMR_72h_fem_gostat_MF <- run.gostat(DMR_TSS_72h_fem, CpG_anno_72h_fem, "MF")

DMR_all_mal_gostat_MF <- run.gostat(DMR_TSS_all_mal, CpG_anno_all_mal, "MF")
DMR_05h_mal_gostat_MF <- run.gostat(DMR_TSS_05h_mal, CpG_anno_05h_mal, "MF")
DMR_1h_mal_gostat_MF <- run.gostat(DMR_TSS_1h_mal, CpG_anno_1h_mal, "MF")
DMR_4h_mal_gostat_MF <- run.gostat(DMR_TSS_4h_mal, CpG_anno_4h_mal, "MF")
DMR_24h_mal_gostat_MF <- run.gostat(DMR_TSS_24h_mal, CpG_anno_24h_mal, "MF")
DMR_72h_mal_gostat_MF <- run.gostat(DMR_TSS_72h_mal, CpG_anno_72h_mal, "MF")

DMR_all_gostat_CC <- run.gostat(DMR_TSS_all, CpG_anno_all, "CC")
DMR_05h_gostat_CC <- run.gostat(DMR_TSS_05h, CpG_anno_05h, "CC")
DMR_1h_gostat_CC <- run.gostat(DMR_TSS_1h, CpG_anno_1h, "CC")
DMR_4h_gostat_CC <- run.gostat(DMR_TSS_4h, CpG_anno_4h, "CC")
DMR_24h_gostat_CC <- run.gostat(DMR_TSS_24h, CpG_anno_24h, "CC")
DMR_72h_gostat_CC<- run.gostat(DMR_TSS_72h, CpG_anno_72h, "CC")

DMR_all_fem_gostat_CC <- run.gostat(DMR_TSS_all_fem, CpG_anno_all_fem, "CC")
DMR_05h_fem_gostat_CC <- run.gostat(DMR_TSS_05h_fem, CpG_anno_05h_fem, "CC")
DMR_1h_fem_gostat_CC <- run.gostat(DMR_TSS_1h_fem, CpG_anno_1h_fem, "CC")
DMR_4h_fem_gostat_CC<- run.gostat(DMR_TSS_4h_fem, CpG_anno_4h_fem, "CC")
DMR_24h_fem_gostat_CC <- run.gostat(DMR_TSS_24h_fem, CpG_anno_24h_fem, "CC")
DMR_72h_fem_gostat_CC <- run.gostat(DMR_TSS_72h_fem, CpG_anno_72h_fem, "CC")

DMR_all_mal_gostat_CC <- run.gostat(DMR_TSS_all_mal, CpG_anno_all_mal, "CC")
DMR_05h_mal_gostat_CC <- run.gostat(DMR_TSS_05h_mal, CpG_anno_05h_mal, "CC")
DMR_1h_mal_gostat_CC <- run.gostat(DMR_TSS_1h_mal, CpG_anno_1h_mal, "CC")
DMR_4h_mal_gostat_CC <- run.gostat(DMR_TSS_4h_mal, CpG_anno_4h_mal, "CC")
DMR_24h_mal_gostat_CC <- run.gostat(DMR_TSS_24h_mal, CpG_anno_24h_mal, "CC")
DMR_72h_mal_gostat_CC <- run.gostat(DMR_TSS_72h_mal, CpG_anno_72h_mal, "CC")

#save results as a table 
write.csv(DMS_all_gostat_BP, "./shortTerm_exp/data/go_res/DMS_all_gostat_BP.csv", row.names = FALSE)
write.csv(DMS_05h_gostat_BP, "./shortTerm_exp/data/go_res/DMS_05h_gostat_BP.csv", row.names = FALSE)
write.csv(DMS_1h_gostat_BP, "./shortTerm_exp/data/go_res/DMS_1h_gostat_BP.csv", row.names = FALSE)
write.csv(DMS_4h_gostat_BP, "./shortTerm_exp/data/go_res/DMS_4h_gostat_BP.csv", row.names = FALSE)
write.csv(DMS_24h_gostat_BP, "./shortTerm_exp/data/go_res/DMS_24h_gostat_BP.csv", row.names = FALSE)
write.csv(DMS_72h_gostat_BP, "./shortTerm_exp/data/go_res/DMS_72h_gostat_BP.csv", row.names = FALSE)

write.csv(DMS_all_fem_gostat_BP, "./shortTerm_exp/data/go_res/DMS_all_fem_gostat_BP.csv", row.names = FALSE)
write.csv(DMS_05h_fem_gostat_BP, "./shortTerm_exp/data/go_res/DMS_05h_fem_gostat_BP.csv", row.names = FALSE)
write.csv(DMS_1h_fem_gostat_BP, "./shortTerm_exp/data/go_res/DMS_1h_fem_gostat_BP.csv", row.names = FALSE)
write.csv(DMS_4h_fem_gostat_BP, "./shortTerm_exp/data/go_res/DMS_4h_fem_gostat_BP.csv", row.names = FALSE)
write.csv(DMS_24h_fem_gostat_BP, "./shortTerm_exp/data/go_res/DMS_24h_fem_gostat_BP.csv", row.names = FALSE)
write.csv(DMS_72h_fem_gostat_BP, "./shortTerm_exp/data/go_res/DMS_72h_fem_gostat_BP.csv", row.names = FALSE)

write.csv(DMS_all_mal_gostat_BP, "./shortTerm_exp/data/go_res/DMS_all_mal_gostat_BP.csv", row.names = FALSE)
write.csv(DMS_05h_mal_gostat_BP, "./shortTerm_exp/data/go_res/DMS_05h_mal_gostat_BP.csv", row.names = FALSE)
write.csv(DMS_1h_mal_gostat_BP, "./shortTerm_exp/data/go_res/DMS_1h_mal_gostat_BP.csv", row.names = FALSE)
write.csv(DMS_4h_mal_gostat_BP, "./shortTerm_exp/data/go_res/DMS_4h_mal_gostat_BP.csv", row.names = FALSE)
write.csv(DMS_24h_mal_gostat_BP, "./shortTerm_exp/data/go_res/DMS_24h_mal_gostat_BP.csv", row.names = FALSE)
write.csv(DMS_72h_mal_gostat_BP, "./shortTerm_exp/data/go_res/DMS_72h_mal_gostat_BP.csv", row.names = FALSE)

write.csv(DMS_all_gostat_MF, "./shortTerm_exp/data/go_res/DMS_all_gostat_MF.csv", row.names = FALSE)
write.csv(DMS_05h_gostat_MF, "./shortTerm_exp/data/go_res/DMS_05h_gostat_MF.csv", row.names = FALSE)
write.csv(DMS_1h_gostat_MF, "./shortTerm_exp/data/go_res/DMS_1h_gostat_MF.csv", row.names = FALSE)
write.csv(DMS_4h_gostat_MF, "./shortTerm_exp/data/go_res/DMS_4h_gostat_MF.csv", row.names = FALSE)
write.csv(DMS_24h_gostat_MF, "./shortTerm_exp/data/go_res/DMS_24h_gostat_MF.csv", row.names = FALSE)
write.csv(DMS_72h_gostat_MF, "./shortTerm_exp/data/go_res/DMS_72h_gostat_MF.csv", row.names = FALSE)

write.csv(DMS_all_fem_gostat_MF, "./shortTerm_exp/data/go_res/DMS_all_fem_gostat_MF.csv", row.names = FALSE)
write.csv(DMS_05h_fem_gostat_MF, "./shortTerm_exp/data/go_res/DMS_05h_fem_gostat_MF.csv", row.names = FALSE)
write.csv(DMS_1h_fem_gostat_MF, "./shortTerm_exp/data/go_res/DMS_1h_fem_gostat_MF.csv", row.names = FALSE)
write.csv(DMS_4h_fem_gostat_MF, "./shortTerm_exp/data/go_res/DMS_4h_fem_gostat_MF.csv", row.names = FALSE)
write.csv(DMS_24h_fem_gostat_MF, "./shortTerm_exp/data/go_res/DMS_24h_fem_gostat_MF.csv", row.names = FALSE)
write.csv(DMS_72h_fem_gostat_MF, "./shortTerm_exp/data/go_res/DMS_72h_fem_gostat_MF.csv", row.names = FALSE)

write.csv(DMS_all_mal_gostat_MF, "./shortTerm_exp/data/go_res/DMS_all_mal_gostat_MF.csv", row.names = FALSE)
write.csv(DMS_05h_mal_gostat_MF, "./shortTerm_exp/data/go_res/DMS_05h_mal_gostat_MF.csv", row.names = FALSE)
write.csv(DMS_1h_mal_gostat_MF, "./shortTerm_exp/data/go_res/DMS_1h_mal_gostat_MF.csv", row.names = FALSE)
write.csv(DMS_4h_mal_gostat_MF, "./shortTerm_exp/data/go_res/DMS_4h_mal_gostat_MF.csv", row.names = FALSE)
write.csv(DMS_24h_mal_gostat_MF, "./shortTerm_exp/data/go_res/DMS_24h_mal_gostat_MF.csv", row.names = FALSE)
write.csv(DMS_72h_mal_gostat_MF, "./shortTerm_exp/data/go_res/DMS_72h_mal_gostat_MF.csv", row.names = FALSE)

write.csv(DMS_all_gostat_CC, "./shortTerm_exp/data/go_res/DMS_all_gostat_CC.csv", row.names = FALSE)
write.csv(DMS_05h_gostat_CC, "./shortTerm_exp/data/go_res/DMS_05h_gostat_CC.csv", row.names = FALSE)
write.csv(DMS_1h_gostat_CC, "./shortTerm_exp/data/go_res/DMS_1h_gostat_CC.csv", row.names = FALSE)
write.csv(DMS_4h_gostat_CC, "./shortTerm_exp/data/go_res/DMS_4h_gostat_CC.csv", row.names = FALSE)
write.csv(DMS_24h_gostat_CC, "./shortTerm_exp/data/go_res/DMS_24h_gostat_CC.csv", row.names = FALSE)
write.csv(DMS_72h_gostat_CC, "./shortTerm_exp/data/go_res/DMS_72h_gostat_CC.csv", row.names = FALSE)

write.csv(DMS_all_fem_gostat_CC, "./shortTerm_exp/data/go_res/DMS_all_fem_gostat_CC.csv", row.names = FALSE)
write.csv(DMS_05h_fem_gostat_CC, "./shortTerm_exp/data/go_res/DMS_05h_fem_gostat_CC.csv", row.names = FALSE)
write.csv(DMS_1h_fem_gostat_CC, "./shortTerm_exp/data/go_res/DMS_1h_fem_gostat_CC.csv", row.names = FALSE)
write.csv(DMS_4h_fem_gostat_CC, "./shortTerm_exp/data/go_res/DMS_4h_fem_gostat_CC.csv", row.names = FALSE)
write.csv(DMS_24h_fem_gostat_CC, "./shortTerm_exp/data/go_res/DMS_24h_fem_gostat_CC.csv", row.names = FALSE)
write.csv(DMS_72h_fem_gostat_CC, "./shortTerm_exp/data/go_res/DMS_72h_fem_gostat_CC.csv", row.names = FALSE)

write.csv(DMS_all_mal_gostat_CC, "./shortTerm_exp/data/go_res/DMS_all_mal_gostat_CC.csv", row.names = FALSE)
write.csv(DMS_05h_mal_gostat_CC, "./shortTerm_exp/data/go_res/DMS_05h_mal_gostat_CC.csv", row.names = FALSE)
write.csv(DMS_1h_mal_gostat_CC, "./shortTerm_exp/data/go_res/DMS_1h_mal_gostat_CC.csv", row.names = FALSE)
write.csv(DMS_4h_mal_gostat_CC, "./shortTerm_exp/data/go_res/DMS_4h_mal_gostat_CC.csv", row.names = FALSE)
write.csv(DMS_24h_mal_gostat_CC, "./shortTerm_exp/data/go_res/DMS_24h_mal_gostat_CC.csv", row.names = FALSE)
write.csv(DMS_72h_mal_gostat_CC, "./shortTerm_exp/data/go_res/DMS_72h_mal_gostat_CC.csv", row.names = FALSE)

#merge dataframse and save
DMS_05h_fem_gostat_BP$time_point <- "0.5h"
DMS_1h_fem_gostat_BP$time_point <- "1h"
DMS_4h_fem_gostat_BP$time_point <- "4h"
DMS_24h_fem_gostat_BP$time_point <- "24h"
DMS_72h_fem_gostat_BP$time_point <- "72h"

DMS_05h_fem_gostat_MF$time_point <- "0.5h"
DMS_1h_fem_gostat_MF$time_point <- "1h"
DMS_4h_fem_gostat_MF$time_point <- "4h"
DMS_24h_fem_gostat_MF$time_point <- "24h"
DMS_72h_fem_gostat_MF$time_point <- "72h"

DMS_05h_fem_gostat_CC$time_point <- "0.5h"
DMS_1h_fem_gostat_CC$time_point <- "1h"
DMS_4h_fem_gostat_CC$time_point <- "4h"
DMS_24h_fem_gostat_CC$time_point <- "24h"
DMS_72h_fem_gostat_CC$time_point <- "72h"

DMS_fem_all_BP <- rbind(DMS_05h_fem_gostat_BP,DMS_1h_fem_gostat_BP,DMS_4h_fem_gostat_BP,DMS_24h_fem_gostat_BP,DMS_72h_fem_gostat_BP)
DMS_fem_all_MF <- rbind(DMS_05h_fem_gostat_MF,DMS_1h_fem_gostat_MF,DMS_4h_fem_gostat_MF,DMS_24h_fem_gostat_MF,DMS_72h_fem_gostat_MF)
DMS_fem_all_CC <- rbind(DMS_05h_fem_gostat_CC,DMS_1h_fem_gostat_CC,DMS_4h_fem_gostat_CC,DMS_24h_fem_gostat_CC,DMS_72h_fem_gostat_CC)

write.csv(DMS_fem_all_BP, "./shortTerm_exp/data/go_res/DMS_fem_combined_gostat_BP.csv", row.names = FALSE)
write.csv(DMS_fem_all_MF, "./shortTerm_exp/data/go_res/DMS_fem_combined_gostat_MF.csv", row.names = FALSE)
write.csv(DMS_fem_all_CC, "./shortTerm_exp/data/go_res/DMS_fem_combined_gostat_CC.csv", row.names = FALSE)

DMS_05h_mal_gostat_BP$time_point <- "0.5h"
DMS_1h_mal_gostat_BP$time_point <- "1h"
DMS_4h_mal_gostat_BP$time_point <- "4h"
DMS_24h_mal_gostat_BP$time_point <- "24h"
DMS_72h_mal_gostat_BP$time_point <- "72h"

DMS_05h_mal_gostat_MF$time_point <- "0.5h"
DMS_1h_mal_gostat_MF$time_point <- "1h"
DMS_4h_mal_gostat_MF$time_point <- "4h"
DMS_24h_mal_gostat_MF$time_point <- "24h"
DMS_72h_mal_gostat_MF$time_point <- "72h"

DMS_05h_mal_gostat_CC$time_point <- "0.5h"
DMS_1h_mal_gostat_CC$time_point <- "1h"
DMS_4h_mal_gostat_CC$time_point <- "4h"
DMS_24h_mal_gostat_CC$time_point <- "24h"
DMS_72h_mal_gostat_CC$time_point <- "72h"


DMS_mal_all_BP <- rbind(DMS_05h_mal_gostat_BP,DMS_1h_mal_gostat_BP,DMS_4h_mal_gostat_BP,DMS_24h_mal_gostat_BP,DMS_72h_mal_gostat_BP)
DMS_mal_all_MF <- rbind(DMS_05h_mal_gostat_MF,DMS_1h_mal_gostat_MF,DMS_4h_mal_gostat_MF,DMS_24h_mal_gostat_MF,DMS_72h_mal_gostat_MF)
DMS_mal_all_CC <- rbind(DMS_05h_mal_gostat_CC,DMS_1h_mal_gostat_CC,DMS_4h_mal_gostat_CC,DMS_24h_mal_gostat_CC,DMS_72h_mal_gostat_CC)

write.csv(DMS_mal_all_BP, "./shortTerm_exp/data/go_res/DMS_mal_combined_gostat_BP.csv", row.names = FALSE)
write.csv(DMS_mal_all_MF, "./shortTerm_exp/data/go_res/DMS_mal_combined_gostat_MF.csv", row.names = FALSE)
write.csv(DMS_mal_all_CC, "./shortTerm_exp/data/go_res/DMS_mal_combined_gostat_CC.csv", row.names = FALSE)

DMR_05h_fem_gostat_BP$time_point <- "0.5h"
DMR_1h_fem_gostat_BP$time_point <- "1h"
DMR_4h_fem_gostat_BP$time_point <- "4h"
DMR_24h_fem_gostat_BP$time_point <- "24h"
DMR_72h_fem_gostat_BP$time_point <- "72h"

DMR_05h_fem_gostat_MF$time_point <- "0.5h"
DMR_1h_fem_gostat_MF$time_point <- "1h"
DMR_4h_fem_gostat_MF$time_point <- "4h"
DMR_24h_fem_gostat_MF$time_point <- "24h"
DMR_72h_fem_gostat_MF$time_point <- "72h"

DMR_05h_fem_gostat_CC$time_point <- "0.5h"
DMR_1h_fem_gostat_CC$time_point <- "1h"
DMR_4h_fem_gostat_CC$time_point <- "4h"
DMR_24h_fem_gostat_CC$time_point <- "24h"
DMR_72h_fem_gostat_CC$time_point <- "72h"

DMR_fem_all_BP <- rbind(DMR_05h_fem_gostat_BP,DMR_1h_fem_gostat_BP,DMR_4h_fem_gostat_BP,DMR_24h_fem_gostat_BP,DMR_72h_fem_gostat_BP)
DMR_fem_all_MF <- rbind(DMR_05h_fem_gostat_MF,DMR_1h_fem_gostat_MF,DMR_4h_fem_gostat_MF,DMR_24h_fem_gostat_MF,DMR_72h_fem_gostat_MF)
DMR_fem_all_CC <- rbind(DMR_05h_fem_gostat_CC,DMR_1h_fem_gostat_CC,DMR_4h_fem_gostat_CC,DMR_24h_fem_gostat_CC,DMR_72h_fem_gostat_CC)

write.csv(DMR_fem_all_BP, "./shortTerm_exp/data/go_res/DMR_fem_combined_gostat_BP.csv", row.names = FALSE)
write.csv(DMR_fem_all_MF, "./shortTerm_exp/data/go_res/DMR_fem_combined_gostat_MF.csv", row.names = FALSE)
write.csv(DMR_fem_all_CC, "./shortTerm_exp/data/go_res/DMR_fem_combined_gostat_CC.csv", row.names = FALSE)

DMR_05h_mal_gostat_BP$time_point <- "0.5h"
DMR_1h_mal_gostat_BP$time_point <- "1h"
DMR_4h_mal_gostat_BP$time_point <- "4h"
DMR_24h_mal_gostat_BP$time_point <- "24h"
DMR_72h_mal_gostat_BP$time_point <- "72h"

DMR_05h_mal_gostat_MF$time_point <- "0.5h"
DMR_1h_mal_gostat_MF$time_point <- "1h"
DMR_4h_mal_gostat_MF$time_point <- "4h"
DMR_24h_mal_gostat_MF$time_point <- "24h"
DMR_72h_mal_gostat_MF$time_point <- "72h"

DMR_05h_mal_gostat_CC$time_point <- "0.5h"
DMR_1h_mal_gostat_CC$time_point <- "1h"
DMR_4h_mal_gostat_CC$time_point <- "4h"
DMR_24h_mal_gostat_CC$time_point <- "24h"
DMR_72h_mal_gostat_CC$time_point <- "72h"

DMR_mal_all_BP <- rbind(DMR_05h_mal_gostat_BP,DMR_1h_mal_gostat_BP,DMR_4h_mal_gostat_BP,DMR_24h_mal_gostat_BP,DMR_72h_mal_gostat_BP)
DMR_mal_all_MF <- rbind(DMR_05h_mal_gostat_MF,DMR_1h_mal_gostat_MF,DMR_4h_mal_gostat_MF,DMR_24h_mal_gostat_MF,DMR_72h_mal_gostat_MF)
DMR_mal_all_CC <- rbind(DMR_05h_mal_gostat_CC,DMR_1h_mal_gostat_CC,DMR_4h_mal_gostat_CC,DMR_24h_mal_gostat_CC,DMR_72h_mal_gostat_CC)

write.csv(DMR_mal_all_BP, "./shortTerm_exp/data/go_res/DMR_mal_combined_gostat_BP.csv", row.names = FALSE)
write.csv(DMR_mal_all_MF, "./shortTerm_exp/data/go_res/DMR_mal_combined_gostat_MF.csv", row.names = FALSE)
write.csv(DMR_mal_all_CC, "./shortTerm_exp/data/go_res/DMR_mal_combined_gostat_CC.csv", row.names = FALSE)

#plot go seq analysis 
plot.funct <- function(data.05h, data.1h, data.4h, data.24h, data.72h, 
                       total.05h, total.1h, total.4h, total.24h, total.72h, 
                       name.plot){
  #calculate percent
  data.05h$perc <- data.05h$Count / total.05h * 100
  data.1h$perc <- data.1h$Count / total.1h * 100
  data.4h$perc <- data.4h$Count / total.4h * 100
  data.24h$perc <- data.24h$Count / total.24h * 100
  data.72h$perc <- data.72h$Count / total.72h * 100
  
  #select top 10 by p value and percent
  data.05h <- data.05h %>% slice_min(tibble(FDR, perc), n = 10, with_ties = FALSE) 
  data.1h <- data.1h %>% slice_min(tibble(FDR, perc), n = 10, with_ties = FALSE) 
  data.4h <- data.4h %>% slice_min(tibble(FDR, perc), n = 10, with_ties = FALSE) 
  data.24h <- data.24h %>% slice_min(tibble(FDR, perc), n = 10, with_ties = FALSE) 
  data.72h <- data.72h %>% slice_min(tibble(FDR, perc), n = 10, with_ties = FALSE) 
  
  #add time point variable 
  data.05h$tp <- "0.5h"
  data.1h$tp <- "1h"
  data.4h$tp <- "4h"
  data.24h$tp <- "24h"
  data.72h$tp <- "72h"
  
  #bind all data
  data <- rbind(data.05h, data.1h, data.4h, data.24h, data.72h)
  
  #convert time point to ordered factor
  data$tp <- factor(data$tp, levels = c("0.5h", "1h", "4h", "24h", "72h"), ordered = TRUE)
  
  #plot 
  p <- data %>% ggplot(aes(x = tp, y = Term, size = perc, color = tp, alpha = FDR)) + geom_point() + theme_bw() +
    scale_color_manual(values = c("#B8DE29FF", "#55C667FF", "#1F968BFF", "#39568CFF", "#481567FF")) +
    scale_size(range = c(.1,10), name = "% of total", breaks = c(5,10,25,50), limits = c(0,50)) + ylab(NULL) + xlab("Time Point") +
    scale_alpha_continuous(name = "Adj. P value", limits = c(0,0.05), range = c(0.8,0.1), 
                           breaks = c(0.0001, 0.001, 0.01, 0.05), labels = c("0.0001", 0.001, 0.01, 0.05)) + guides(color="none") + 
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12), legend.text = element_text(size = 12),
          legend.title = element_text(size = 14)) 
  
  ggsave(filename = name.plot, plot = p, width = 10, height = 6, units = "in", dpi = 300)
  
  return(p)
}

DMS.fem.BP <- plot.funct(DMS_05h_fem_gostat_BP, DMS_1h_fem_gostat_BP, DMS_4h_fem_gostat_BP, DMS_24h_fem_gostat_BP, DMS_72h_fem_gostat_BP,
               nrow(DMS_TSS_05h_fem), nrow(DMS_TSS_1h_fem), nrow(DMS_TSS_4h_fem), nrow(DMS_TSS_24h_fem), nrow(DMS_TSS_72h_fem),
               "./shortTerm_exp/plots/finalized_tiff/go_plots/DMS_fem_goplot_BP.tiff")

DMS.fem.MF <- plot.funct(DMS_05h_fem_gostat_MF, DMS_1h_fem_gostat_MF, DMS_4h_fem_gostat_MF, DMS_24h_fem_gostat_MF, DMS_72h_fem_gostat_MF,
           nrow(DMS_TSS_05h_fem), nrow(DMS_TSS_1h_fem), nrow(DMS_TSS_4h_fem), nrow(DMS_TSS_24h_fem), nrow(DMS_TSS_72h_fem),
           "./shortTerm_exp/plots/finalized_tiff/go_plots/DMS_fem_goplot_MF.tiff")

DMR.fem.BP <- plot.funct(DMR_05h_fem_gostat_BP, DMR_1h_fem_gostat_BP, DMR_4h_fem_gostat_BP, DMR_24h_fem_gostat_BP, DMR_72h_fem_gostat_BP,
           nrow(DMR_TSS_05h_fem), nrow(DMR_TSS_1h_fem), nrow(DMR_TSS_4h_fem), nrow(DMR_TSS_24h_fem), nrow(DMR_TSS_72h_fem),
           "./shortTerm_exp/plots/finalized_tiff/go_plots/DMR_fem_goplot_BP.tiff")

DMR.fem.MF <- plot.funct(DMR_05h_fem_gostat_MF, DMR_1h_fem_gostat_MF, DMR_4h_fem_gostat_MF, DMR_24h_fem_gostat_MF, DMR_72h_fem_gostat_MF,
           nrow(DMR_TSS_05h_fem), nrow(DMR_TSS_1h_fem), nrow(DMR_TSS_4h_fem), nrow(DMR_TSS_24h_fem), nrow(DMR_TSS_72h_fem),
           "./shortTerm_exp/plots/finalized_tiff/go_plots/DMR_fem_goplot_MF.tiff")

DMS.mal.BP <- plot.funct(DMS_05h_mal_gostat_BP, DMS_1h_mal_gostat_BP, DMS_4h_mal_gostat_BP, DMS_24h_mal_gostat_BP, DMS_72h_mal_gostat_BP,
           nrow(DMS_TSS_05h_mal), nrow(DMS_TSS_1h_mal), nrow(DMS_TSS_4h_mal), nrow(DMS_TSS_24h_mal), nrow(DMS_TSS_72h_mal),
           "./shortTerm_exp/plots/finalized_tiff/go_plots/DMS_mal_goplot_BP.tiff")

DMS.mal.MF <- plot.funct(DMS_05h_mal_gostat_MF, DMS_1h_mal_gostat_MF, DMS_4h_mal_gostat_MF, DMS_24h_mal_gostat_MF, DMS_72h_mal_gostat_MF,
           nrow(DMS_TSS_05h_mal), nrow(DMS_TSS_1h_mal), nrow(DMS_TSS_4h_mal), nrow(DMS_TSS_24h_mal), nrow(DMS_TSS_72h_mal),
           "./shortTerm_exp/plots/finalized_tiff/go_plots/DMS_mal_goplot_MF.tiff")

DMR.mal.BP <- plot.funct(DMR_05h_mal_gostat_BP, DMR_1h_mal_gostat_BP, DMR_4h_mal_gostat_BP, DMR_24h_mal_gostat_BP, DMR_72h_mal_gostat_BP,
           nrow(DMR_TSS_05h_mal), nrow(DMR_TSS_1h_mal), nrow(DMR_TSS_4h_mal), nrow(DMR_TSS_24h_mal), nrow(DMR_TSS_72h_mal),
           "./shortTerm_exp/plots/finalized_tiff/go_plots/DMR_mal_goplot_BP.tiff")

DMR.mal.MF <- plot.funct(DMR_05h_mal_gostat_MF, DMR_1h_mal_gostat_MF, DMR_4h_mal_gostat_MF, DMR_24h_mal_gostat_MF, DMR_72h_mal_gostat_MF,
           nrow(DMR_TSS_05h_mal), nrow(DMR_TSS_1h_mal), nrow(DMR_TSS_4h_mal), nrow(DMR_TSS_24h_mal), nrow(DMR_TSS_72h_mal),
           "./shortTerm_exp/plots/finalized_tiff/go_plots/DMR_mal_goplot_MF.tiff")

## Get lists of differentially methylated genes ## 
gene.list.funct <- function(data, name.csv){
  #limit to genes within 10kb 
  data <- subset(data, dist.to.feature >= -10000 & dist.to.feature <= 10000)
  #get list to seach in biomart 
  val <- data$feature.name
  #search biomart
  gene.data <- getBM(mart = bm, attributes = c('ensembl_transcript_id','external_gene_name','go_id', "name_1006", 
                                  "namespace_1003", "go_linkage_type"), filters = "ensembl_transcript_id", values = val)
  gene.data <- gene.data[!duplicated(gene.data),]
  write.csv(gene.data, name.csv, row.names = FALSE)
}

gene.list.funct(DMS_TSS_all, "./shortTerm_exp/data/gene_lists/DMS_allTP_allSex_genelist.csv")
gene.list.funct(DMS_TSS_05h, "./shortTerm_exp/data/gene_lists/DMS_05h_allSex_genelist.csv")
gene.list.funct(DMS_TSS_1h, "./shortTerm_exp/data/gene_lists/DMS_1h_allSex_genelist.csv")
gene.list.funct(DMS_TSS_4h, "./shortTerm_exp/data/gene_lists/DMS_4h_allSex_genelist.csv")
gene.list.funct(DMS_TSS_24h, "./shortTerm_exp/data/gene_lists/DMS_24h_allSex_genelist.csv")
gene.list.funct(DMS_TSS_72h, "./shortTerm_exp/data/gene_lists/DMS_72h_allSex_genelist.csv")

gene.list.funct(DMS_TSS_all_fem, "./shortTerm_exp/data/gene_lists/DMS_allTP_fem_genelist.csv")
gene.list.funct(DMS_TSS_05h_fem, "./shortTerm_exp/data/gene_lists/DMS_05h_fem_genelist.csv")
gene.list.funct(DMS_TSS_1h_fem, "./shortTerm_exp/data/gene_lists/DMS_1h_fem_genelist.csv")
gene.list.funct(DMS_TSS_4h_fem, "./shortTerm_exp/data/gene_lists/DMS_4h_fem_genelist.csv")
gene.list.funct(DMS_TSS_24h_fem, "./shortTerm_exp/data/gene_lists/DMS_24h_fem_genelist.csv")
gene.list.funct(DMS_TSS_72h_fem, "./shortTerm_exp/data/gene_lists/DMS_72h_fem_genelist.csv")

gene.list.funct(DMS_TSS_all_mal, "./shortTerm_exp/data/gene_lists/DMS_allTP_mal_genelist.csv")
gene.list.funct(DMS_TSS_05h_mal, "./shortTerm_exp/data/gene_lists/DMS_05h_mal_genelist.csv")
gene.list.funct(DMS_TSS_1h_mal, "./shortTerm_exp/data/gene_lists/DMS_1h_mal_genelist.csv")
gene.list.funct(DMS_TSS_4h_mal, "./shortTerm_exp/data/gene_lists/DMS_4h_mal_genelist.csv")
gene.list.funct(DMS_TSS_24h_mal, "./shortTerm_exp/data/gene_lists/DMS_24h_mal_genelist.csv")
gene.list.funct(DMS_TSS_72h_mal, "./shortTerm_exp/data/gene_lists/DMS_72h_mal_genelist.csv")

gene.list.funct(DMR_TSS_all, "./shortTerm_exp/data/gene_lists/DMR_allTP_allSex_genelist.csv")
gene.list.funct(DMR_TSS_05h, "./shortTerm_exp/data/gene_lists/DMR_05h_allSex_genelist.csv")
gene.list.funct(DMR_TSS_1h, "./shortTerm_exp/data/gene_lists/DMR_1h_allSex_genelist.csv")
gene.list.funct(DMR_TSS_4h, "./shortTerm_exp/data/gene_lists/DMR_4h_allSex_genelist.csv")
gene.list.funct(DMR_TSS_24h, "./shortTerm_exp/data/gene_lists/DMR_24h_allSex_genelist.csv")
gene.list.funct(DMR_TSS_72h, "./shortTerm_exp/data/gene_lists/DMR_72h_allSex_genelist.csv")

gene.list.funct(DMR_TSS_all_fem, "./shortTerm_exp/data/gene_lists/DMR_allTP_fem_genelist.csv")
gene.list.funct(DMR_TSS_05h_fem, "./shortTerm_exp/data/gene_lists/DMR_05h_fem_genelist.csv")
gene.list.funct(DMR_TSS_1h_fem, "./shortTerm_exp/data/gene_lists/DMR_1h_fem_genelist.csv")
gene.list.funct(DMR_TSS_4h_fem, "./shortTerm_exp/data/gene_lists/DMR_4h_fem_genelist.csv")
gene.list.funct(DMR_TSS_24h_fem, "./shortTerm_exp/data/gene_lists/DMR_24h_fem_genelist.csv")
gene.list.funct(DMR_TSS_72h_fem, "./shortTerm_exp/data/gene_lists/DMR_72h_fem_genelist.csv")

gene.list.funct(DMR_TSS_all_mal, "./shortTerm_exp/data/gene_lists/DMR_allTP_mal_genelist.csv")
gene.list.funct(DMR_TSS_05h_mal, "./shortTerm_exp/data/gene_lists/DMR_05h_mal_genelist.csv")
gene.list.funct(DMR_TSS_1h_mal, "./shortTerm_exp/data/gene_lists/DMR_1h_mal_genelist.csv")
gene.list.funct(DMR_TSS_4h_mal, "./shortTerm_exp/data/gene_lists/DMR_4h_mal_genelist.csv")
gene.list.funct(DMR_TSS_24h_mal, "./shortTerm_exp/data/gene_lists/DMR_24h_mal_genelist.csv")
gene.list.funct(DMR_TSS_72h_mal, "./shortTerm_exp/data/gene_lists/DMR_72h_mal_genelist.csv")

