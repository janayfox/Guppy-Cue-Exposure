#######################################################
### Goal: Plot and analyze annotation stats, run go seq 
### Author: Janay Fox
### R script
#######################################################

## Set up ##
#install packages
#install.packages("ggplot2")
#install.packages("dplyr")

#load packages
library(ggplot2)
library(dplyr)
library(genomation)
library(DescTools)
library(tibble)
library(ggpubr)
library(GenomicRanges)
library(methylKit)
library(GOstats)
library(biomaRt)
library(GSEABase)
library(ggpubr)
library(stringr)

#load in data
DMS_anno_num_all <- readRDS("./dev_exp/data/methylkit_res/anno_res/DMS_annStats_num_all.RDS")
DMS_anno_num_fem <- readRDS("./dev_exp/data/methylkit_res/anno_res/DMS_annStats_num_fem.RDS")
DMS_anno_num_mal <- readRDS("./dev_exp/data/methylkit_res/anno_res/DMS_annStats_num_mal.RDS")

DMR_anno_num_all <- readRDS("./dev_exp/data/methylkit_res/anno_res/DMR_annStats_num_all.RDS")
DMR_anno_num_fem <- readRDS("./dev_exp/data/methylkit_res/anno_res/DMR_annStats_num_fem.RDS")
DMR_anno_num_mal <- readRDS("./dev_exp/data/methylkit_res/anno_res/DMR_annStats_num_mal.RDS")

CpG_anno_num_all <- readRDS("./dev_exp/data/methylkit_res/anno_res/CpG_annStats_num_all.RDS")
CpG_anno_num_fem <- readRDS("./dev_exp/data/methylkit_res/anno_res/CpG_annStats_num_fem.RDS")
CpG_anno_num_mal <- readRDS("./dev_exp/data/methylkit_res/anno_res/CpG_annStats_num_mal.RDS")

regions_anno_num_all <- readRDS("./dev_exp/data/methylkit_res/anno_res/regions_annStats_num_all.RDS")
regions_anno_num_fem <- readRDS("./dev_exp/data/methylkit_res/anno_res/regions_annStats_num_fem.RDS")
regions_anno_num_mal <- readRDS("./dev_exp/data/methylkit_res/anno_res/regions_annStats_num_mal.RDS")

DMS_TSS_all <- readRDS("./dev_exp/data/methylkit_res/anno_res/DMS_TSS_all.RDS")
DMS_TSS_all_hyper <- readRDS("./dev_exp/data/methylkit_res/anno_res/DMS_TSS_all_hyper.RDS")
DMS_TSS_all_hypo <- readRDS("./dev_exp/data/methylkit_res/anno_res/DMS_TSS_all_hypo.RDS")

DMS_TSS_fem <- readRDS("./dev_exp/data/methylkit_res/anno_res/DMS_TSS_fem.RDS")
DMS_TSS_fem_hyper <- readRDS("./dev_exp/data/methylkit_res/anno_res/DMS_TSS_fem_hyper.RDS")
DMS_TSS_fem_hypo <- readRDS("./dev_exp/data/methylkit_res/anno_res/DMS_TSS_fem_hypo.RDS")

DMS_TSS_mal <- readRDS("./dev_exp/data/methylkit_res/anno_res/DMS_TSS_mal.RDS")
DMS_TSS_mal_hyper <- readRDS("./dev_exp/data/methylkit_res/anno_res/DMS_TSS_mal_hyper.RDS")
DMS_TSS_mal_hypo <- readRDS("./dev_exp/data/methylkit_res/anno_res/DMS_TSS_mal_hypo.RDS")

DMR_TSS_all <- readRDS("./dev_exp/data/methylkit_res/anno_res/DMR_TSS_all.RDS")
DMR_TSS_all_hyper <- readRDS("./dev_exp/data/methylkit_res/anno_res/DMR_TSS_all_hyper.RDS")
DMR_TSS_all_hypo <- readRDS("./dev_exp/data/methylkit_res/anno_res/DMR_TSS_all_hypo.RDS")

DMR_TSS_fem <- readRDS("./dev_exp/data/methylkit_res/anno_res/DMR_TSS_fem.RDS")
DMR_TSS_fem_hyper <- readRDS("./dev_exp/data/methylkit_res/anno_res/DMR_TSS_fem_hyper.RDS")
DMR_TSS_fem_hypo <- readRDS("./dev_exp/data/methylkit_res/anno_res/DMR_TSS_fem_hypo.RDS")

DMR_TSS_mal <- readRDS("./dev_exp/data/methylkit_res/anno_res/DMR_TSS_mal.RDS")
DMR_TSS_mal_hyper <- readRDS("./dev_exp/data/methylkit_res/anno_res/DMR_TSS_mal_hyper.RDS")
DMR_TSS_mal_hypo <- readRDS("./dev_exp/data/methylkit_res/anno_res/DMR_TSS_mal_hypo.RDS")

CpG_anno_all <- readRDS("./dev_exp/data/methylkit_res/anno_res/CpG_anno_all.RDS")
CpG_anno_fem <- readRDS("./dev_exp/data/methylkit_res/anno_res/CpG_anno_fem.RDS")
CpG_anno_mal <- readRDS("./dev_exp/data/methylkit_res/anno_res/CpG_anno_mal.RDS")

regions_anno_all <- readRDS("./dev_exp/data/methylkit_res/anno_res/regions_anno_all.RDS")
regions_anno_fem <- readRDS("./dev_exp/data/methylkit_res/anno_res/regions_anno_fem.RDS")
regions_anno_mal <- readRDS("./dev_exp/data/methylkit_res/anno_res/regions_anno_mal.RDS")

#set functions
#function that runs Gtests and post hoc tests
gtest.funct <- function(DMS.data, cpg.data){
  #run initial G test
  print("Initial G test")
  init.test <- GTest(DMS.data, p = cpg.data/(cpg.data[1]+cpg.data[2]+cpg.data[3]+cpg.data[4]))
  print(init.test)
  
  #run post hoc tests
  print("Promoters")
  prom.test <- GTest(c(DMS.data[1], (DMS.data[2]+DMS.data[3]+DMS.data[4])), 
                     p = c(cpg.data[1], (cpg.data[2]+cpg.data[3]+cpg.data[4]))/
                       (cpg.data[1]+cpg.data[2]+cpg.data[3]+cpg.data[4]))
  print(prom.test)
  
  print("Exons")
  ex.test <- GTest(c(DMS.data[2], (DMS.data[1]+DMS.data[3]+DMS.data[4])), 
                   p = c(cpg.data[2], (cpg.data[1]+cpg.data[3]+cpg.data[4]))/
                     (cpg.data[1]+cpg.data[2]+cpg.data[3]+cpg.data[4]))
  print(ex.test)
  
  print("Introns")
  int.test <- GTest(c(DMS.data[3], (DMS.data[1]+DMS.data[2]+DMS.data[4])), 
                    p = c(cpg.data[3], (cpg.data[1]+cpg.data[2]+cpg.data[4]))/
                      (cpg.data[1]+cpg.data[2]+cpg.data[3]+cpg.data[4]))
  print(int.test)
  
  print("Intergenic")
  inter.test <- GTest(c(DMS.data[4], (DMS.data[1]+DMS.data[2]+DMS.data[3])), 
                      p = c(cpg.data[4], (cpg.data[1]+cpg.data[2]+cpg.data[3]))/
                        (cpg.data[1]+cpg.data[2]+cpg.data[3]+cpg.data[4]))
  print(inter.test)
  
  #adjsut p values 
  gtest.p <- c(prom.test$p.value,ex.test$p.value,int.test$p.value,inter.test$p.value)
  p.adjust(gtest.p, method = "hommel")
}

#function that plots annotation distribution 
plot.dist <- function(DMS.data, DMR.data, cpg.data, plotname){
  #format data for plotting 
  DMS.plot.data <- data.frame(DMS.data)
  DMS.plot.data$dist <- "DMS"
  colnames(DMS.plot.data)[1] <- "count"
  DMS.plot.data <- rownames_to_column(DMS.plot.data, "feature")
  
  DMR.plot.data <- data.frame(DMR.data)
  DMR.plot.data$dist <- "DMR"
  colnames(DMR.plot.data)[1] <- "count"
  DMR.plot.data <- rownames_to_column(DMR.plot.data, "feature")
  
  cpg.plot.data <- data.frame(cpg.data)
  cpg.plot.data$dist <- "CpG"
  colnames(cpg.plot.data)[1] <- "count"
  cpg.plot.data <- rownames_to_column(cpg.plot.data, "feature")
  
  plot.data <- rbind(DMS.plot.data, DMR.plot.data, cpg.plot.data)
  
  p <- ggplot(plot.data, aes(x = dist, y = count, fill = feature)) +
    geom_col(colour = "black", position = "fill") + scale_y_continuous(labels = scales::percent) + theme_bw() + xlab(NULL) +
    ylab("Percent") + scale_fill_manual(values = c("exon" = "lightgreen", "intron" = "skyblue", 
                                                   "intergenic" = "orange", "promoter" = "yellow"),
                                        labels = c("exon" = "Exon", "intron" = "Intron", "intergenic" = "Intergenic",
                                                   "promoter" = "Promoter"),
                                        name = "Feature") + 
    theme(axis.title = element_text(size = 16), axis.text = element_text(size = 16), legend.text = element_text(size = 14),
          legend.title = element_text(size = 15)) 
  
  ggsave(filename = plotname, plot = p, width = 5, height = 5, units = "in", dpi = 300)
  
  return(p)
}

plot.dist.v2 <- function(DMS.data, cpg.data, plotname){
  #format data for plotting 
  DMS.plot.data <- data.frame(DMS.data)
  DMS.plot.data$dist <- "DMS"
  colnames(DMS.plot.data)[1] <- "count"
  DMS.plot.data <- rownames_to_column(DMS.plot.data, "feature")
  
  cpg.plot.data <- data.frame(cpg.data)
  cpg.plot.data$dist <- "CpG"
  colnames(cpg.plot.data)[1] <- "count"
  cpg.plot.data <- rownames_to_column(cpg.plot.data, "feature")
  
  plot.data <- rbind(DMS.plot.data, cpg.plot.data)
  
  p <- ggplot(plot.data, aes(x = dist, y = count, fill = feature)) +
    geom_col(colour = "black", position = "fill") + scale_y_continuous(labels = scales::percent) + theme_bw() + xlab(NULL) +
    ylab("Percent") + scale_fill_manual(values = c("exon" = "lightgreen", "intron" = "skyblue", 
                                                   "intergenic" = "orange", "promoter" = "yellow"),
                                        labels = c("exon" = "Exon", "intron" = "Intron", "intergenic" = "Intergenic",
                                                   "promoter" = "Promoter"),
                                        name = "Feature") + 
    theme(axis.title = element_text(size = 16), axis.text = element_text(size = 16), legend.text = element_text(size = 14),
          legend.title = element_text(size = 15)) 
  
  ggsave(filename = plotname, plot = p, width = 5, height = 5, units = "in", dpi = 300)
  
  return(p)
}

#function that runs go stats
run.gostat <- function(test.data, uni.data, ont){
  #get background gene sets from all CpGs
  uni.data <- getAssociationWithTSS(uni.data)
  
  #limit to genes within 10kb 
  test.data <- subset(test.data, dist.to.feature >= -10000 & dist.to.feature <= 10000)
  
  universe <- uni.data$feature.name
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

#function that plots go seq analysis
plot.go <- function(data.all, data.fem, data.mal, 
                       total.all, total.fem, total.mal, 
                       name.plot){
  #calculate percent
  data.all$perc <- data.all$Count / total.all * 100
  data.fem$perc <- data.fem$Count / total.fem * 100
  data.mal$perc <- data.mal$Count / total.mal * 100

  #select top 10 by p value and percent
  data.all <- data.all %>% slice_min(tibble(FDR, perc), n = 10, with_ties = FALSE)
  data.fem <- data.fem %>% slice_min(tibble(FDR, perc), n = 10, with_ties = FALSE)
  data.mal <- data.mal %>% slice_min(tibble(FDR, perc), n = 10, with_ties = FALSE)
  
  #bind all data
  data <- rbind(data.all, data.fem, data.mal)
  
  #convert time point to ordered factor
  data$sex <- factor(data$sex, levels = c("all", "fem", "mal"), ordered = TRUE)
  
  #plot
  p <- data %>% ggplot(aes(x = sex, y = Term, size = perc, color = sex, alpha = FDR)) + geom_point() + theme_bw() +
    scale_color_manual(values = c("lightgreen", "#FFE17B", "skyblue"), name = "Sex", 
                       labels = c("fem" = "Females", "mal" = "Males", "all" = "All")) +
    scale_size(range = c(3,10), name = "% of total", breaks = c(5,10,25), limits = c(0,50)) + ylab(NULL) + xlab("Sex") +
    scale_alpha_continuous(name = "Adj. P value", limits = c(0,0.05), range = c(0.8,0.1),
                           breaks = c(0.0001, 0.001, 0.01, 0.05), labels = c("0.0001", 0.001, 0.01, 0.05)) + guides(color="none") +
    scale_y_discrete(labels = function(Term) str_wrap(Term, width = 100)) +
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 10), legend.text = element_text(size = 12),
          legend.title = element_text(size = 14))
  
  ggsave(filename = name.plot, plot = p, width = 8, height = 10, units = "in", dpi = 300)
  
  return(p)
}

#function that plots go seq analysis without all
plot.go.noAll <- function(data.fem, data.mal, 
                    total.fem, total.mal, 
                    name.plot){
  #calculate percent
  data.fem$perc <- data.fem$Count / total.fem * 100
  data.mal$perc <- data.mal$Count / total.mal * 100
  
  #select top 10 by p value and percent
  data.fem <- data.fem %>% slice_min(tibble(FDR, perc), n = 10, with_ties = FALSE)
  data.mal <- data.mal %>% slice_min(tibble(FDR, perc), n = 10, with_ties = FALSE)
  
  #bind all data
  data <- rbind(data.fem, data.mal)
  
  #convert time point to ordered factor
  data$sex <- factor(data$sex, levels = c("fem", "mal"), ordered = TRUE)
  
  #plot
  p <- data %>% ggplot(aes(x = sex, y = Term, size = perc, color = sex, alpha = FDR)) + geom_point() + theme_bw() +
    scale_color_manual(values = c("#FFE17B", "skyblue"), name = "Sex", 
                       labels = c("fem" = "Females", "mal" = "Males")) +
    scale_size(range = c(3,10), name = "% of total", breaks = c(5,10,25), limits = c(0,50)) + ylab(NULL) + xlab("Sex") +
    scale_alpha_continuous(name = "Adj. P value", limits = c(0,0.05), range = c(0.8,0.1),
                           breaks = c(0.0001, 0.001, 0.01, 0.05), labels = c("0.0001", 0.001, 0.01, 0.05)) + guides(color="none") +
    scale_y_discrete(labels = function(Term) str_wrap(Term, width = 100)) +
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 10), legend.text = element_text(size = 12),
          legend.title = element_text(size = 14))
  
  ggsave(filename = name.plot, plot = p, width = 8, height = 10, units = "in", dpi = 300)
  
  return(p)
}

#function that edits go stat results
edit.go.res <- function(data, category.name, type.name, direction.name, sex.name){
  #add in extra info
  data$category <- category.name
  data$type <- type.name
  data$direction <- direction.name
  data$sex <- sex.name
  
  #edit col name
  colnames(data)[1] <- "GO_ID"
  
  return(data)
  }

## Plot annotation stats ##
dist_all <- plot.dist(DMS_anno_num_all, DMR_anno_num_all, CpG_anno_num_all, 
                      "./dev_exp/plots/finalized_tiff/anno_barplots/all_anno_barplot.tiff")

dist_fem <- plot.dist(DMS_anno_num_fem, DMR_anno_num_fem, CpG_anno_num_fem, 
                      "./dev_exp/plots/finalized_tiff/anno_barplots/fem_anno_barplot.tiff")

dist_mal <- plot.dist(DMS_anno_num_mal, DMR_anno_num_mal, CpG_anno_num_mal, 
                      "./dev_exp/plots/finalized_tiff/anno_barplots/mal_anno_barplot.tiff")

#plot without DMRs 
dist_fem <- plot.dist.v2(DMS_anno_num_fem, CpG_anno_num_fem, 
                      "./dev_exp/plots/finalized_tiff/anno_barplots/fem_anno_barplot_v2.tiff")

dist_mal <- plot.dist.v2(DMS_anno_num_mal, CpG_anno_num_mal, 
                      "./dev_exp/plots/finalized_tiff/anno_barplots/mal_anno_barplot_v2.tiff")

#run g tests
gtest.funct(DMS_anno_num_all,CpG_anno_num_all)
gtest.funct(DMS_anno_num_fem,CpG_anno_num_fem)
gtest.funct(DMS_anno_num_mal,CpG_anno_num_mal)

gtest.funct(DMR_anno_num_all,regions_anno_num_all)
gtest.funct(DMR_anno_num_fem,regions_anno_num_fem)
gtest.funct(DMR_anno_num_mal,regions_anno_num_mal)

## Run GO stats ##
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

#run on all DMS gene sets
DMS_all_gostat_BP <- run.gostat(DMS_TSS_all, CpG_anno_all, "BP")
DMS_all_hyper_gostat_BP <- run.gostat(DMS_TSS_all_hyper, CpG_anno_all, "BP")
DMS_all_hypo_gostat_BP <- run.gostat(DMS_TSS_all_hypo, CpG_anno_all, "BP")

DMS_all_gostat_MF <- run.gostat(DMS_TSS_all, CpG_anno_all, "MF")
DMS_all_hyper_gostat_MF <- run.gostat(DMS_TSS_all_hyper, CpG_anno_all, "MF")
DMS_all_hypo_gostat_MF <- run.gostat(DMS_TSS_all_hypo, CpG_anno_all, "MF")

DMS_all_gostat_CC <- run.gostat(DMS_TSS_all, CpG_anno_all, "CC")
DMS_all_hyper_gostat_CC <- run.gostat(DMS_TSS_all_hyper, CpG_anno_all, "CC")
DMS_all_hypo_gostat_CC <- run.gostat(DMS_TSS_all_hypo, CpG_anno_all, "CC")

DMS_fem_gostat_BP <- run.gostat(DMS_TSS_fem, CpG_anno_fem, "BP")
DMS_fem_hyper_gostat_BP <- run.gostat(DMS_TSS_fem_hyper, CpG_anno_fem, "BP")
DMS_fem_hypo_gostat_BP <- run.gostat(DMS_TSS_fem_hypo, CpG_anno_fem, "BP")

DMS_fem_gostat_MF <- run.gostat(DMS_TSS_fem, CpG_anno_fem, "MF")
DMS_fem_hyper_gostat_MF <- run.gostat(DMS_TSS_fem_hyper, CpG_anno_fem, "MF")
DMS_fem_hypo_gostat_MF <- run.gostat(DMS_TSS_fem_hypo, CpG_anno_fem, "MF")

DMS_fem_gostat_CC <- run.gostat(DMS_TSS_fem, CpG_anno_fem, "CC")
DMS_fem_hyper_gostat_CC <- run.gostat(DMS_TSS_fem_hyper, CpG_anno_fem, "CC")
DMS_fem_hypo_gostat_CC <- run.gostat(DMS_TSS_fem_hypo, CpG_anno_fem, "CC")

DMS_mal_gostat_BP <- run.gostat(DMS_TSS_mal, CpG_anno_mal, "BP")
DMS_mal_hyper_gostat_BP <- run.gostat(DMS_TSS_mal_hyper, CpG_anno_mal, "BP")
DMS_mal_hypo_gostat_BP <- run.gostat(DMS_TSS_mal_hypo, CpG_anno_mal, "BP")

DMS_mal_gostat_MF <- run.gostat(DMS_TSS_mal, CpG_anno_mal, "MF")
DMS_mal_hyper_gostat_MF <- run.gostat(DMS_TSS_mal_hyper, CpG_anno_mal, "MF")
DMS_mal_hypo_gostat_MF <- run.gostat(DMS_TSS_mal_hypo, CpG_anno_mal, "MF")

DMS_mal_gostat_CC <- run.gostat(DMS_TSS_mal, CpG_anno_mal, "CC")
DMS_mal_hyper_gostat_CC <- run.gostat(DMS_TSS_mal_hyper, CpG_anno_mal, "CC")
DMS_mal_hypo_gostat_CC <- run.gostat(DMS_TSS_mal_hypo, CpG_anno_mal, "CC")

#run on all DMR gene sets
DMR_all_gostat_BP <- run.gostat(DMR_TSS_all, regions_anno_all, "BP")
DMR_all_hyper_gostat_BP <- run.gostat(DMR_TSS_all_hyper, regions_anno_all, "BP")
DMR_all_hypo_gostat_BP <- run.gostat(DMR_TSS_all_hypo, regions_anno_all, "BP")

DMR_all_gostat_MF <- run.gostat(DMR_TSS_all, regions_anno_all, "MF")
DMR_all_hyper_gostat_MF <- run.gostat(DMR_TSS_all_hyper, regions_anno_all, "MF")
DMR_all_hypo_gostat_MF <- run.gostat(DMR_TSS_all_hypo, regions_anno_all, "MF")

DMR_all_gostat_CC <- run.gostat(DMR_TSS_all, regions_anno_all, "CC")
DMR_all_hyper_gostat_CC <- run.gostat(DMR_TSS_all_hyper, regions_anno_all, "CC")
DMR_all_hypo_gostat_CC <- run.gostat(DMR_TSS_all_hypo, regions_anno_all, "CC")

DMR_fem_gostat_BP <- run.gostat(DMR_TSS_fem, regions_anno_fem, "BP")
DMR_fem_hyper_gostat_BP <- run.gostat(DMR_TSS_fem_hyper, regions_anno_fem, "BP")
DMR_fem_hypo_gostat_BP <- run.gostat(DMR_TSS_fem_hypo, regions_anno_fem, "BP")

DMR_fem_gostat_MF <- run.gostat(DMR_TSS_fem, regions_anno_fem, "MF")
DMR_fem_hyper_gostat_MF <- run.gostat(DMR_TSS_fem_hyper, regions_anno_fem, "MF")
DMR_fem_hypo_gostat_MF <- run.gostat(DMR_TSS_fem_hypo, regions_anno_fem, "MF")

DMR_fem_gostat_CC <- run.gostat(DMR_TSS_fem, regions_anno_fem, "CC")
DMR_fem_hyper_gostat_CC <- run.gostat(DMR_TSS_fem_hyper, regions_anno_fem, "CC")
DMR_fem_hypo_gostat_CC <- run.gostat(DMR_TSS_fem_hypo, regions_anno_fem, "CC")

DMR_mal_gostat_BP <- run.gostat(DMR_TSS_mal, regions_anno_mal, "BP")
DMR_mal_hyper_gostat_BP <- run.gostat(DMR_TSS_mal_hyper, regions_anno_mal, "BP")
DMR_mal_hypo_gostat_BP <- run.gostat(DMR_TSS_mal_hypo, regions_anno_mal, "BP")

DMR_mal_gostat_MF <- run.gostat(DMR_TSS_mal, regions_anno_mal, "MF")
DMR_mal_hyper_gostat_MF <- run.gostat(DMR_TSS_mal_hyper, regions_anno_mal, "MF")
DMR_mal_hypo_gostat_MF <- run.gostat(DMR_TSS_mal_hypo, regions_anno_mal, "MF")

DMR_mal_gostat_CC <- run.gostat(DMR_TSS_mal, regions_anno_mal, "CC")
DMR_mal_hyper_gostat_CC <- run.gostat(DMR_TSS_mal_hyper, regions_anno_mal, "CC")
DMR_mal_hypo_gostat_CC <- run.gostat(DMR_TSS_mal_hypo, regions_anno_mal, "CC")

#add extra data
DMS_all_gostat_BP <- edit.go.res(DMS_all_gostat_BP, "BP", "DMS", "both", "all")
DMS_all_hyper_gostat_BP <- edit.go.res(DMS_all_hyper_gostat_BP, "BP", "DMS", "hyper", "all")
DMS_all_hypo_gostat_BP <- edit.go.res(DMS_all_hypo_gostat_BP, "BP", "DMS", "hypo", "all")

DMS_fem_gostat_BP <- edit.go.res(DMS_fem_gostat_BP, "BP", "DMS", "both", "fem")
DMS_fem_hyper_gostat_BP <- edit.go.res(DMS_fem_hyper_gostat_BP, "BP", "DMS", "hyper", "fem")
DMS_fem_hypo_gostat_BP <- edit.go.res(DMS_fem_hypo_gostat_BP, "BP", "DMS", "hypo", "fem")

DMS_mal_gostat_BP <- edit.go.res(DMS_mal_gostat_BP, "BP", "DMS", "both", "mal")
DMS_mal_hyper_gostat_BP <- edit.go.res(DMS_mal_hyper_gostat_BP, "BP", "DMS", "hyper", "mal")
DMS_mal_hypo_gostat_BP <- edit.go.res(DMS_mal_hypo_gostat_BP, "BP", "DMS", "hypo", "mal")

DMS_all_gostat_CC <- edit.go.res(DMS_all_gostat_CC, "CC", "DMS", "both", "all")
DMS_all_hyper_gostat_CC <- edit.go.res(DMS_all_hyper_gostat_CC, "CC", "DMS", "hyper", "all")
DMS_all_hypo_gostat_CC <- edit.go.res(DMS_all_hypo_gostat_CC, "CC", "DMS", "hypo", "all")

DMS_fem_gostat_CC <- edit.go.res(DMS_fem_gostat_CC, "CC", "DMS", "both", "fem")
DMS_fem_hyper_gostat_CC <- edit.go.res(DMS_fem_hyper_gostat_CC, "CC", "DMS", "hyper", "fem")
DMS_fem_hypo_gostat_CC <- edit.go.res(DMS_fem_hypo_gostat_CC, "CC", "DMS", "hypo", "fem")

DMS_mal_gostat_CC <- edit.go.res(DMS_mal_gostat_CC, "CC", "DMS", "both", "mal")
DMS_mal_hyper_gostat_CC <- edit.go.res(DMS_mal_hyper_gostat_CC, "CC", "DMS", "hyper", "mal")
DMS_mal_hypo_gostat_CC <- edit.go.res(DMS_mal_hypo_gostat_CC, "CC", "DMS", "hypo", "mal")

DMS_all_gostat_MF <- edit.go.res(DMS_all_gostat_MF, "MF", "DMS", "both", "all")
DMS_all_hyper_gostat_MF <- edit.go.res(DMS_all_hyper_gostat_MF, "MF", "DMS", "hyper", "all")
DMS_all_hypo_gostat_MF <- edit.go.res(DMS_all_hypo_gostat_MF, "MF", "DMS", "hypo", "all")

DMS_fem_gostat_MF <- edit.go.res(DMS_fem_gostat_MF, "MF", "DMS", "both", "fem")
DMS_fem_hyper_gostat_MF <- edit.go.res(DMS_fem_hyper_gostat_MF, "MF", "DMS", "hyper", "fem")
DMS_fem_hypo_gostat_MF <- edit.go.res(DMS_fem_hypo_gostat_MF, "MF", "DMS", "hypo", "fem")

DMS_mal_gostat_MF <- edit.go.res(DMS_mal_gostat_MF, "MF", "DMS", "both", "mal")
DMS_mal_hyper_gostat_MF <- edit.go.res(DMS_mal_hyper_gostat_MF, "MF", "DMS", "hyper", "mal")
DMS_mal_hypo_gostat_MF <- edit.go.res(DMS_mal_hypo_gostat_MF, "MF", "DMS", "hypo", "mal")

DMR_all_gostat_BP <- edit.go.res(DMR_all_gostat_BP, "BP", "DMR", "both", "all")
DMR_all_hyper_gostat_BP <- edit.go.res(DMR_all_hyper_gostat_BP, "BP", "DMR", "hyper", "all")
DMR_all_hypo_gostat_BP <- edit.go.res(DMR_all_hypo_gostat_BP, "BP", "DMR", "hypo", "all")

DMR_fem_gostat_BP <- edit.go.res(DMR_fem_gostat_BP, "BP", "DMR", "both", "fem")
DMR_fem_hyper_gostat_BP <- edit.go.res(DMR_fem_hyper_gostat_BP, "BP", "DMR", "hyper", "fem")
DMR_fem_hypo_gostat_BP <- edit.go.res(DMR_fem_hypo_gostat_BP, "BP", "DMR", "hypo", "fem")

DMR_mal_gostat_BP <- edit.go.res(DMR_mal_gostat_BP, "BP", "DMR", "both", "mal")
DMR_mal_hyper_gostat_BP <- edit.go.res(DMR_mal_hyper_gostat_BP, "BP", "DMR", "hyper", "mal")
DMR_mal_hypo_gostat_BP <- edit.go.res(DMR_mal_hypo_gostat_BP, "BP", "DMR", "hypo", "mal")

DMR_all_gostat_CC <- edit.go.res(DMR_all_gostat_CC, "CC", "DMR", "both", "all")
DMR_all_hyper_gostat_CC <- edit.go.res(DMR_all_hyper_gostat_CC, "CC", "DMR", "hyper", "all")
DMR_all_hypo_gostat_CC <- edit.go.res(DMR_all_hypo_gostat_CC, "CC", "DMR", "hypo", "all")

DMR_fem_gostat_CC <- edit.go.res(DMR_fem_gostat_CC, "CC", "DMR", "both", "fem")
DMR_fem_hyper_gostat_CC <- edit.go.res(DMR_fem_hyper_gostat_CC, "CC", "DMR", "hyper", "fem")
DMR_fem_hypo_gostat_CC <- edit.go.res(DMR_fem_hypo_gostat_CC, "CC", "DMR", "hypo", "fem")

DMR_mal_gostat_CC <- edit.go.res(DMR_mal_gostat_CC, "CC", "DMR", "both", "mal")
DMR_mal_hyper_gostat_CC <- edit.go.res(DMR_mal_hyper_gostat_CC, "CC", "DMR", "hyper", "mal")
DMR_mal_hypo_gostat_CC <- edit.go.res(DMR_mal_hypo_gostat_CC, "CC", "DMR", "hypo", "mal")

DMR_all_gostat_MF <- edit.go.res(DMR_all_gostat_MF, "MF", "DMR", "both", "all")
DMR_all_hyper_gostat_MF <- edit.go.res(DMR_all_hyper_gostat_MF, "MF", "DMR", "hyper", "all")
DMR_all_hypo_gostat_MF <- edit.go.res(DMR_all_hypo_gostat_MF, "MF", "DMR", "hypo", "all")

DMR_fem_gostat_MF <- edit.go.res(DMR_fem_gostat_MF, "MF", "DMR", "both", "fem")
DMR_fem_hyper_gostat_MF <- edit.go.res(DMR_fem_hyper_gostat_MF, "MF", "DMR", "hyper", "fem")
DMR_fem_hypo_gostat_MF <- edit.go.res(DMR_fem_hypo_gostat_MF, "MF", "DMR", "hypo", "fem")

DMR_mal_gostat_MF <- edit.go.res(DMR_mal_gostat_MF, "MF", "DMR", "both", "mal")
DMR_mal_hyper_gostat_MF <- edit.go.res(DMR_mal_hyper_gostat_MF, "MF", "DMR", "hyper", "mal")
DMR_mal_hypo_gostat_MF <- edit.go.res(DMR_mal_hypo_gostat_MF, "MF", "DMR", "hypo", "mal")

#combine datasets and save
all_DMS_data <- rbind(DMS_all_gostat_BP, DMS_all_hyper_gostat_BP, DMS_all_hypo_gostat_BP,
                  DMS_all_gostat_CC, DMS_all_hyper_gostat_CC, DMS_all_hypo_gostat_CC,
                  DMS_all_gostat_MF, DMS_all_hyper_gostat_MF, DMS_all_hypo_gostat_MF)
write.csv(all_DMS_data, "./dev_exp/data/go_res/go_res_DMS_all.csv", row.names = FALSE)

all_DMR_data <- rbind(DMR_all_gostat_BP, DMR_all_hyper_gostat_BP, DMR_all_hypo_gostat_BP,
                      DMR_all_gostat_CC, DMR_all_hyper_gostat_CC, DMR_all_hypo_gostat_CC,
                      DMR_all_gostat_MF, DMR_all_hyper_gostat_MF, DMR_all_hypo_gostat_MF)
write.csv(all_DMR_data, "./dev_exp/data/go_res/go_res_DMR_all.csv", row.names = FALSE)

fem_DMS_data <- rbind(DMS_fem_gostat_BP, DMS_fem_hyper_gostat_BP, DMS_fem_hypo_gostat_BP,
                      DMS_fem_gostat_CC, DMS_fem_hyper_gostat_CC, DMS_fem_hypo_gostat_CC,
                      DMS_fem_gostat_MF, DMS_fem_hyper_gostat_MF, DMS_fem_hypo_gostat_MF)
write.csv(fem_DMS_data, "./dev_exp/data/go_res/go_res_DMS_fem.csv", row.names = FALSE)

fem_DMR_data <- rbind(DMR_fem_gostat_BP, DMR_fem_hyper_gostat_BP, DMR_fem_hypo_gostat_BP,
                      DMR_fem_gostat_CC, DMR_fem_hyper_gostat_CC, DMR_fem_hypo_gostat_CC,
                      DMR_fem_gostat_MF, DMR_fem_hyper_gostat_MF, DMR_fem_hypo_gostat_MF)
write.csv(fem_DMR_data, "./dev_exp/data/go_res/go_res_DMR_fem.csv", row.names = FALSE)

mal_DMS_data <- rbind(DMS_mal_gostat_BP, DMS_mal_hyper_gostat_BP, DMS_mal_hypo_gostat_BP,
                      DMS_mal_gostat_CC, DMS_mal_hyper_gostat_CC, DMS_mal_hypo_gostat_CC,
                      DMS_mal_gostat_MF, DMS_mal_hyper_gostat_MF, DMS_mal_hypo_gostat_MF)
write.csv(mal_DMS_data, "./dev_exp/data/go_res/go_res_DMS_mal.csv", row.names = FALSE)

mal_DMR_data <- rbind(DMR_mal_gostat_BP, DMR_mal_hyper_gostat_BP, DMR_mal_hypo_gostat_BP,
                      DMR_mal_gostat_CC, DMR_mal_hyper_gostat_CC, DMR_mal_hypo_gostat_CC,
                      DMR_mal_gostat_MF, DMR_mal_hyper_gostat_MF, DMR_mal_hypo_gostat_MF)
write.csv(mal_DMR_data, "./dev_exp/data/go_res/go_res_DMR_mal.csv", row.names = FALSE)

## Plot go seq results ##
DMS.BP.hyper.plot <- plot.go(DMS_all_hyper_gostat_BP, DMS_fem_hyper_gostat_BP, DMS_mal_hyper_gostat_BP,
                                nrow(DMS_TSS_all_hyper), nrow(DMS_TSS_fem_hyper), nrow(DMS_TSS_mal_hyper),
                                "./dev_exp/plots/finalized_tiff/go_seq/DMS_hyper_BP_go_plot.tiff")

DMS.BP.hypo.plot <- plot.go(DMS_all_hypo_gostat_BP, DMS_fem_hypo_gostat_BP, DMS_mal_hypo_gostat_BP,
                             nrow(DMS_TSS_all_hypo), nrow(DMS_TSS_fem_hypo), nrow(DMS_TSS_mal_hypo),
                             "./dev_exp/plots/finalized_tiff/go_seq/DMS_hypo_BP_go_plot.tiff")

DMR.BP.hyper.plot <- plot.go(DMR_all_hyper_gostat_BP, DMR_fem_hyper_gostat_BP, DMR_mal_hyper_gostat_BP,
                             nrow(DMR_TSS_all_hyper), nrow(DMR_TSS_fem_hyper), nrow(DMR_TSS_mal_hyper),
                             "./dev_exp/plots/finalized_tiff/go_seq/DMR_hyper_BP_go_plot.tiff")

DMR.BP.hypo.plot <- plot.go(DMR_all_hypo_gostat_BP, DMR_fem_hypo_gostat_BP, DMR_mal_hypo_gostat_BP,
                            nrow(DMR_TSS_all_hypo), nrow(DMR_TSS_fem_hypo), nrow(DMR_TSS_mal_hypo),
                            "./dev_exp/plots/finalized_tiff/go_seq/DMR_hypo_BP_go_plot.tiff")

DMS.MF.hyper.plot <- plot.go(DMS_all_hyper_gostat_MF, DMS_fem_hyper_gostat_MF, DMS_mal_hyper_gostat_MF,
                             nrow(DMS_TSS_all_hyper), nrow(DMS_TSS_fem_hyper), nrow(DMS_TSS_mal_hyper),
                             "./dev_exp/plots/finalized_tiff/go_seq/DMS_hyper_MF_go_plot.tiff")

DMS.MF.hypo.plot <- plot.go(DMS_all_hypo_gostat_MF, DMS_fem_hypo_gostat_MF, DMS_mal_hypo_gostat_MF,
                            nrow(DMS_TSS_all_hypo), nrow(DMS_TSS_fem_hypo), nrow(DMS_TSS_mal_hypo),
                            "./dev_exp/plots/finalized_tiff/go_seq/DMS_hypo_MF_go_plot.tiff")

DMR.MF.hyper.plot <- plot.go(DMR_all_hyper_gostat_MF, DMR_fem_hyper_gostat_MF, DMR_mal_hyper_gostat_MF,
                             nrow(DMR_TSS_all_hyper), nrow(DMR_TSS_fem_hyper), nrow(DMR_TSS_mal_hyper),
                             "./dev_exp/plots/finalized_tiff/go_seq/DMR_hyper_MF_go_plot.tiff")

DMR.MF.hypo.plot <- plot.go(DMR_all_hypo_gostat_MF, DMR_fem_hypo_gostat_MF, DMR_mal_hypo_gostat_MF,
                            nrow(DMR_TSS_all_hypo), nrow(DMR_TSS_fem_hypo), nrow(DMR_TSS_mal_hypo),
                            "./dev_exp/plots/finalized_tiff/go_seq/DMR_hypo_MF_go_plot.tiff")

DMS.CC.hyper.plot <- plot.go(DMS_all_hyper_gostat_CC, DMS_fem_hyper_gostat_CC, DMS_mal_hyper_gostat_CC,
                             nrow(DMS_TSS_all_hyper), nrow(DMS_TSS_fem_hyper), nrow(DMS_TSS_mal_hyper),
                             "./dev_exp/plots/finalized_tiff/go_seq/DMS_hyper_CC_go_plot.tiff")

DMS.CC.hypo.plot <- plot.go(DMS_all_hypo_gostat_CC, DMS_fem_hypo_gostat_CC, DMS_mal_hypo_gostat_CC,
                            nrow(DMS_TSS_all_hypo), nrow(DMS_TSS_fem_hypo), nrow(DMS_TSS_mal_hypo),
                            "./dev_exp/plots/finalized_tiff/go_seq/DMS_hypo_CC_go_plot.tiff")

DMR.CC.hyper.plot <- plot.go(DMR_all_hyper_gostat_CC, DMR_fem_hyper_gostat_CC, DMR_mal_hyper_gostat_CC,
                             nrow(DMR_TSS_all_hyper), nrow(DMR_TSS_fem_hyper), nrow(DMR_TSS_mal_hyper),
                             "./dev_exp/plots/finalized_tiff/go_seq/DMR_hyper_CC_go_plot.tiff")

DMR.CC.hypo.plot <- plot.go(DMR_all_hypo_gostat_CC, DMR_fem_hypo_gostat_CC, DMR_mal_hypo_gostat_CC,
                            nrow(DMR_TSS_all_hypo), nrow(DMR_TSS_fem_hypo), nrow(DMR_TSS_mal_hypo),
                            "./dev_exp/plots/finalized_tiff/go_seq/DMR_hypo_CC_go_plot.tiff")

## Plot go seq results without all##
DMS.BP.hyper.plot <- plot.go.noAll(DMS_fem_hyper_gostat_BP, DMS_mal_hyper_gostat_BP,
                             nrow(DMS_TSS_fem_hyper), nrow(DMS_TSS_mal_hyper),
                             "./dev_exp/plots/finalized_tiff/go_seq/DMS_hyper_BP_go_plot_noAll.tiff")

DMS.BP.hypo.plot <- plot.go.noAll(DMS_fem_hypo_gostat_BP, DMS_mal_hypo_gostat_BP,
                            nrow(DMS_TSS_fem_hypo), nrow(DMS_TSS_mal_hypo),
                            "./dev_exp/plots/finalized_tiff/go_seq/DMS_hypo_BP_go_plot_noAll.tiff")

DMR.BP.hyper.plot <- plot.go.noAll(DMR_fem_hyper_gostat_BP, DMR_mal_hyper_gostat_BP,
                             nrow(DMR_TSS_fem_hyper), nrow(DMR_TSS_mal_hyper),
                             "./dev_exp/plots/finalized_tiff/go_seq/DMR_hyper_BP_go_plot_noAll.tiff")

DMR.BP.hypo.plot <- plot.go.noAll(DMR_fem_hypo_gostat_BP, DMR_mal_hypo_gostat_BP,
                            nrow(DMR_TSS_fem_hypo), nrow(DMR_TSS_mal_hypo),
                            "./dev_exp/plots/finalized_tiff/go_seq/DMR_hypo_BP_go_plot_noAll.tiff")

DMS.MF.hyper.plot <- plot.go.noAll(DMS_fem_hyper_gostat_MF, DMS_mal_hyper_gostat_MF,
                             nrow(DMS_TSS_fem_hyper), nrow(DMS_TSS_mal_hyper),
                             "./dev_exp/plots/finalized_tiff/go_seq/DMS_hyper_MF_go_plot_noAll.tiff")

DMS.MF.hypo.plot <- plot.go.noAll(DMS_fem_hypo_gostat_MF, DMS_mal_hypo_gostat_MF,
                            nrow(DMS_TSS_fem_hypo), nrow(DMS_TSS_mal_hypo),
                            "./dev_exp/plots/finalized_tiff/go_seq/DMS_hypo_MF_go_plot_noAll.tiff")

DMR.MF.hyper.plot <- plot.go.noAll(DMR_fem_hyper_gostat_MF, DMR_mal_hyper_gostat_MF,
                             nrow(DMR_TSS_fem_hyper), nrow(DMR_TSS_mal_hyper),
                             "./dev_exp/plots/finalized_tiff/go_seq/DMR_hyper_MF_go_plot_noAll.tiff")

DMR.MF.hypo.plot <- plot.go.noAll(DMR_fem_hypo_gostat_MF, DMR_mal_hypo_gostat_MF,
                            nrow(DMR_TSS_fem_hypo), nrow(DMR_TSS_mal_hypo),
                            "./dev_exp/plots/finalized_tiff/go_seq/DMR_hypo_MF_go_plot_noAll.tiff")

DMS.CC.hyper.plot <- plot.go.noAll(DMS_fem_hyper_gostat_CC, DMS_mal_hyper_gostat_CC,
                             nrow(DMS_TSS_fem_hyper), nrow(DMS_TSS_mal_hyper),
                             "./dev_exp/plots/finalized_tiff/go_seq/DMS_hyper_CC_go_plot_noAll.tiff")

DMS.CC.hypo.plot <- plot.go.noAll(DMS_fem_hypo_gostat_CC, DMS_mal_hypo_gostat_CC,
                            nrow(DMS_TSS_fem_hypo), nrow(DMS_TSS_mal_hypo),
                            "./dev_exp/plots/finalized_tiff/go_seq/DMS_hypo_CC_go_plot_noAll.tiff")

DMR.CC.hyper.plot <- plot.go.noAll(DMR_fem_hyper_gostat_CC, DMR_mal_hyper_gostat_CC,
                             nrow(DMR_TSS_fem_hyper), nrow(DMR_TSS_mal_hyper),
                             "./dev_exp/plots/finalized_tiff/go_seq/DMR_hyper_CC_go_plot_noAll.tiff")

DMR.CC.hypo.plot <- plot.go.noAll(DMR_fem_hypo_gostat_CC, DMR_mal_hypo_gostat_CC,
                            nrow(DMR_TSS_fem_hypo), nrow(DMR_TSS_mal_hypo),
                            "./dev_exp/plots/finalized_tiff/go_seq/DMR_hypo_CC_go_plot_noAll.tiff")
