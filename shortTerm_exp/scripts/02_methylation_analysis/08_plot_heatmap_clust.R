#######################################################
### Goal: Plot heatmaps and do cluster analysis
### Author: Janay Fox
### R script
#######################################################

## Set up ##
#install packages
#install.packages("ggplot2")
#install.packages("dplyr")
#install.packages("ggpubr")

#load packaes
library(ggplot2)
library(dplyr)
library(ggpubr)
library(methylKit)
library(pheatmap)
library(stringr)

#write functions
#create function that makes a column with a name for each DMS/DMR
nameSite <- function(data){
  data$site_name <- paste(data$chr, data$start, sep = "_")
  return(data)
}

#function that generates a percent methylation matrix
get_percMeth_matrix <- function(diffMethfile, methFile){
  
  #read data files
  diffMeth <- readRDS(diffMethfile)
  meth <- readRDS(methFile)
  
  #calculate percent methylation 
  perc_meth <- as.data.frame(percMethylation(meth))
  
  meth <- nameSite(meth)
  diffMeth <- nameSite(diffMeth)
  
  #add rownames to perc meth matrix  
  rownames(perc_meth) <- meth$site_name
  
  #subset only DMS from meth matrix 
  perc_meth <- perc_meth[rownames(perc_meth) %in% diffMeth$site_name, ]
  
  return(perc_meth)
  
}

#make function to plot heatmaps of individual time points
plot_heatmap <- function(data, plotname){
  #cluster samples
  sample.clust <- hclust(dist(t(data), method = "euclidean"), method = "ward.D2")
  
  #generate color pallette for heatmap
  cor.color <- colorRampPalette(c("#FFFF00", "black", "#0000FF"))(100)
  
  #make dataframe for sample annotations
  sample.info <- data.frame(sample = colnames(data))
  sample.info[substr(sample.info$sample,4,4)=="A", "Cue"] <- "Alarm Cue"
  sample.info[substr(sample.info$sample,4,4)=="C", "Cue"] <- "Control"
  rownames(sample.info) <- sample.info[,1]
  sample.info[,1] <- NULL
  
  #set annotation colors
  anno.colors <- list(Cue = c("Alarm Cue" = "#E14D2A", "Control" = "#31C6D4"))
  
  #plot heatmap
  p <-  pheatmap(data, cluster_rows = TRUE, cluster_cols = sample.clust,
                 show_colnames = FALSE, show_rownames = FALSE, scale = "row",
                 color = cor.color, annotation_col = sample.info,
                 annotation_colors = anno.colors, treeheight_row = 0)
  #save plot
  tiff(plotname, units="in", width = 6, height = 6, res = 600)
  print(p)
  dev.off()
}

#make function for plotting all time points together
plot_heatmap_combined <- function(data, plotname){
  #cluster samples
  sample.clust <- hclust(dist(t(data), method = "euclidean"), method = "ward.D2")
  
  #generate color pallette for heatmap
  cor.color <- colorRampPalette(c("#FFFF00", "black", "#0000FF"))(100)
  
  #make dataframe for sample annotations
  sample.info <- data.frame(sample = colnames(data))
  sample.info[substr(sample.info$sample,4,4)=="A", "Cue"] <- "Alarm Cue"
  sample.info[substr(sample.info$sample,4,4)=="C", "Cue"] <- "Control"
  sample.info$time_point <- "NA"
  sample.info <- within(sample.info, time_point[sample %in% IDs_05h] <- "05h")
  sample.info <- within(sample.info, time_point[sample %in% IDs_1h] <- "1h")
  sample.info <- within(sample.info, time_point[sample %in% IDs_4h] <- "4h")
  sample.info <- within(sample.info, time_point[sample %in% IDs_24h] <- "24h")
  sample.info <- within(sample.info, time_point[sample %in% IDs_72h] <- "72h")
  
  rownames(sample.info) <- sample.info[,1]
  sample.info[,1] <- NULL
  
  #set annotation colors
  anno.colors <- list(Cue = c("Alarm Cue" = "#E14D2A", "Control" = "#31C6D4"),
                      time_point = c("05h" = "#B8DE29FF", "1h" = "#55C667FF", "4h" = "#1F968BFF", "24h" = "#39568CFF", "72h" = "#481567FF"))
  
  #plot heatmap
  p <-  pheatmap(data, cluster_rows = TRUE, cluster_cols = sample.clust,
                 show_colnames = FALSE, show_rownames = FALSE, scale = "row",
                 color = cor.color, annotation_col = sample.info,
                 annotation_colors = anno.colors, treeheight_row = 0)
  #save plot
  tiff(plotname, units="in", width = 6, height = 6, res = 600)
  print(p)
  dev.off()
}
#get percent methylation marices for all datasets
DMS_fem_all_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res/DMS_res_all/DMSdiffmeth_all_fem_5X.RDS",
                                            "./shortTerm_exp/data/methylkit_res/DMS_res_all/DMSmeth_all_fem_5X.RDS")

DMS_fem_05h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res/DMS_res_05h/DMSdiffmeth_05h_fem_5X.RDS",
                                           "./shortTerm_exp/data/methylkit_res/DMS_res_05h/DMSmeth_05h_fem_5X.RDS")

DMS_fem_1h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res/DMS_res_1h/DMSdiffmeth_1h_fem_5X.RDS",
                                           "./shortTerm_exp/data/methylkit_res/DMS_res_1h/DMSmeth_1h_fem_5X.RDS")

DMS_fem_4h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res/DMS_res_4h/DMSdiffmeth_4h_fem_5X.RDS",
                                           "./shortTerm_exp/data/methylkit_res/DMS_res_4h/DMSmeth_4h_fem_5X.RDS")

DMS_fem_24h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res/DMS_res_24h/DMSdiffmeth_24h_fem_5X.RDS",
                                           "./shortTerm_exp/data/methylkit_res/DMS_res_24h/DMSmeth_24h_fem_5X.RDS")

DMS_fem_72h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res/DMS_res_72h/DMSdiffmeth_72h_fem_5X.RDS",
                                           "./shortTerm_exp/data/methylkit_res/DMS_res_72h/DMSmeth_72h_fem_5X.RDS")

DMS_mal_all_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res/DMS_res_all/DMSdiffmeth_all_mal_5X.RDS",
                                            "./shortTerm_exp/data/methylkit_res/DMS_res_all/DMSmeth_all_mal_5X.RDS")

DMS_mal_05h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res/DMS_res_05h/DMSdiffmeth_05h_mal_5X.RDS",
                                            "./shortTerm_exp/data/methylkit_res/DMS_res_05h/DMSmeth_05h_mal_5X.RDS")

DMS_mal_1h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res/DMS_res_1h/DMSdiffmeth_1h_mal_5X.RDS",
                                           "./shortTerm_exp/data/methylkit_res/DMS_res_1h/DMSmeth_1h_mal_5X.RDS")

DMS_mal_4h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res/DMS_res_4h/DMSdiffmeth_4h_mal_5X.RDS",
                                           "./shortTerm_exp/data/methylkit_res/DMS_res_4h/DMSmeth_4h_mal_5X.RDS")

DMS_mal_24h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res/DMS_res_24h/DMSdiffmeth_24h_mal_5X.RDS",
                                            "./shortTerm_exp/data/methylkit_res/DMS_res_24h/DMSmeth_24h_mal_5X.RDS")

DMS_mal_72h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res/DMS_res_72h/DMSdiffmeth_72h_mal_5X.RDS",
                                            "./shortTerm_exp/data/methylkit_res/DMS_res_72h/DMSmeth_72h_mal_5X.RDS")

DMR_fem_all_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res/DMR_res_all/DMRdiffmeth_all_fem_5X.RDS",
                                            "./shortTerm_exp/data/methylkit_res/DMR_res_all/DMRmeth_all_fem_5X.RDS")

DMR_fem_05h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res/DMR_res_05h/DMRdiffmeth_05h_fem_5X.RDS",
                                            "./shortTerm_exp/data/methylkit_res/DMR_res_05h/DMR_tile_meth_05h_fem_5X.RDS")

DMR_fem_1h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res/DMR_res_1h/DMRdiffmeth_1h_fem_5X.RDS",
                                           "./shortTerm_exp/data/methylkit_res/DMR_res_1h/DMR_tile_meth_1h_fem_5X.RDS")

DMR_fem_4h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res/DMR_res_4h/DMRdiffmeth_4h_fem_5X.RDS",
                                           "./shortTerm_exp/data/methylkit_res/DMR_res_4h/DMR_tile_meth_4h_fem_5X.RDS")

DMR_fem_24h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res/DMR_res_24h/DMRdiffmeth_24h_fem_5X.RDS",
                                            "./shortTerm_exp/data/methylkit_res/DMR_res_24h/DMR_tile_meth_24h_fem_5X.RDS")

DMR_fem_72h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res/DMR_res_72h/DMRdiffmeth_72h_fem_5X.RDS",
                                            "./shortTerm_exp/data/methylkit_res/DMR_res_72h/DMR_tile_meth_72h_fem_5X.RDS")

DMR_mal_all_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res/DMR_res_all/DMRdiffmeth_all_mal_5X.RDS",
                                            "./shortTerm_exp/data/methylkit_res/DMR_res_all/DMR_tile_meth_all_mal_5X.RDS")

DMR_mal_05h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res/DMR_res_05h/DMRdiffmeth_05h_mal_5X.RDS",
                                            "./shortTerm_exp/data/methylkit_res/DMR_res_05h/DMR_tile_meth_05h_mal_5X.RDS")

DMR_mal_1h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res/DMR_res_1h/DMRdiffmeth_1h_mal_5X.RDS",
                                           "./shortTerm_exp/data/methylkit_res/DMR_res_1h/DMR_tile_meth_1h_mal_5X.RDS")

DMR_mal_4h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res/DMR_res_4h/DMRdiffmeth_4h_mal_5X.RDS",
                                           "./shortTerm_exp/data/methylkit_res/DMR_res_4h/DMR_tile_meth_4h_mal_5X.RDS")

DMR_mal_24h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res/DMR_res_24h/DMRdiffmeth_24h_mal_5X.RDS",
                                            "./shortTerm_exp/data/methylkit_res/DMR_res_24h/DMR_tile_meth_24h_mal_5X.RDS")

DMR_mal_72h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res/DMR_res_72h/DMRdiffmeth_72h_mal_5X.RDS",
                                            "./shortTerm_exp/data/methylkit_res/DMR_res_72h/DMR_tile_meth_72h_mal_5X.RDS")

#get methylation matrix for all sites and regions (time points combined)
#read data files
meth_allsites_fem <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMS_res_all/DMSmeth_all_fem_5X.RDS")
meth_allsites_mal <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMS_res_all/DMSmeth_all_mal_5X.RDS")

meth_allregions_fem <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMR_res_all/DMR_tile_meth_all_fem_5X.RDS")
meth_allregions_mal <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMR_res_all/DMR_tile_meth_all_mal_5X.RDS")

#calculate percent methylation 
perc_meth_allsites_fem <- as.data.frame(percMethylation(meth_allsites_fem))
perc_meth_allsites_mal <- as.data.frame(percMethylation(meth_allsites_mal))

perc_meth_allregions_fem <- as.data.frame(percMethylation(meth_allregions_fem))
perc_meth_allregions_mal <- as.data.frame(percMethylation(meth_allregions_mal))

#add site names
meth_allsites_fem <- nameSite(meth_allsites_fem)
meth_allsites_mal <- nameSite(meth_allsites_mal)

meth_allregions_fem <- nameSite(meth_allregions_fem)
meth_allregions_mal <- nameSite(meth_allregions_mal)

#add rownames to perc meth matrix  
rownames(perc_meth_allsites_fem) <- meth_allsites_fem$site_name
rownames(perc_meth_allsites_mal) <- meth_allsites_mal$site_name

rownames(perc_meth_allregions_fem) <- meth_allregions_fem$site_name
rownames(perc_meth_allregions_mal) <- meth_allregions_mal$site_name

#remove extra file 
rm(meth_allsites_fem)
rm(meth_allsites_mal)

rm(meth_allregions_fem)
rm(meth_allregions_mal)

##Make heat maps ##
#plot for each time point/sex
plot_heatmap(DMS_fem_all_percMeth, "./shortTerm_exp/plots/finalized_tiff/heatmaps/allTP_DMS_fem_heatmap.tiff")
plot_heatmap(DMS_fem_05h_percMeth, "./shortTerm_exp/plots/finalized_tiff/heatmaps/05h_DMS_fem_heatmap.tiff")
plot_heatmap(DMS_fem_1h_percMeth, "./shortTerm_exp/plots/finalized_tiff/heatmaps/1h_DMS_fem_heatmap.tiff")
plot_heatmap(DMS_fem_4h_percMeth, "./shortTerm_exp/plots/finalized_tiff/heatmaps/4h_DMS_fem_heatmap.tiff")
plot_heatmap(DMS_fem_24h_percMeth, "./shortTerm_exp/plots/finalized_tiff/heatmaps/24h_DMS_fem_heatmap.tiff")
plot_heatmap(DMS_fem_72h_percMeth, "./shortTerm_exp/plots/finalized_tiff/heatmaps/72h_DMS_fem_heatmap.tiff")

plot_heatmap(DMS_mal_all_percMeth, "./shortTerm_exp/plots/finalized_tiff/heatmaps/allTP_DMS_mal_heatmap.tiff")
plot_heatmap(DMS_mal_05h_percMeth, "./shortTerm_exp/plots/finalized_tiff/heatmaps/05h_DMS_mal_heatmap.tiff")
plot_heatmap(DMS_mal_1h_percMeth, "./shortTerm_exp/plots/finalized_tiff/heatmaps/1h_DMS_mal_heatmap.tiff")
plot_heatmap(DMS_mal_4h_percMeth, "./shortTerm_exp/plots/finalized_tiff/heatmaps/4h_DMS_mal_heatmap.tiff")
plot_heatmap(DMS_mal_24h_percMeth, "./shortTerm_exp/plots/finalized_tiff/heatmaps/24h_DMS_mal_heatmap.tiff")
plot_heatmap(DMS_mal_72h_percMeth, "./shortTerm_exp/plots/finalized_tiff/heatmaps/72h_DMS_mal_heatmap.tiff")

#DMRs
plot_heatmap(DMR_fem_05h_percMeth, "./shortTerm_exp/plots/finalized_tiff/heatmaps/05h_DMR_fem_heatmap.tiff")
plot_heatmap(DMR_fem_1h_percMeth, "./shortTerm_exp/plots/finalized_tiff/heatmaps/1h_DMR_fem_heatmap.tiff")
plot_heatmap(DMR_fem_4h_percMeth, "./shortTerm_exp/plots/finalized_tiff/heatmaps/4h_DMR_fem_heatmap.tiff")
plot_heatmap(DMR_fem_24h_percMeth, "./shortTerm_exp/plots/finalized_tiff/heatmaps/24h_DMR_fem_heatmap.tiff")
plot_heatmap(DMR_fem_72h_percMeth, "./shortTerm_exp/plots/finalized_tiff/heatmaps/72h_DMR_fem_heatmap.tiff")

plot_heatmap(DMR_mal_05h_percMeth, "./shortTerm_exp/plots/finalized_tiff/heatmaps/05h_DMR_mal_heatmap.tiff")
plot_heatmap(DMR_mal_1h_percMeth, "./shortTerm_exp/plots/finalized_tiff/heatmaps/1h_DMR_mal_heatmap.tiff")
plot_heatmap(DMR_mal_4h_percMeth, "./shortTerm_exp/plots/finalized_tiff/heatmaps/4h_DMR_mal_heatmap.tiff")
plot_heatmap(DMR_mal_24h_percMeth, "./shortTerm_exp/plots/finalized_tiff/heatmaps/24h_DMR_mal_heatmap.tiff")
plot_heatmap(DMR_mal_72h_percMeth, "./shortTerm_exp/plots/finalized_tiff/heatmaps/72h_DMR_mal_heatmap.tiff")

#combined time points
plot_heatmap_combined(DMS_combined_fem_meth, "./shortTerm_exp/plots/finalized_tiff/heatmaps/DMS_combined_fem_heatmap.tiff")
plot_heatmap_combined(DMS_combined_mal_meth, "./shortTerm_exp/plots/finalized_tiff/heatmaps/DMS_combined_mal_heatmap.tiff")

plot_heatmap_combined(DMR_combined_fem_meth, "./shortTerm_exp/plots/finalized_tiff/heatmaps/DMR_combined_fem_heatmap.tiff")
plot_heatmap_combined(DMR_combined_mal_meth, "./shortTerm_exp/plots/finalized_tiff/heatmaps/DMR_combined_mal_heatmap.tiff")

