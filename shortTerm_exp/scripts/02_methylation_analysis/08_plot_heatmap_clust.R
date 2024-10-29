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

#run for all datasets
DMS_fem_all_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res_od/DMS_res_all/DMSdiffmeth_all_fem_5X.RDS",
                                            "./shortTerm_exp/data/methylkit_res_od/DMS_res_all/DMSmeth_all_fem_5X.RDS")

DMS_fem_05h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res_od/DMS_res_05h/DMSdiffmeth_05h_fem_5X.RDS",
                                           "./shortTerm_exp/data/methylkit_res_od/DMS_res_05h/DMSmeth_05h_fem_5X.RDS")

DMS_fem_1h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res_od/DMS_res_1h/DMSdiffmeth_1h_fem_5X.RDS",
                                           "./shortTerm_exp/data/methylkit_res_od/DMS_res_1h/DMSmeth_1h_fem_5X.RDS")

DMS_fem_4h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res_od/DMS_res_4h/DMSdiffmeth_4h_fem_5X.RDS",
                                           "./shortTerm_exp/data/methylkit_res_od/DMS_res_4h/DMSmeth_4h_fem_5X.RDS")

DMS_fem_24h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res_od/DMS_res_24h/DMSdiffmeth_24h_fem_5X.RDS",
                                           "./shortTerm_exp/data/methylkit_res_od/DMS_res_24h/DMSmeth_24h_fem_5X.RDS")

DMS_fem_72h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res_od/DMS_res_72h/DMSdiffmeth_72h_fem_5X.RDS",
                                           "./shortTerm_exp/data/methylkit_res_od/DMS_res_72h/DMSmeth_72h_fem_5X.RDS")

DMS_mal_all_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res_od/DMS_res_all/DMSdiffmeth_all_mal_5X.RDS",
                                            "./shortTerm_exp/data/methylkit_res_od/DMS_res_all/DMSmeth_all_mal_5X.RDS")

DMS_mal_05h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res_od/DMS_res_05h/DMSdiffmeth_05h_mal_5X.RDS",
                                            "./shortTerm_exp/data/methylkit_res_od/DMS_res_05h/DMSmeth_05h_mal_5X.RDS")

DMS_mal_1h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res_od/DMS_res_1h/DMSdiffmeth_1h_mal_5X.RDS",
                                           "./shortTerm_exp/data/methylkit_res_od/DMS_res_1h/DMSmeth_1h_mal_5X.RDS")

DMS_mal_4h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res_od/DMS_res_4h/DMSdiffmeth_4h_mal_5X.RDS",
                                           "./shortTerm_exp/data/methylkit_res_od/DMS_res_4h/DMSmeth_4h_mal_5X.RDS")

DMS_mal_24h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res_od/DMS_res_24h/DMSdiffmeth_24h_mal_5X.RDS",
                                            "./shortTerm_exp/data/methylkit_res_od/DMS_res_24h/DMSmeth_24h_mal_5X.RDS")

DMS_mal_72h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res_od/DMS_res_72h/DMSdiffmeth_72h_mal_5X.RDS",
                                            "./shortTerm_exp/data/methylkit_res_od/DMS_res_72h/DMSmeth_72h_mal_5X.RDS")

DMR_fem_all_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res_od/DMR_res_all/DMRdiffmeth_all_fem_5X.RDS",
                                            "./shortTerm_exp/data/methylkit_res_od/DMR_res_all/DMRmeth_all_fem_5X.RDS")

DMR_fem_05h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res_od/DMR_res_05h/DMRdiffmeth_05h_fem_5X.RDS",
                                            "./shortTerm_exp/data/methylkit_res_od/DMR_res_05h/DMR_tile_meth_05h_fem_5X.RDS")

DMR_fem_1h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res_od/DMR_res_1h/DMRdiffmeth_1h_fem_5X.RDS",
                                           "./shortTerm_exp/data/methylkit_res_od/DMR_res_1h/DMR_tile_meth_1h_fem_5X.RDS")

DMR_fem_4h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res_od/DMR_res_4h/DMRdiffmeth_4h_fem_5X.RDS",
                                           "./shortTerm_exp/data/methylkit_res_od/DMR_res_4h/DMR_tile_meth_4h_fem_5X.RDS")

DMR_fem_24h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res_od/DMR_res_24h/DMRdiffmeth_24h_fem_5X.RDS",
                                            "./shortTerm_exp/data/methylkit_res_od/DMR_res_24h/DMR_tile_meth_24h_fem_5X.RDS")

DMR_fem_72h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res_od/DMR_res_72h/DMRdiffmeth_72h_fem_5X.RDS",
                                            "./shortTerm_exp/data/methylkit_res_od/DMR_res_72h/DMR_tile_meth_72h_fem_5X.RDS")

DMR_mal_all_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res_od/DMR_res_all/DMRdiffmeth_all_mal_5X.RDS",
                                            "./shortTerm_exp/data/methylkit_res_od/DMR_res_all/DMR_tile_meth_all_mal_5X.RDS")

DMR_mal_05h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res_od/DMR_res_05h/DMRdiffmeth_05h_mal_5X.RDS",
                                            "./shortTerm_exp/data/methylkit_res_od/DMR_res_05h/DMR_tile_meth_05h_mal_5X.RDS")

DMR_mal_1h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res_od/DMR_res_1h/DMRdiffmeth_1h_mal_5X.RDS",
                                           "./shortTerm_exp/data/methylkit_res_od/DMR_res_1h/DMR_tile_meth_1h_mal_5X.RDS")

DMR_mal_4h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res_od/DMR_res_4h/DMRdiffmeth_4h_mal_5X.RDS",
                                           "./shortTerm_exp/data/methylkit_res_od/DMR_res_4h/DMR_tile_meth_4h_mal_5X.RDS")

DMR_mal_24h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res_od/DMR_res_24h/DMRdiffmeth_24h_mal_5X.RDS",
                                            "./shortTerm_exp/data/methylkit_res_od/DMR_res_24h/DMR_tile_meth_24h_mal_5X.RDS")

DMR_mal_72h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res_od/DMR_res_72h/DMRdiffmeth_72h_mal_5X.RDS",
                                            "./shortTerm_exp/data/methylkit_res_od/DMR_res_72h/DMR_tile_meth_72h_mal_5X.RDS")

#get methylation matrix for all sites and regions 
#read data files
meth_allsites <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMS_res_all/DMSmeth_all_5X.RDS")
meth_allsites_fem <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMS_res_all/DMSmeth_all_fem_5X.RDS")
meth_allsites_mal <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMS_res_all/DMSmeth_all_mal_5X.RDS")

meth_allregions <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMR_res_all/DMR_tile_meth_all_5X.RDS")
meth_allregions_fem <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMR_res_all/DMR_tile_meth_all_fem_5X.RDS")
meth_allregions_mal <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMR_res_all/DMR_tile_meth_all_mal_5X.RDS")

#calculate percent methylation 
perc_meth_allsites <- as.data.frame(percMethylation(meth_allsites))
perc_meth_allsites_fem <- as.data.frame(percMethylation(meth_allsites_fem))
perc_meth_allsites_mal <- as.data.frame(percMethylation(meth_allsites_mal))

perc_meth_allregions <- as.data.frame(percMethylation(meth_allregions))
perc_meth_allregions_fem <- as.data.frame(percMethylation(meth_allregions_fem))
perc_meth_allregions_mal <- as.data.frame(percMethylation(meth_allregions_mal))

#add site names
meth_allsites <- nameSite(meth_allsites)
meth_allsites_fem <- nameSite(meth_allsites_fem)
meth_allsites_mal <- nameSite(meth_allsites_mal)

meth_allregions <- nameSite(meth_allregions)
meth_allregions_fem <- nameSite(meth_allregions_fem)
meth_allregions_mal <- nameSite(meth_allregions_mal)

#add rownames to perc meth matrix  
rownames(perc_meth_allsites) <- meth_allsites$site_name
rownames(perc_meth_allsites_fem) <- meth_allsites_fem$site_name
rownames(perc_meth_allsites_mal) <- meth_allsites_mal$site_name

rownames(perc_meth_allregions) <- meth_allregions$site_name
rownames(perc_meth_allregions_fem) <- meth_allregions_fem$site_name
rownames(perc_meth_allregions_mal) <- meth_allregions_mal$site_name

#remove extra file 
rm(meth_allsites)
rm(meth_allsites_fem)
rm(meth_allsites_mal)

rm(meth_allregions)
rm(meth_allregions_fem)
rm(meth_allregions_mal)

##Run PCAs## 
#get lists of all DMSs for each sex
DMS_combined_all <- c(row.names(DMS_fem_05h_percMeth), row.names(DMS_fem_1h_percMeth), row.names(DMS_fem_4h_percMeth), 
                  row.names(DMS_fem_24h_percMeth), row.names(DMS_fem_72h_percMeth),
                  row.names(DMS_mal_05h_percMeth), row.names(DMS_mal_1h_percMeth), row.names(DMS_mal_4h_percMeth), 
                  row.names(DMS_mal_24h_percMeth), row.names(DMS_mal_72h_percMeth))

DMS_combined_all <- unique(DMS_combined_all) #remove duplicates

DMS_combined_fem <- c(row.names(DMS_fem_05h_percMeth), row.names(DMS_fem_1h_percMeth), row.names(DMS_fem_4h_percMeth), 
                      row.names(DMS_fem_24h_percMeth), row.names(DMS_fem_72h_percMeth))

DMS_combined_fem <- unique(DMS_combined_fem) #remove duplicates

DMS_combined_mal <- c(row.names(DMS_mal_05h_percMeth), row.names(DMS_mal_1h_percMeth), row.names(DMS_mal_4h_percMeth), 
                      row.names(DMS_mal_24h_percMeth), row.names(DMS_mal_72h_percMeth))

DMS_combined_mal <- unique(DMS_combined_mal) #remove duplicates

DMR_combined_all <- c(row.names(DMR_fem_05h_percMeth), row.names(DMR_fem_1h_percMeth), row.names(DMR_fem_4h_percMeth), 
                  row.names(DMR_fem_24h_percMeth), row.names(DMR_fem_72h_percMeth),
                  row.names(DMR_mal_05h_percMeth), row.names(DMR_mal_1h_percMeth), row.names(DMR_mal_4h_percMeth), 
                  row.names(DMR_mal_24h_percMeth), row.names(DMR_mal_72h_percMeth))

DMR_combined_all <- unique(DMR_combined_all) #remove duplicates

DMR_combined_fem <- c(row.names(DMR_fem_05h_percMeth), row.names(DMR_fem_1h_percMeth), row.names(DMR_fem_4h_percMeth), 
                  row.names(DMR_fem_24h_percMeth), row.names(DMR_fem_72h_percMeth))

DMR_combined_fem <- unique(DMR_combined_fem) #remove duplicates

DMR_combined_mal <- c(row.names(DMR_mal_05h_percMeth), row.names(DMR_mal_1h_percMeth), row.names(DMR_mal_4h_percMeth), 
                  row.names(DMR_mal_24h_percMeth), row.names(DMR_mal_72h_percMeth))

DMR_combined_mal <- unique(DMR_combined_mal) #remove duplicates

#extract DMSs and DMRs
DMS_combined_all_meth <- perc_meth_allsites %>% filter(row.names(perc_meth_allsites) %in% DMS_combined_all)
DMS_combined_fem_meth <- perc_meth_allsites_fem %>% filter(row.names(perc_meth_allsites_fem) %in% DMS_combined_fem)
DMS_combined_mal_meth <- perc_meth_allsites_mal %>% filter(row.names(perc_meth_allsites_mal) %in% DMS_combined_mal)

DMR_combined_all_meth <- perc_meth_allregions %>% filter(row.names(perc_meth_allregions) %in% DMR_combined_all)
DMR_combined_fem_meth <- perc_meth_allregions_fem %>% filter(row.names(perc_meth_allregions_fem) %in% DMR_combined_fem)
DMR_combined_mal_meth <- perc_meth_allregions_mal %>% filter(row.names(perc_meth_allregions_mal) %in% DMR_combined_mal)

#make function to transpose data 
prep_data <- function(df){
  #transpose data
  t.data <- t(df)
  rownames(t.data) <- colnames(df)
  colnames(t.data) <- rownames(df)
  return(as.data.frame(t.data))
}

#for all DMSs/DMRs for all sexes together and separate
t_DMS_combined_meth_all <- prep_data(DMS_combined_all_meth)
t_DMR_combined_meth_all <- prep_data(DMR_combined_all_meth)

t_DMS_combined_meth_fem <- prep_data(DMS_combined_fem_meth)
t_DMR_combined_meth_fem <- prep_data(DMR_combined_fem_meth)

t_DMS_combined_meth_mal <- prep_data(DMS_combined_mal_meth)
t_DMR_combined_meth_mal <- prep_data(DMR_combined_mal_meth)

#for all CpGs
t_DMS_fem_all <- prep_data(DMS_fem_all_percMeth)
t_DMS_mal_all <- prep_data(DMS_mal_all_percMeth)

#run PCA
pca_DMS_combined_all <- prcomp(t_DMS_combined_meth_all, center = TRUE)
pca_DMR_combined_all <- prcomp(t_DMR_combined_meth_all, center = TRUE)

pca_DMS_combined_fem <- prcomp(t_DMS_combined_meth_fem, center = TRUE)
pca_DMR_combined_fem <- prcomp(t_DMR_combined_meth_fem, center = TRUE)

pca_DMS_combined_mal <- prcomp(t_DMS_combined_meth_mal, center = TRUE)
pca_DMR_combined_mal <- prcomp(t_DMR_combined_meth_mal, center = TRUE)

pca_DMS_fem_all <- prcomp(t_DMS_fem_all, center = TRUE)
pca_DMS_mal_all <- prcomp(t_DMS_mal_all, center = TRUE)

#view results 
summary(pca_DMS_combined_all)
summary(pca_DMR_combined_all)
summary(pca_DMS_combined_fem)
summary(pca_DMR_combined_fem)
summary(pca_DMS_combined_mal)
summary(pca_DMR_combined_mal)
summary(pca_DMS_fem_all)
summary(pca_DMS_mal_all)

#make lists of IDs for each time point 
IDs_05h <- c(colnames(DMS_fem_05h_percMeth), colnames(DMS_mal_05h_percMeth))
IDs_1h <- c(colnames(DMS_fem_1h_percMeth), colnames(DMS_mal_1h_percMeth))
IDs_4h <- c(colnames(DMS_fem_4h_percMeth), colnames(DMS_mal_4h_percMeth))
IDs_24h <- c(colnames(DMS_fem_24h_percMeth), colnames(DMS_mal_24h_percMeth))
IDs_72h <- c(colnames(DMS_fem_72h_percMeth), colnames(DMS_mal_72h_percMeth))

#add on tank, sex, and cue 
add_var <- function(data, pca) {
  #add first four PCs to df
  data <- cbind(data, pca$x[,1:4])
  
  #convert rownames to ID column 
  data.new <- tibble::rownames_to_column(data, "ID")
  
  #add sex
  data.new$sex <- data.new$ID
  data.new$sex <- str_sub(data.new$sex, -1)
  
  #add cue
  data.new$cue <- data.new$ID
  data.new$cue <- str_sub(data.new$cue, 4, 4)
  
  #add tank
  data.new$tank <- data.new$ID
  data.new$tank <- str_sub(data.new$tank, 4, 7)
  data.new$tank <- gsub("F", "", data.new$tank)
  data.new$tank <- gsub("M", "", data.new$tank)
  
  #add time point 
  data.new$time_point <- "NA"
  data.new <- within(data.new, time_point[ID %in% IDs_05h] <- "05h")
  data.new <- within(data.new, time_point[ID %in% IDs_1h] <- "1h")
  data.new <- within(data.new, time_point[ID %in% IDs_4h] <- "4h")
  data.new <- within(data.new, time_point[ID %in% IDs_24h] <- "24h")
  data.new <- within(data.new, time_point[ID %in% IDs_72h] <- "72h")
  
  #put rownames back on 
  rownames(data.new) <- rownames(data)
  return(data.new)
}

t_DMS_combined_meth_all <- add_var(t_DMS_combined_meth_all, pca_DMS_combined_all)
t_DMR_combined_meth_all <- add_var(t_DMR_combined_meth_all, pca_DMR_combined_all)
t_DMS_combined_meth_fem <- add_var(t_DMS_combined_meth_fem, pca_DMS_combined_fem)
t_DMR_combined_meth_fem <- add_var(t_DMR_combined_meth_fem, pca_DMR_combined_fem)
t_DMS_combined_meth_mal <- add_var(t_DMS_combined_meth_mal, pca_DMS_combined_mal)
t_DMR_combined_meth_mal <- add_var(t_DMR_combined_meth_mal, pca_DMR_combined_mal)
t_DMS_fem_all <- add_var(t_DMS_fem_all, pca_DMS_fem_all)
t_DMS_mal_all <- add_var(t_DMS_mal_all, pca_DMS_mal_all)

#plot PCAs
plot_pca <- function(pca.data, pca.x, pca.y, xlabel, ylabel, plotname){
  p <- ggplot(data = pca.data, aes(x = pca.x, y = pca.y, color = cue, shape = time_point)) + theme_bw() + geom_point() +
        scale_color_manual(name = "Cue",
                           values = c("#E14D2A", "#31C6D4")) +
        geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
        xlab(xlabel) + ylab(ylabel) +
        theme(axis.title = element_text(size = 18), axis.text = element_text(size = 16), legend.text = element_text(size=18),
              legend.title = element_text(size = 18))
  
  ggsave(filename = plotname, plot = p, width = 6, height = 5, units = "in", dpi = 300)
}

plot_pca(t_DMS_combined_meth_all, t_DMS_combined_meth_all$PC1, t_DMS_combined_meth_all$PC2, "PC1 (20.25%)", "PC2 (10.94%)",
         "./shortTerm_exp/plots/finalized_tiff/pca_plots_od/DMS_combined_all_pca.tiff")
plot_pca(t_DMR_combined_meth_all, t_DMR_combined_meth_all$PC1, t_DMR_combined_meth_all$PC2, "PC1 (9.11%)", "PC2 (5.02%)",
         "./shortTerm_exp/plots/finalized_tiff/pca_plots_od/DMR_combined_all_pca.tiff")

plot_pca(t_DMS_combined_meth_fem, t_DMS_combined_meth_fem$PC1, t_DMS_combined_meth_fem$PC2, "PC1 (29.62%)", "PC2 (12.17%)",
         "./shortTerm_exp/plots/finalized_tiff/pca_plots_od/DMS_combined_fem_pca.tiff")
plot_pca(t_DMR_combined_meth_fem, t_DMR_combined_meth_fem$PC1, t_DMR_combined_meth_fem$PC2, "PC1 (13.21%)", "PC2 (8.49%)",
         "./shortTerm_exp/plots/finalized_tiff/pca_plots_od/DMR_combined_fem_pca.tiff")

plot_pca(t_DMS_combined_meth_mal, t_DMS_combined_meth_mal$PC1, t_DMS_combined_meth_mal$PC2, "PC1 (20.48%)", "PC2 (11.42%)",
         "./shortTerm_exp/plots/finalized_tiff/pca_plots_od/DMS_combined_mal_pca.tiff")
plot_pca(t_DMR_combined_meth_mal, t_DMR_combined_meth_mal$PC1, t_DMR_combined_meth_mal$PC2, "PC1 (10.68%)", "PC2 (7.04%)",
         "./shortTerm_exp/plots/finalized_tiff/pca_plots_od/DMR_combined_mal_pca.tiff")

# plot_pca(t_DMS_combined_meth_all, t_DMS_combined_meth_all$PC1, t_DMS_combined_meth_all$PC2, "PC1 (26.31%)", "PC2 (13.46%)",
#          "./shortTerm_exp/plots/finalized_tiff/pca_plots/DMS_combined_all_pca.tiff")
# plot_pca(t_DMR_combined_meth_all, t_DMR_combined_meth_all$PC1, t_DMR_combined_meth_all$PC2, "PC1 (12.78%)", "PC2 (6.58%)",
#          "./shortTerm_exp/plots/finalized_tiff/pca_plots/DMR_combined_all_pca.tiff")
# 
# plot_pca(t_DMS_combined_meth_fem, t_DMS_combined_meth_fem$PC1, t_DMS_combined_meth_fem$PC2, "PC1 (33.84%)", "PC2 (14.16%)",
#          "./shortTerm_exp/plots/finalized_tiff/pca_plots/DMS_combined_fem_pca.tiff")
# plot_pca(t_DMR_combined_meth_fem, t_DMR_combined_meth_fem$PC1, t_DMR_combined_meth_fem$PC2, "PC1 (16.55%)", "PC2 (8.37%)",
#          "./shortTerm_exp/plots/finalized_tiff/pca_plots/DMR_combined_fem_pca.tiff")
# 
# plot_pca(t_DMS_combined_meth_mal, t_DMS_combined_meth_mal$PC1, t_DMS_combined_meth_mal$PC2, "PC1 (24.88%)", "PC2 (16.09%)",
#          "./shortTerm_exp/plots/finalized_tiff/pca_plots/DMS_combined_mal_pca.tiff")
# plot_pca(t_DMR_combined_meth_mal, t_DMR_combined_meth_mal$PC1, t_DMR_combined_meth_mal$PC2, "PC1 (13.74%)", "PC2 (8.67%)",
#          "./shortTerm_exp/plots/finalized_tiff/pca_plots/DMR_combined_mal_pca.tiff")

# plot_pca(t_DMS_fem_all, t_DMS_fem_all$PC1, t_DMS_fem_all$PC2, "PC1 (30.9%)", "PC2 (13.08%)",
#          "./shortTerm_exp/plots/finalized_tiff/pca_plots/DMS_all_fem_pca.tiff")
# plot_pca(t_DMS_mal_all, t_DMS_mal_all$PC1, t_DMS_mal_all$PC2, "PC1 (36.9%)", "PC2 (11.5%)",
#          "./shortTerm_exp/plots/finalized_tiff/pca_plots/DMS_all_mal_pca.tiff")

##Make heat maps ##
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

#plot for each time point/sex
plot_heatmap(DMS_fem_all_percMeth, "./shortTerm_exp/plots/finalized_tiff/heatmaps_od/allTP_DMS_fem_heatmap.tiff")
plot_heatmap(DMS_fem_05h_percMeth, "./shortTerm_exp/plots/finalized_tiff/heatmaps_od/05h_DMS_fem_heatmap.tiff")
plot_heatmap(DMS_fem_1h_percMeth, "./shortTerm_exp/plots/finalized_tiff/heatmaps_od/1h_DMS_fem_heatmap.tiff")
plot_heatmap(DMS_fem_4h_percMeth, "./shortTerm_exp/plots/finalized_tiff/heatmaps_od/4h_DMS_fem_heatmap.tiff")
plot_heatmap(DMS_fem_24h_percMeth, "./shortTerm_exp/plots/finalized_tiff/heatmaps_od/24h_DMS_fem_heatmap.tiff")
plot_heatmap(DMS_fem_72h_percMeth, "./shortTerm_exp/plots/finalized_tiff/heatmaps_od/72h_DMS_fem_heatmap.tiff")

plot_heatmap(DMS_mal_all_percMeth, "./shortTerm_exp/plots/finalized_tiff/heatmaps_od/allTP_DMS_mal_heatmap.tiff")
plot_heatmap(DMS_mal_05h_percMeth, "./shortTerm_exp/plots/finalized_tiff/heatmaps_od/05h_DMS_mal_heatmap.tiff")
plot_heatmap(DMS_mal_1h_percMeth, "./shortTerm_exp/plots/finalized_tiff/heatmaps_od/1h_DMS_mal_heatmap.tiff")
plot_heatmap(DMS_mal_4h_percMeth, "./shortTerm_exp/plots/finalized_tiff/heatmaps_od/4h_DMS_mal_heatmap.tiff")
plot_heatmap(DMS_mal_24h_percMeth, "./shortTerm_exp/plots/finalized_tiff/heatmaps_od/24h_DMS_mal_heatmap.tiff")
plot_heatmap(DMS_mal_72h_percMeth, "./shortTerm_exp/plots/finalized_tiff/heatmaps_od/72h_DMS_mal_heatmap.tiff")

#DMRs
plot_heatmap(DMR_fem_05h_percMeth, "./shortTerm_exp/plots/finalized_tiff/heatmaps_od/05h_DMR_fem_heatmap.tiff")
plot_heatmap(DMR_fem_1h_percMeth, "./shortTerm_exp/plots/finalized_tiff/heatmaps_od/1h_DMR_fem_heatmap.tiff")
plot_heatmap(DMR_fem_4h_percMeth, "./shortTerm_exp/plots/finalized_tiff/heatmaps_od/4h_DMR_fem_heatmap.tiff")
plot_heatmap(DMR_fem_24h_percMeth, "./shortTerm_exp/plots/finalized_tiff/heatmaps_od/24h_DMR_fem_heatmap.tiff")
plot_heatmap(DMR_fem_72h_percMeth, "./shortTerm_exp/plots/finalized_tiff/heatmaps_od/72h_DMR_fem_heatmap.tiff")

plot_heatmap(DMR_mal_05h_percMeth, "./shortTerm_exp/plots/finalized_tiff/heatmaps_od/05h_DMR_mal_heatmap.tiff")
plot_heatmap(DMR_mal_1h_percMeth, "./shortTerm_exp/plots/finalized_tiff/heatmaps_od/1h_DMR_mal_heatmap.tiff")
plot_heatmap(DMR_mal_4h_percMeth, "./shortTerm_exp/plots/finalized_tiff/heatmaps_od/4h_DMR_mal_heatmap.tiff")
plot_heatmap(DMR_mal_24h_percMeth, "./shortTerm_exp/plots/finalized_tiff/heatmaps_od/24h_DMR_mal_heatmap.tiff")
plot_heatmap(DMR_mal_72h_percMeth, "./shortTerm_exp/plots/finalized_tiff/heatmaps_od/72h_DMR_mal_heatmap.tiff")

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

plot_heatmap_combined(DMS_combined_fem_meth, "./shortTerm_exp/plots/finalized_tiff/heatmaps_od/DMS_combined_fem_heatmap.tiff")
plot_heatmap_combined(DMS_combined_mal_meth, "./shortTerm_exp/plots/finalized_tiff/heatmaps_od/DMS_combined_mal_heatmap.tiff")

plot_heatmap_combined(DMR_combined_fem_meth, "./shortTerm_exp/plots/finalized_tiff/heatmaps_od/DMR_combined_fem_heatmap.tiff")
plot_heatmap_combined(DMR_combined_mal_meth, "./shortTerm_exp/plots/finalized_tiff/heatmaps_od/DMR_combined_mal_heatmap.tiff")

