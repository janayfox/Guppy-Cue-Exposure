DMS_mydiff_all_fem <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMS_res_all/DMSmydiff_all_fem_5X_data.RDS")
DMS_mydiff_05h_fem <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMS_res_05h/DMSmydiff_05h_fem_5X_data.RDS")
DMS_mydiff_1h_fem <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMS_res_1h/DMSmydiff_1h_fem_5X_data.RDS")
DMS_mydiff_4h_fem <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMS_res_4h/DMSmydiff_4h_fem_5X_data.RDS")
DMS_mydiff_24h_fem <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMS_res_24h/DMSmydiff_24h_fem_5X_data.RDS")
DMS_mydiff_72h_fem <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMS_res_72h/DMSmydiff_72h_fem_5X_data.RDS")

DMS_mydiff_all_mal <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMS_res_all/DMSmydiff_all_mal_5X_data.RDS")
DMS_mydiff_05h_mal <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMS_res_05h/DMSmydiff_05h_mal_5X_data.RDS")
DMS_mydiff_1h_mal <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMS_res_1h/DMSmydiff_1h_mal_5X_data.RDS")
DMS_mydiff_4h_mal <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMS_res_4h/DMSmydiff_4h_mal_5X_data.RDS")
DMS_mydiff_24h_mal <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMS_res_24h/DMSmydiff_24h_mal_5X_data.RDS")
DMS_mydiff_72h_mal <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMS_res_72h/DMSmydiff_72h_mal_5X_data.RDS")

DMR_mydiff_all_fem <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMR_res_all/DMRmydiff_all_fem_5X_data.RDS")
DMR_mydiff_05h_fem <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMR_res_05h/DMRmydiff_05h_fem_5X_data.RDS")
DMR_mydiff_1h_fem <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMR_res_1h/DMRmydiff_1h_fem_5X_data.RDS")
DMR_mydiff_4h_fem <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMR_res_4h/DMRmydiff_4h_fem_5X_data.RDS")
DMR_mydiff_24h_fem <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMR_res_24h/DMRmydiff_24h_fem_5X_data.RDS")
DMR_mydiff_72h_fem <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMR_res_72h/DMRmydiff_72h_fem_5X_data.RDS")

DMR_mydiff_all_mal <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMR_res_all/DMRmydiff_all_mal_5X_data.RDS")
DMR_mydiff_05h_mal <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMR_res_05h/DMRmydiff_05h_mal_5X_data.RDS")
DMR_mydiff_1h_mal <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMR_res_1h/DMRmydiff_1h_mal_5X_data.RDS")
DMR_mydiff_4h_mal <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMR_res_4h/DMRmydiff_4h_mal_5X_data.RDS")
DMR_mydiff_24h_mal <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMR_res_24h/DMRmydiff_24h_mal_5X_data.RDS")
DMR_mydiff_72h_mal <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMR_res_72h/DMRmydiff_72h_mal_5X_data.RDS")

library(methylKit)

#call significant methylation
DMS_1h_fem <- DMS_mydiff_1h_fem[DMS_mydiff_1h_fem$qvalue <= 0.0125 & abs(DMS_mydiff_1h_fem$meth.diff) > 15 ,]
DMS_1h_mal <- DMS_mydiff_1h_mal[DMS_mydiff_1h_mal$qvalue <= 0.0125 & abs(DMS_mydiff_1h_mal$meth.diff) > 15 ,]


#get prec. methylation for each feature 
#read in anno results and extract
cpg_anno_fem <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_all/CpG_anno_all_fem.RDS")
cpg_anno_fem_res <- as.data.frame(cpg_anno_fem@members)

cpg_anno_mal <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_all/CpG_anno_all_mal.RDS")
cpg_anno_mal_res <- as.data.frame(cpg_anno_mal@members)

#set label 
cpg_anno_fem_res$feature <- ifelse(cpg_anno_fem_res$prom == 1, "promoter",
                                   ifelse(cpg_anno_fem_res$exon == 1, "exon",
                                          ifelse(cpg_anno_fem_res$intron == 1, "intron", "intergenic")))

cpg_anno_mal_res$feature <- ifelse(cpg_anno_mal_res$prom == 1, "promoter",
                                   ifelse(cpg_anno_mal_res$exon == 1, "exon",
                                          ifelse(cpg_anno_mal_res$intron == 1, "intron", "intergenic")))

#get percent methylation
library(methylKit)
meth_fem <- readRDS("./shortTerm_exp/data/methylkit_res/DMS_res_all/DMSmeth_all_fem_5X.RDS")
meth_mal <- readRDS("./shortTerm_exp/data/methylkit_res/DMS_res_all/DMSmeth_all_mal_5X.RDS")

percmeth_fem <- as.data.frame(percMethylation(meth_fem))
percmeth_mal <- as.data.frame(percMethylation(meth_mal))

#add feature col on 
percmeth_fem$feature <- cpg_anno_fem_res$feature
percmeth_mal$feature <- cpg_anno_mal_res$feature

library(tidyr)
library(ggplot2)

# Reshape the data from wide to long format
fem_long <- percmeth_fem %>%
  pivot_longer(cols = -feature, names_to = "individual", values_to = "percent_methylation")

mal_long <- percmeth_mal %>%
  pivot_longer(cols = -feature, names_to = "individual", values_to = "percent_methylation")

#calculate averages
fem_long %>%
  group_by(feature) %>%
  summarise(
    average_percent_methylation = mean(percent_methylation, na.rm = TRUE),
    sd_percent_methylation = sd(percent_methylation, na.rm = TRUE)
  )

mal_long %>%
  group_by(feature) %>%
  summarise(
    average_percent_methylation = mean(percent_methylation, na.rm = TRUE),
    sd_percent_methylation = sd(percent_methylation, na.rm = TRUE)
  )
#make plotss 

# Sample 1000 rows per feature
fem_sample <- fem_long %>%
  group_by(feature) %>%
  slice_sample(n = 1000) %>%
  ungroup()

mal_sample <- mal_long %>%
  group_by(feature) %>%
  slice_sample(n = 1000) %>%
  ungroup()

#make boxplot 
fem.plot <- ggplot(fem_sample, aes(x = feature, y = percent_methylation)) +
                      geom_boxplot(fill = "lightblue") +
                      labs(x = "Feature", y = "Percent Methylation") + 
                      theme_bw() 
fem.plot

mal.plot <- ggplot(mal_sample, aes(x = feature, y = percent_methylation)) +
  geom_boxplot(fill = "lightblue") +
  labs(x = "Feature", y = "Percent Methylation") + 
  theme_bw() 
mal.plot

panel <- ggarrange(fem.plot,mal.plot, labels = c("A", "B"))
panel

ggsave("percmeth_boxplot.tiff", plot = panel, device = "tiff", width = 6, height = 4, units = "in")
