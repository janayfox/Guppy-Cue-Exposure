#######################################################
### Goal: Cluster samples based on CpGs, DMS, and DMR
### Author: Janay Fox
### R script
#######################################################

## Set up ## 
#load packages
library(ggplot2)
library(factoextra)
library(tidyverse)
library(vegan)
library(goeveg)
library(ggfortify)
library(ggforce)
library(cluster)
library(dendextend)
library(RColorBrewer)
library(tidyr)
library(dplyr)

set.seed(10)

#read in data
methAllData21L <- readRDS("./gup_cue_exp/dev_exp/data/methylkit_res/methAllData21L.RDS")
diffMeth21L <- readRDS("./gup_cue_exp/dev_exp/data/methylkit_res/DMSmethAllData21L.RDS")

#add in month of assay data for each tank
methAllData21L$monthAssay <- NA
methAllData21L$monthAssay <- ifelse(methAllData21L$tank %in% c("DAC2", "DC2", "DAC3", "DC3"), "May", methAllData21L$monthAssay)
methAllData21L$monthAssay <- ifelse(methAllData21L$tank %in% c("DAC4", "DC4", "DAC5", "DC5"), "August", methAllData21L$monthAssay)
methAllData21L$monthAssay <- ifelse(methAllData21L$tank %in% c("DAC6", "DC6", "DAC7", "DC7"), "September", methAllData21L$monthAssay)

diffMeth21L$monthAssay <- NA
diffMeth21L$monthAssay <- ifelse(diffMeth21L$tank %in% c("DAC2", "DC2", "DAC3", "DC3"), "May", diffMeth21L$monthAssay)
diffMeth21L$monthAssay <- ifelse(diffMeth21L$tank %in% c("DAC4", "DC4", "DAC5", "DC5"), "August", diffMeth21L$monthAssay)
diffMeth21L$monthAssay <- ifelse(diffMeth21L$tank %in% c("DAC6", "DC6", "DAC7", "DC7"), "September", diffMeth21L$monthAssay)

#calculate which quartile fish are in for shoaling behaviour and open field 
methAllData21L <- mutate(methAllData21L, quantile_tight_diff = ntile(methAllData21L$tight_diff_s, 4))
methAllData21L <- mutate(methAllData21L, quantile_PC1 = ntile(methAllData21L$PC1, 4))
methAllData21L <- mutate(methAllData21L, quantile_dist = ntile(methAllData21L$dist_cm, 4))

diffMeth21L <- mutate(diffMeth21L, quantile_tight_diff = ntile(diffMeth21L$tight_diff_s, 4))
diffMeth21L <- mutate(diffMeth21L, quantile_PC1 = ntile(diffMeth21L$PC1, 4))
diffMeth21L <- mutate(diffMeth21L, quantile_dist = ntile(diffMeth21L$dist_cm, 4))

#convert all character columns to factors 
methAllData21L[sapply(methAllData21L, is.character)] <- lapply(methAllData21L[sapply(methAllData21L, is.character)], as.factor)

diffMeth21L[sapply(diffMeth21L, is.character)] <- lapply(diffMeth21L[sapply(diffMeth21L, is.character)], as.factor)

## CpG clustering ##
#remove uninformative CpGs (sd =< 10%)
methCpG <- methAllData21L[,2:29593] #make dataset of only CpGs
methCpG <- methCpG[,sapply(methCpG, function(x) {sd(x) >= 20})] #remove uninformative CpGs
keep.CpG <- colnames(methCpG) #make list of column names to keep
methAllData21L.CpG <- select(methAllData21L, keep.CpG) #new dataset of perc values of retained CpGs
methAllData21L <- methAllData21L[,-2:-29593] #remove CpGs
methAllData21L <- cbind(methAllData21L, methAllData21L.CpG) #bind with retained CpGs

#run PCA on all CpGs
methPCA21L <- prcomp(methAllData21L[,40:1263], scale. = TRUE, center = TRUE)

#check for impact of sex on first 6 PCs
sex.pca.plot.PC1.2 <- fviz_pca_ind(methPCA21L,
                             col.ind = methAllData21L$sex,
                             palette = c("#00AFBB",  "#FC4E07"),
                             addEllipses = TRUE, 
                             ellipse.type = "confidence",
                             legend.title = "Sex",
                             repel = TRUE     # Avoid text overlapping
)
sex.pca.plot.PC1.2

sex.pca.plot.PC3.4 <- fviz_pca_ind(methPCA21L,
                                   col.ind = methAllData21L$sex,
                                   palette = c("#00AFBB",  "#FC4E07"),
                                   addEllipses = TRUE, 
                                   ellipse.type = "confidence",
                                   legend.title = "Sex",
                                   repel = TRUE, # Avoid text overlapping
                                   axes = c(3,4)
)
sex.pca.plot.PC3.4

sex.pca.plot.PC5.6 <- fviz_pca_ind(methPCA21L,
                                   col.ind = methAllData21L$sex,
                                   palette = c("#00AFBB",  "#FC4E07"),
                                   addEllipses = TRUE, 
                                   ellipse.type = "confidence",
                                   legend.title = "Sex",
                                   repel = TRUE, # Avoid text overlapping
                                   axes = c(5,6)
)
sex.pca.plot.PC5.6 

#check for impact of tank 
tank.pca.plot.PC1.2 <- fviz_pca_ind(methPCA21L,
                              col.ind = methAllData21L$tank, 
                              addEllipses = TRUE, 
                              ellipse.type = "confidence",
                              legend.title = "Tank",
                              repel = TRUE     # Avoid text overlapping
)
tank.pca.plot.PC1.2

#extract first 2 PCs
methAllData21L <- cbind(methAllData21L, methPCA21L$x[,1:2])
colnames(methAllData21L)[1264] = "PC1_cpg"
colnames(methAllData21L)[1265] = "PC2_cpg"

#plot cpg PCA 
cpg.pca.plot <- ggplot(data = methAllData21L, aes(x = PC1_cpg, y = PC2_cpg, color = tank, shape = cue)) + theme_bw() + geom_point() +
                            scale_color_manual(values = c("#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#084594",
                                               "#FDD0A2", "#FDAE6B", "#FD8D3C", "#F16913", "#D94801", "#8C2D04")) + labs(color = "Tank", shape = "Cue", lty = "Cue") +
                            stat_ellipse(aes(group = cue, lty = cue)) +
                            xlab("PC1 (23.5%)") + ylab("PC2 (5.7%)")

cpg.pca.plot

## so plotting this is not easy, should I do tank by color or shape? too many tanks for shapes but maybe could just do diff shapes within a cue and then reuse between cues
## and put the cues as different colors? should I do ellipses around tanks or cue? If I do tank then DAC2 only has 3 individuals and so it causes an issue for some reason
## also its just not very easy to interpret if I do tanks anyways 

## RUN NMDS 
# create matrix of just perc meth values 
cpg.nmds.dat <- as.matrix(methAllData21L[,40:1263])

#run NMDS with Bray-Curtis 
cpg.nmds = metaMDS(cpg.nmds.dat, distance = "bray", autotransform = FALSE, trymax = 100, maxit = 500)

#extract NMDS scores
cpg.nmds.scores = as.data.frame(scores(cpg.nmds)$sites)

#check dimensions for projection of NMDS 
cpg.dim.plot <- dimcheckMDS(cpg.nmds.dat) #three dimensions with stress < 0.1 (for when I run on only informative CpGs)
cpg.goodness <- goodness(cpg.nmds)
cpg.stress.plot <- stressplot(cpg.nmds)

#add other data to nmds datafram 
cpg.nmds.scores$fish_ID <- methAllData21L$fish_ID
cpg.nmds.scores$cue <- methAllData21L$cue
cpg.nmds.scores$sex <- methAllData21L$sex
cpg.nmds.scores$tank <- methAllData21L$tank
cpg.nmds.scores$monthAssay <- methAllData21L$monthAssay


cpg.nmds.plot <- ggplot(cpg.nmds.scores, aes(x = NMDS1, y = NMDS2, color = tank, shape = cue)) + geom_point() + theme_bw() +
                        scale_color_manual(values = c("#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#084594",
                                                      "#FDD0A2", "#FDAE6B", "#FD8D3C", "#F16913", "#D94801", "#8C2D04")) + 
                        labs(color = "Tank", shape = "Cue", lty = "Cue") + stat_ellipse(aes(group = cue, lty = cue))
cpg.nmds.plot

## Determine cluster groups ##
#use average silhoutte method to classify individuals into clusters 
#create dataset of just perc meth 
percMeth <- methAllData21L[,40:1263]

#function to compute average silhoutte for k clusters
avg_sil <- function(k) {
  km.res <- kmeans(percMeth, centers = k, nstart = 25)
  ss <- silhouette(km.res$cluster, dist(percMeth))
  mean(ss[, 3])
}

#compute and plot for k = 2 to k = 12 (number of tanks)
k.values <- 2:12

#extract avg silouette for 2-12 clusters 
avg_sil_values <- map_dbl(k.values, avg_sil)

#plot results
plot(k.values, avg_sil_values, type = "b", pch = 19, frame = FALSE, xlab = "Number of clusters K", ylab = "Average Silhouettes")

#compute k-means clustering with k = 2
final.cluster <- kmeans(percMeth, 2, nstart = 25, iter.max = 100)
methAllData21L$cpgCluster <- as.factor(final.cluster$cluster) #add onto data
cpg.nmds.scores$cpgCluster <- as.factor(final.cluster$cluster)

## replot PCA and NMDS ##
cpg.nmds.plot.cluster <- ggplot(cpg.nmds.scores, aes(x = NMDS1, y = NMDS2, color = tank, shape = cue)) + geom_point() + theme_bw() +
  scale_color_manual(values = c("#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#084594",
                                "#FDD0A2", "#FDAE6B", "#FD8D3C", "#F16913", "#D94801", "#8C2D04")) + 
  labs(color = "Tank", shape = "Cue", lty = "Cluster") + stat_ellipse(aes(group = cpgCluster, lty = cpgCluster))
cpg.nmds.plot.cluster

cpg.pca.plot.cluster <- ggplot(data = methAllData21L, aes(x = PC1_cpg, y = PC2_cpg, color = tank, shape = cue)) + theme_bw() + geom_point() +
  scale_color_manual(values = c("#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#084594",
                                "#FDD0A2", "#FDAE6B", "#FD8D3C", "#F16913", "#D94801", "#8C2D04")) + labs(color = "Tank", shape = "Cue", lty = "Cluster") +
  stat_ellipse(aes(group = cpgCluster, lty = cpgCluster)) + xlab("PC1 (23.5%)") + ylab("PC2 (5.7%)")

cpg.pca.plot.cluster

#check variables against the clustering 
cpg.pca.plot.cluster.sex <- ggplot(data = methAllData21L, aes(x = PC1_cpg, y = PC2_cpg, color = sex, shape = cue)) + theme_bw() + geom_point() +
  labs(color = "Sex", shape = "Cue", lty = "Cluster") +
  stat_ellipse(aes(group = cpgCluster, lty = cpgCluster)) + xlab("PC1") + ylab("PC2")

cpg.pca.plot.cluster.sex

cpg.pca.plot.cluster.month <- ggplot(data = methAllData21L, aes(x = PC1_cpg, y = PC2_cpg, color = monthAssay, shape = cue)) + theme_bw() + geom_point() +
  labs(color = "Month Assay", shape = "Cue", lty = "Cluster") +
  stat_ellipse(aes(group = cpgCluster, lty = cpgCluster)) + xlab("PC1") + ylab("PC2")

cpg.pca.plot.cluster.month

cpg.pca.plot.cluster.time <- ggplot(data = methAllData21L, aes(x = PC1_cpg, y = PC2_cpg, color = time, shape = cue)) + theme_bw() + geom_point() +
  labs(color = "Time of Sampling", shape = "Cue", lty = "Cluster") +
  stat_ellipse(aes(group = cpgCluster, lty = cpgCluster)) + xlab("PC1") + ylab("PC2")

cpg.pca.plot.cluster.time

cpg.pca.plot.cluster.shoal <- ggplot(data = methAllData21L, aes(x = PC1_cpg, y = PC2_cpg, color = quantile_tight_diff, shape = cue)) + theme_bw() + geom_point() +
  labs(color = "Quantile tight shoal", shape = "Cue", lty = "Cluster") +
  stat_ellipse(aes(group = cpgCluster, lty = cpgCluster)) + xlab("PC1") + ylab("PC2")

cpg.pca.plot.cluster.shoal

cpg.pca.plot.cluster.PC1 <- ggplot(data = methAllData21L, aes(x = PC1_cpg, y = PC2_cpg, color = quantile_PC1, shape = cue)) + theme_bw() + geom_point() +
  labs(color = "Quantile PC1 Openfield", shape = "Cue", lty = "Cluster") +
  stat_ellipse(aes(group = cpgCluster, lty = cpgCluster)) + xlab("PC1") + ylab("PC2")

cpg.pca.plot.cluster.PC1

cpg.pca.plot.cluster.dist <- ggplot(data = methAllData21L, aes(x = PC1_cpg, y = PC2_cpg, color = quantile_dist, shape = cue)) + theme_bw() + geom_point() +
  labs(color = "Quantile distance travelled", shape = "Cue", lty = "Cluster") +
  stat_ellipse(aes(group = cpgCluster, lty = cpgCluster)) + xlab("PC1") + ylab("PC2")

cpg.pca.plot.cluster.dist

## Run hierarchical clustering ## 
#create dataset of just perc meth but with sample ID as row name 
hclust.dat <- methAllData21L[,c(1,40:1263)]
rownames(hclust.dat) <- hclust.dat[,1]
hclust.dat[,1] <- NULL

#calculate a Euclidean distance matrix 
hclust.dist <- vegdist(hclust.dat, method = "euclidean")

#fit cluster dendrogram
cpg.hclust <- hclust(hclust.dist, method="ward.D2")
cpg.den <- as.dendrogram(cpg.hclust)

#plot dendrogram
#get color pallette values 
brewer.pal(12, "Set3")

#make bars for color bars 
sex_bar <- ifelse(methAllData21L$sex == "female", "#c51b7d", "blue")
tanks <- methAllData21L[,"tank"]
cols_12 <- RColorBrewer::brewer.pal(12, "Paired")
tank_bar <- cols_12[tanks]

#bind bars together
bars <- cbind(sex_bar, tank_bar)

#create function to color label nodes to cue type 
lab_col <- as.numeric(methAllData21L[,"cue"])
lab_col <- lab_col[order.dendrogram(cpg.den)] #sort based on order in dendrogram
lab_col_use <- ifelse(lab_col == 1, "red", "black")
labels_colors(cpg.den) <- lab_col_use

cpg.den %>% set("branches_k_color", value = c("darkgrey", "black"), k = 2) %>% plot()
colored_bars(cpg.den, colors = bars, rowLabels = c("Sex", "Tank"), sort_by_labels_order = TRUE)

#or another way to visualize 
labels_colors(cpg.den) <- "white"
cpg.den %>% set("branches_k_color", value = c("darkgrey", "black"), k = 2) %>% set("leaves_pch", 16) %>%
  set("leaves_col", lab_col_use) %>% set("labels_cex", 0.1) %>% plot()
colored_bars(cpg.den, colors = bars, rowLabels = c("Sex", "Tank"), sort_by_labels_order = TRUE)

##Run clustering analysis for DMS ##
#run PCA on all DMS
DMSPCA21L <- prcomp(diffMeth21L[,2:120], scale. = TRUE, center = TRUE)

#check for impact of sex on first 6 PCs
sex.pca.plot.PC1.2 <- fviz_pca_ind(DMSPCA21L,
                                   col.ind = diffMeth21L$sex,
                                   palette = c("#00AFBB",  "#FC4E07"),
                                   addEllipses = TRUE, 
                                   ellipse.type = "confidence",
                                   legend.title = "Sex",
                                   repel = TRUE     # Avoid text overlapping
)
sex.pca.plot.PC1.2

sex.pca.plot.PC3.4 <- fviz_pca_ind(DMSPCA21L,
                                   col.ind = diffMeth21L$sex,
                                   palette = c("#00AFBB",  "#FC4E07"),
                                   addEllipses = TRUE, 
                                   ellipse.type = "confidence",
                                   legend.title = "Sex",
                                   repel = TRUE, # Avoid text overlapping
                                   axes = c(3,4)
)
sex.pca.plot.PC3.4

sex.pca.plot.PC5.6 <- fviz_pca_ind(DMSPCA21L,
                                   col.ind = diffMeth21L$sex,
                                   palette = c("#00AFBB",  "#FC4E07"),
                                   addEllipses = TRUE, 
                                   ellipse.type = "confidence",
                                   legend.title = "Sex",
                                   repel = TRUE, # Avoid text overlapping
                                   axes = c(5,6)
)
sex.pca.plot.PC5.6 

#check for impact of tank 
tank.pca.plot.PC1.2 <- fviz_pca_ind(methPCA21L,
                                    col.ind = methAllData21L$tank, 
                                    addEllipses = TRUE, 
                                    ellipse.type = "confidence",
                                    legend.title = "Tank",
                                    repel = TRUE     # Avoid text overlapping
)
tank.pca.plot.PC1.2

#extract first 2 PCs
diffMeth21L <- cbind(diffMeth21L, DMSPCA21L$x[,1:2])
colnames(diffMeth21L)[159] = "PC1_DMS"
colnames(diffMeth21L)[160] = "PC2_DMS"

#plot cpg PCA 
DMS.pca.plot <- ggplot(data = diffMeth21L, aes(x = PC1_DMS, y = PC2_DMS, color = tank, shape = cue)) + theme_bw() + geom_point() +
  scale_color_manual(values = c("#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#084594",
                                "#FDD0A2", "#FDAE6B", "#FD8D3C", "#F16913", "#D94801", "#8C2D04")) + labs(color = "Tank", shape = "Cue", lty = "Cue") +
  stat_ellipse(aes(group = cue, lty = cue)) +
  xlab("PC1 (34.2%)") + ylab("PC2 (7.3%)")

DMS.pca.plot

## RUN NMDS 
# create matrix of just perc meth values 
DMS.nmds.dat <- as.matrix(diffMeth21L[,2:120])

#run NMDS with Bray-Curtis 
DMS.nmds = metaMDS(DMS.nmds.dat, distance = "bray", autotransform = FALSE, trymax = 100, maxit = 2000)

#extract NMDS scores
DMS.nmds.scores = as.data.frame(scores(DMS.nmds)$sites)

#check dimensions for projection of NMDS 
cpg.dim.plot <- dimcheckMDS(DMS.nmds.dat) #three dimensions with stress < 0.1 (for when I run on only informative CpGs)
cpg.goodness <- goodness(DMS.nmds)
cpg.stress.plot <- stressplot(DMS.nmds)

#add other data to nmds datafram 
DMS.nmds.scores$fish_ID <- diffMeth21L$fish_ID
DMS.nmds.scores$cue <- diffMeth21L$cue
DMS.nmds.scores$sex <- diffMeth21L$sex
DMS.nmds.scores$tank <- diffMeth21L$tank
DMS.nmds.scores$monthAssay <- diffMeth21L$monthAssay


DMS.nmds.plot <- ggplot(DMS.nmds.scores, aes(x = NMDS1, y = NMDS2, color = tank, shape = cue)) + geom_point() + theme_bw() +
  scale_color_manual(values = c("#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#084594",
                                "#FDD0A2", "#FDAE6B", "#FD8D3C", "#F16913", "#D94801", "#8C2D04")) + 
  labs(color = "Tank", shape = "Cue", lty = "Cue") + stat_ellipse(aes(group = cue, lty = cue))
DMS.nmds.plot

## Determine cluster groups ##
#use average silhoutte method to classify individuals into clusters 
#create dataset of just perc meth 
percMethDMS <- diffMeth21L[,2:120]

#function to compute average silhoutte for k clusters
avg_sil_DMS <- function(k) {
  km.res <- kmeans(percMethDMS, centers = k, nstart = 25)
  ss <- silhouette(km.res$cluster, dist(percMethDMS))
  mean(ss[, 3])
}

#extract avg silouette for 2-12 clusters 
avg_sil_values_DMS <- map_dbl(k.values, avg_sil_DMS)

#plot results
plot(k.values, avg_sil_values_DMS, type = "b", pch = 19, frame = FALSE, xlab = "Number of clusters K", ylab = "Average Silhouettes")

#compute k-means clustering with k = 2
final.cluster.DMS <- kmeans(percMethDMS, 2, nstart = 25, iter.max = 100)
diffMeth21L$dmsCluster <- as.factor(final.cluster.DMS$cluster) #add onto data
DMS.nmds.scores$dmsCluster <- as.factor(final.cluster.DMS$cluster)

## replot PCA and NMDS ##
DMS.nmds.plot.cluster <- ggplot(DMS.nmds.scores, aes(x = NMDS1, y = NMDS2, color = tank, shape = cue)) + geom_point() + theme_bw() +
  scale_color_manual(values = c("#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#084594",
                                "#FDD0A2", "#FDAE6B", "#FD8D3C", "#F16913", "#D94801", "#8C2D04")) + 
  labs(color = "Tank", shape = "Cue", lty = "Cluster") + stat_ellipse(aes(group = dmsCluster, lty = dmsCluster))
DMS.nmds.plot.cluster

DMS.pca.plot.cluster <- ggplot(data = diffMeth21L, aes(x = PC1_DMS, y = PC2_DMS, color = tank, shape = cue)) + theme_bw() + geom_point() +
  scale_color_manual(values = c("#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#084594",
                                "#FDD0A2", "#FDAE6B", "#FD8D3C", "#F16913", "#D94801", "#8C2D04")) + labs(color = "Tank", shape = "Cue", lty = "Cluster") +
  stat_ellipse(aes(group = dmsCluster, lty = dmsCluster)) + xlab("PC1 (34.2%)") + ylab("PC2 (7.3%)")

DMS.pca.plot.cluster

#check variables against the clustering 
DMS.pca.plot.cluster.sex <- ggplot(data = diffMeth21L, aes(x = PC1_DMS, y = PC2_DMS, color = sex, shape = cue)) + theme_bw() + geom_point() +
  labs(color = "Sex", shape = "Cue", lty = "Cluster") +
  stat_ellipse(aes(group = dmsCluster, lty = dmsCluster)) + xlab("PC1") + ylab("PC2")

DMS.pca.plot.cluster.sex

DMS.pca.plot.cluster.month <- ggplot(data = diffMeth21L, aes(x = PC1_DMS, y = PC2_DMS, color = monthAssay, shape = cue)) + theme_bw() + geom_point() +
  labs(color = "Month Assay", shape = "Cue", lty = "Cluster") +
  stat_ellipse(aes(group = dmsCluster, lty = dmsCluster)) + xlab("PC1") + ylab("PC2")

DMS.pca.plot.cluster.month

DMS.pca.plot.cluster.time <- ggplot(data = diffMeth21L, aes(x = PC1_DMS, y = PC2_DMS, color = time, shape = cue)) + theme_bw() + geom_point() +
  labs(color = "Time of Sampling", shape = "Cue", lty = "Cluster") +
  stat_ellipse(aes(group = dmsCluster, lty = dmsCluster)) + xlab("PC1") + ylab("PC2")

DMS.pca.plot.cluster.time

DMS.pca.plot.cluster.shoal <- ggplot(data = diffMeth21L, aes(x = PC1_DMS, y = PC2_DMS, color = quantile_tight_diff, shape = cue)) + theme_bw() + geom_point() +
  labs(color = "Quantile tight shoal", shape = "Cue", lty = "Cluster") +
  stat_ellipse(aes(group = dmsCluster, lty = dmsCluster)) + xlab("PC1") + ylab("PC2")

DMS.pca.plot.cluster.shoal

DMS.pca.plot.cluster.PC1 <- ggplot(data = diffMeth21L, aes(x = PC1_DMS, y = PC2_DMS, color = quantile_PC1, shape = cue)) + theme_bw() + geom_point() +
  labs(color = "Quantile PC1 Openfield", shape = "Cue", lty = "Cluster") +
  stat_ellipse(aes(group = dmsCluster, lty = dmsCluster)) + xlab("PC1") + ylab("PC2")

DMS.pca.plot.cluster.PC1

DMS.pca.plot.cluster.dist <- ggplot(data = diffMeth21L, aes(x = PC1_DMS, y = PC2_DMS, color = quantile_dist, shape = cue)) + theme_bw() + geom_point() +
  labs(color = "Quantile distance travelled", shape = "Cue", lty = "Cluster") +
  stat_ellipse(aes(group = dmsCluster, lty = dmsCluster)) + xlab("PC1") + ylab("PC2")

DMS.pca.plot.cluster.dist

## Run hierarchical clustering ## 
#create dataset of just perc meth but with sample ID as row name 
dms.hclust.dat <- diffMeth21L[,1:120]
rownames(dms.hclust.dat) <- dms.hclust.dat[,1]
dms.hclust.dat[,1] <- NULL

#calculate a Euclidean distance matrix 
dms.hclust.dist <- vegdist(dms.hclust.dat, method = "euclidean")

#fit cluster dendrogram
dms.hclust <- hclust(dms.hclust.dist, method="ward.D2")
dms.den <- as.dendrogram(dms.hclust)

#plot dendrogram
#get color pallette values 
brewer.pal(12, "Set3")

#make bars for color bars 
sex_bar <- ifelse(diffMeth21L$sex == "female", "#c51b7d", "blue")
tanks <- diffMeth21L[,"tank"]
cols_12 <- RColorBrewer::brewer.pal(12, "Paired")
tank_bar <- cols_12[tanks]

#bind bars together
bars <- cbind(sex_bar, tank_bar)

#create function to color label nodes to cue type 
lab_col <- as.numeric(diffMeth21L[,"cue"])
lab_col <- lab_col[order.dendrogram(dms.den)] #sort based on order in dendrogram
lab_col_use <- ifelse(lab_col == 1, "red", "black")
labels_colors(dms.den) <- lab_col_use

dms.den %>% set("branches_k_color", value = c("darkgrey", "black"), k = 2) %>% plot()
colored_bars(dms.den, colors = bars, rowLabels = c("Sex", "Tank"), sort_by_labels_order = TRUE)

#or another way to visualize 
labels_colors(dms.den) <- "white"
dms.den %>% set("branches_k_color", value = c("darkgrey", "black"), k = 2) %>% set("leaves_pch", 16) %>%
  set("leaves_col", lab_col_use) %>% set("labels_cex", 0.1) %>% plot()
colored_bars(dms.den, colors = bars, rowLabels = c("Sex", "Tank"), sort_by_labels_order = TRUE)
