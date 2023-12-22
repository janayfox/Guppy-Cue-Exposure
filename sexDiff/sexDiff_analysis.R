#######################################################
### Goal: Plot and analyze sex diff results data
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
library(methylKit)
library(genomation)
library(DescTools)
library(tibble)

# load data 
diffmethChr.dev <- readRDS("./sexDiff/methylKit_res/DMSdiffmethchr5x_dev.RDS")
diffmethChr.st <- readRDS("./sexDiff/methylKit_res/DMSdiffmethchr5x_st.RDS")

diffmeth.dev <- readRDS("./sexDiff/methylKit_res/DMSdiffmeth5x_dev_data.RDS")
diffmeth.st <- readRDS("./sexDiff/methylKit_res/DMSdiffmeth5x_st_data.RDS")

anno.dev <- readRDS("./sexDiff/methylKit_res/DMS_annStats_num_dev.RDS")
anno.st <- readRDS("./sexDiff/methylKit_res/DMS_annStats_num_st.RDS")

cpg.anno.dev <- readRDS("./sexDiff/methylKit_res/CpG_annStats_num_dev.RDS")
cpg.anno.st <- readRDS("./sexDiff/methylKit_res/CpG_annStats_num_st.RDS")

tss.dev <- readRDS("./sexDiff/methylKit_res/DMS_TSS_dev.RDS")
tss.st <- readRDS("./sexDiff/methylKit_res/DMS_TSS_st.RDS")

## Analyse and plot Manhattans ##
#calculate percentage of DMS on sex chromosome
table(diffmeth.dev$chr=="NC_024342.1")
396/(2357+396) * 100

table(diffmeth.st$chr=="NC_024342.1")
1480/(4451+1480) * 100

#rename chromosomes
rename.chr <- function(df){
              df$chr[df$chr == "NC_024331.1"] <- 1
              df$chr[df$chr == "NC_024332.1"] <- 2
              df$chr[df$chr == "NC_024333.1"] <- 3
              df$chr[df$chr == "NC_024334.1"] <- 4
              df$chr[df$chr == "NC_024335.1"] <- 5
              df$chr[df$chr == "NC_024336.1"] <- 6
              df$chr[df$chr == "NC_024337.1"] <- 7
              df$chr[df$chr == "NC_024338.1"] <- 8
              df$chr[df$chr == "NC_024339.1"] <- 9
              df$chr[df$chr == "NC_024340.1"] <- 10
              df$chr[df$chr == "NC_024341.1"] <- 11
              df$chr[df$chr == "NC_024342.1"] <- 12
              df$chr[df$chr == "NC_024343.1"] <- 13
              df$chr[df$chr == "NC_024344.1"] <- 14
              df$chr[df$chr == "NC_024345.1"] <- 15
              df$chr[df$chr == "NC_024346.1"] <- 16
              df$chr[df$chr == "NC_024347.1"] <- 17
              df$chr[df$chr == "NC_024348.1"] <- 18
              df$chr[df$chr == "NC_024349.1"] <- 19
              df$chr[df$chr == "NC_024350.1"] <- 20
              df$chr[df$chr == "NC_024351.1"] <- 21
              df$chr[df$chr == "NC_024352.1"] <- 22
              df$chr[df$chr == "NC_024353.1"] <- 23
              return(df)
}

diffmeth.dev <- rename.chr(diffmeth.dev)
diffmeth.st <- rename.chr(diffmeth.st)

#convert to numeric
diffmeth.dev$chr <- as.numeric(diffmeth.dev$chr)
diffmeth.st$chr <- as.numeric(diffmeth.st$chr)

#add a column with a name for each DMS 
diffmeth.dev$DMS <- paste(diffmeth.dev$chr, diffmeth.dev$start, sep = "_")
diffmeth.st$DMS <- paste(diffmeth.st$chr, diffmeth.st$start, sep = "_")

#compute cumulative position of SNP
dev.don <- diffmeth.dev %>%
  #compute chr size
  group_by(chr) %>% 
  summarise(chr_len=max(start)) %>%
  #calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  dplyr::select(-chr_len) %>%
  #add this info to the initial dataset
  left_join(diffmeth.dev, ., by=c("chr"="chr")) %>%
  #add a cumulative position of each SNP
  arrange(chr, start) %>%
  mutate(DMScum=start+tot)

st.don <- diffmeth.st %>%
  #compute chr size
  group_by(chr) %>% 
  summarise(chr_len=max(start)) %>%
  #calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  dplyr::select(-chr_len) %>%
  #add this info to the initial dataset
  left_join(diffmeth.st, ., by=c("chr"="chr")) %>%
  #add a cumulative position of each SNP
  arrange(chr, start) %>%
  mutate(DMScum=start+tot)

#create axis
axisDev <- dev.don %>% group_by(chr) %>% summarize(center=(max(DMScum) + min(DMScum))/2)
axisSt <- st.don %>% group_by(chr) %>% summarize(center=(max(DMScum) + min(DMScum))/2)

#create lists of colors for chromosomes
color_val <- rep(c("grey", "skyblue"), length.out = 23)
color_val[12] <- "orange"

#make manhattan plots 
plot.manhat <- function(data, axisDF){
ggplot(data, aes(x=DMScum, y=meth.diff)) +
  geom_point(aes(color=as.factor(chr)), alpha=0.8, size=1.3) +
  scale_color_manual(values = color_val) +
  scale_x_continuous(label = axisDF$chr, breaks = axisDF$center) +
  theme_light() + ylab("Change in % methylation") + xlab("Chromosome") +
  geom_hline(yintercept = 0, color = "red", linetype="dashed") +
  theme( 
    legend.position="none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
}

dev.man <- plot.manhat(dev.don, axisDev)
dev.man 

tiff("sexDevManPlot.tiff", units="in", width = 7, height = 4, res = 600)
dev.man
dev.off()

st.man <- plot.manhat(st.don, axisSt)
st.man

tiff("sexStManPlot.tiff", units="in", width = 7, height = 4, res = 600)
st.man
dev.off()

## Analyze annotation stats ## 
GTest(anno.dev, p = cpg.anno.dev/(cpg.anno.dev[1]+cpg.anno.dev[2]+cpg.anno.dev[3]+cpg.anno.dev[4]))

GTest(anno.st, p = cpg.anno.st/(cpg.anno.st[1]+cpg.anno.st[2]+cpg.anno.st[3]+cpg.anno.st[4]))$expected

#post hoc tests 
prom.test.dev <- GTest(c(anno.dev[1], (anno.dev[2]+anno.dev[3]+anno.dev[4])), 
                  p = c(cpg.anno.dev[1], (cpg.anno.dev[2]+cpg.anno.dev[3]+cpg.anno.dev[4]))/
                  (cpg.anno.dev[1]+cpg.anno.dev[2]+cpg.anno.dev[3]+cpg.anno.dev[4]))

ex.test.dev <- GTest(c(anno.dev[2], (anno.dev[1]+anno.dev[3]+anno.dev[4])), 
                    p = c(cpg.anno.dev[2], (cpg.anno.dev[1]+cpg.anno.dev[3]+cpg.anno.dev[4]))/
                      (cpg.anno.dev[1]+cpg.anno.dev[2]+cpg.anno.dev[3]+cpg.anno.dev[4]))

int.test.dev <- GTest(c(anno.dev[3], (anno.dev[1]+anno.dev[2]+anno.dev[4])), 
                p = c(cpg.anno.dev[3], (cpg.anno.dev[1]+cpg.anno.dev[2]+cpg.anno.dev[4]))/
                  (cpg.anno.dev[1]+cpg.anno.dev[2]+cpg.anno.dev[3]+cpg.anno.dev[4]))

inter.test.dev <- GTest(c(anno.dev[4], (anno.dev[1]+anno.dev[2]+anno.dev[3])), 
                      p = c(cpg.anno.dev[4], (cpg.anno.dev[1]+cpg.anno.dev[2]+cpg.anno.dev[3]))/
                        (cpg.anno.dev[1]+cpg.anno.dev[2]+cpg.anno.dev[3]+cpg.anno.dev[4]))


prom.test.st <- GTest(c(anno.st[1], (anno.st[2]+anno.st[3]+anno.st[4])), 
                       p = c(cpg.anno.st[1], (cpg.anno.st[2]+cpg.anno.st[3]+cpg.anno.st[4]))/
                         (cpg.anno.st[1]+cpg.anno.st[2]+cpg.anno.st[3]+cpg.anno.st[4]))

ex.test.st <- GTest(c(anno.st[2], (anno.st[1]+anno.st[3]+anno.st[4])), 
                     p = c(cpg.anno.st[2], (cpg.anno.st[1]+cpg.anno.st[3]+cpg.anno.st[4]))/
                       (cpg.anno.st[1]+cpg.anno.st[2]+cpg.anno.st[3]+cpg.anno.st[4]))

int.test.st <- GTest(c(anno.st[3], (anno.st[1]+anno.st[2]+anno.st[4])), 
                      p = c(cpg.anno.st[3], (cpg.anno.st[1]+cpg.anno.st[2]+cpg.anno.st[4]))/
                        (cpg.anno.st[1]+cpg.anno.st[2]+cpg.anno.st[3]+cpg.anno.st[4]))

inter.test.st <- GTest(c(anno.st[4], (anno.st[1]+anno.st[2]+anno.st[3])), 
                        p = c(cpg.anno.st[4], (cpg.anno.st[1]+cpg.anno.st[2]+cpg.anno.st[3]))/
                          (cpg.anno.st[1]+cpg.anno.st[2]+cpg.anno.st[3]+cpg.anno.st[4]))

#adjsut p values 
dev.gtest.p <- c(prom.test.dev$p.value,ex.test.dev$p.value,int.test.dev$p.value,inter.test.dev$p.value)
dev.gtest.p <- p.adjust(dev.gtest.p, method = "hommel")

st.gtest.p <- c(prom.test.st$p.value,ex.test.st$p.value,int.test.st$p.value,inter.test.st$p.value)
st.gtest.p <- p.adjust(st.gtest.p, method = "hommel")

#format data for plotting 
anno.data.dev <- data.frame(anno.dev)
anno.data.dev$dist <- "DMS"
colnames(anno.data.dev)[1] <- "count"
anno.data.dev <- rownames_to_column(anno.data.dev, "feature")

cpg.anno.data.dev <- data.frame(cpg.anno.dev)
cpg.anno.data.dev$dist <- "CpG"
colnames(cpg.anno.data.dev)[1] <- "count"
cpg.anno.data.dev <- rownames_to_column(cpg.anno.data.dev, "feature")

dev.plot.data <- rbind(anno.data.dev, cpg.anno.data.dev)

anno.data.st <- data.frame(anno.st)
anno.data.st$dist <- "DMS"
colnames(anno.data.st)[1] <- "count"
anno.data.st <- rownames_to_column(anno.data.st, "feature")

cpg.anno.data.st <- data.frame(cpg.anno.st)
cpg.anno.data.st$dist <- "CpG"
colnames(cpg.anno.data.st)[1] <- "count"
cpg.anno.data.st <- rownames_to_column(cpg.anno.data.st, "feature")

st.plot.data <- rbind(anno.data.st, cpg.anno.data.st)

bar.plot.func <- function(data){
          ggplot(data, aes(x = dist, y = count, fill = feature)) +
            geom_col(colour = "black", position = "fill") + scale_y_continuous(labels = scales::percent) + theme_bw() + xlab(NULL) +
            ylab("Percent") + scale_fill_manual(values = c("exon" = "lightgreen", "intron" = "skyblue", 
                                                           "intergenic" = "orange", "promoter" = "yellow"),
                                                name = "Feature") + 
            theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12), legend.text = element_text(size = 12),
                  legend.title = element_text(size = 14)) 

}

dev.barplot <- bar.plot.func(dev.plot.data)
dev.barplot

tiff("devBarPlot.tiff", units="in", width = 4, height = 5, res = 600)
dev.barplot
dev.off()

st.barplot <- bar.plot.func(st.plot.data)
st.barplot

tiff("stBarPlot.tiff", units="in", width = 4, height = 5, res = 600)
st.barplot
dev.off()
