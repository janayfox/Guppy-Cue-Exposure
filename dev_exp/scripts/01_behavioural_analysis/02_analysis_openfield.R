############################################
### Goal: Plot and analyze open field data
### Author: Janay Fox
### R script
############################################

## Set up ##
#install packages
#install.packages("ggplot2")
#install.packages("dplyr")
#install.packages("here")
#install.packages("forcats")
#install.packages("reshape2")
#install.packages("factoextra")
#install.packages("tidyverse")
#install.packages("lme4")
#install.packages("lmerTest")
#install.packages("ggpubr")
#install.packages("MuMIn")
#install.packages("DHARMa")

#load packaes
library(ggplot2)
library(dplyr)
library(here)
library(forcats)
library(reshape2)
library(factoextra)
library(tidyverse)
library(lme4)
library(lmerTest)
library(ggpubr)
library(MuMIn)
library(DHARMa)

#read data files
data <- read_csv(here("gup_cue_exp", "data", "clean", "clean_openfield.csv"))

#convert character columns to factors
data[sapply(data, is.character)] <- lapply(data[sapply(data, is.character)], as.factor)

## Exploratory plots ##
#Distance travelled plot
dist.plot <-data %>% ggplot(aes(x=cue, y=dist_cm, fill=sex)) + 
  geom_boxplot() + labs(y="Distance travelled (cm)", x="Cue") + 
  scale_x_discrete(labels=c("ac" = "Alarm cue", "c" = "Control cue")) + 
  theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=15)) + 
  scale_fill_discrete(labels = c("Females", "Males"))
dist.plot

#Time in shelter plot
shelt.plot <-data %>% ggplot(aes(x=cue, y=time_shelter_s, fill=sex)) + 
  geom_boxplot() + labs(y="Time spent in shelter (s)", x="Cue") + 
  scale_x_discrete(labels=c("ac" = "Alarm cue", "c" = "Control cue")) + 
  theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=15)) + 
  scale_fill_discrete(labels = c("Females", "Males"))
shelt.plot

#Time in outer edge
edge.plot <-data %>% ggplot(aes(x=cue, y=time_outer_s, fill=sex)) + 
  geom_boxplot() + labs(y="Time spent in outer edge (s)", x="Cue") + 
  scale_x_discrete(labels=c("ac" = "Alarm cue", "c" = "Control cue")) + 
  theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=15)) + 
  scale_fill_discrete(labels = c("Females", "Males"))
edge.plot

#Time frozen
frozen.plot <-data %>% ggplot(aes(x=cue, y=frozen_s, fill=sex)) + 
  geom_boxplot() + labs(y="Time spent frozen (s)", x="Cue") + 
  scale_x_discrete(labels=c("ac" = "Alarm cue", "c" = "Control cue")) + 
  theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=15)) + 
  scale_fill_discrete(labels = c("Females", "Males"))
frozen.plot

### Run PCAs ###

#subset data into sex specific datasets
fem.data <- subset(data, sex == "f")
mal.data <- subset(data, sex == "m")

#run PCA
#for all samples
of.pca <- prcomp(data[,c(5,8,11,14)], scale = TRUE)
of.eig <- get_eigenvalue(of.pca) #get eignenvalues
of.var <- get_pca_var(of.pca) #get variance
of.var$contrib
cue.groups <- as.factor(data$cue) #group by cue
sex.groups <- as.factor(data$sex) #also want to see groups by sex
time.groups <- as.factor(data$time) #check for grouping by time of behavioural assay

#plot using cue groups
cue.pca.plot <- fviz_pca_ind(of.pca,
                             col.ind = cue.groups, # Color by the quality of representation
                             palette = c("#00AFBB",  "#FC4E07"),
                             addEllipses = TRUE, # Concentration ellipses
                             ellipse.type = "confidence",
                             legend.title = "Cue",
                             repel = TRUE     # Avoid text overlapping
)
cue.pca.plot

#plot using sex groups
sex.pca.plot <- fviz_pca_ind(of.pca,
                             col.ind = sex.groups, # Color by the quality of representation
                             palette = c("#00AFBB",  "#FC4E07"),
                             addEllipses = TRUE, # Concentration ellipses
                             ellipse.type = "confidence",
                             legend.title = "Sex",
                             repel = TRUE     # Avoid text overlapping
)
sex.pca.plot

#plot using time groups
time.pca.plot <- fviz_pca_ind(of.pca,
                             col.ind = time.groups, # Color by the quality of representation
                             palette = c("#00AFBB",  "#FC4E07"),
                             addEllipses = TRUE, # Concentration ellipses
                             ellipse.type = "confidence",
                             legend.title = "Cue",
                             repel = TRUE     # Avoid text overlapping
)
time.pca.plot
#no grouping by time of behavioural trial

# plot variance 
of.pca.var.plot <- fviz_pca_var(of.pca,
                                 col.var = "contrib", # Color by contributions to the PC
                                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                 repel = TRUE     # Avoid text overlapping
)
of.pca.var.plot

#extract first 3 PCs
data <- cbind(data, of.pca$x[,1:3])

### run LMMs ###
wt.data <- read_csv(here("gup_cue_exp", "data", "clean", "dev_size.csv")) #import weight data 
data <- merge(data, wt.data[,1:4], by = "ID") #merge onto dataset

#check for impact of cue on weight
size.plot <- data %>% ggplot(aes(x=sex, y=weight_g, fill=cue)) + 
  geom_boxplot() + labs(y="Weight (g)", x="Sex") + 
  scale_x_discrete(labels=c("f" = "Females", "m" = "Males")) + geom_hline(yintercept=0, linetype='dashed', size = 1) +
  theme(axis.text = element_text(size=12), legend.text = element_text(size=12), axis.title = element_text(size=15)) + 
  stat_compare_means(method="t.test", label = "p.signif") + 
  scale_fill_discrete(name = "Cue", labels = c("ac" = "Alarm", "c" = "Control"))
size.plot

#for PC1
#check for need to include tank as a random effect
lm.test <- lm(PC1 ~ cue * sex * weight_g, data = data)
lm.test.resid <- rstandard(lm.test)
plot(lm.test.resid ~ data$tank, xlab = "Tank", ylab = "Standardized residuals")


#run LMMs
PC1.lmm <- lmer(PC1 ~ cue * sex * weight_g + (1|tank), data = data)
PC1.lmm.nowt <- lmer(PC1 ~ cue * sex + (1|tank), data = data)
PC1.lmm.cue <- lmer(PC1 ~ cue + (1|tank), data = data)

#compare models
AIC.table.PC1 <- model.sel(PC1.lmm, PC1.lmm.nowt, PC1.lmm.cue)
(AIC.table.PC1 <- AIC.table.PC1[, c("df", "logLik", "AICc", "delta")])

#model validation
#homogeneity of variance
par(mar = c(4, 4, 0.5, 0.5))
plot(resid(PC1.lmm) ~ fitted(PC1.lmm), xlab = "Predicted values", ylab = "Normalized residuals")
abline(h = 0, lty = 2)

#independence of model residuals
par(mfrow = c(1, 3), mar = c(4, 4, 0.5, 0.5))
plot(resid(PC1.lmm) ~ data$weight_g, xlab = "Weight (g)", ylab = "Normalized residuals")
abline(h = 0, lty = 2)

boxplot(resid(PC1.lmm) ~ data$cue, xlab = "Cue",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)

boxplot(resid(PC1.lmm) ~ data$sex, xlab = "Sex",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)

#normality of residuals
hist(resid(PC1.lmm))

# Residual diagnostics
simulationOutput <- simulateResiduals(fittedModel = PC1.lmm)
plot(simulationOutput)

#interpret results 
summary(PC1.lmm)

#for PC2
#check for need to include tank as a random effect
lm.test <- lm(PC2 ~ cue * sex * weight_g, data = data)
lm.test.resid <- rstandard(lm.test)
plot(lm.test.resid ~ data$tank, xlab = "Tank", ylab = "Standardized residuals")

#run LMMs
PC2.lmm <- lmer(PC2 ~ cue * sex * weight_g + (1|tank), data = data)
PC2.lmm.nowt <- lmer(PC2 ~ cue * sex + (1|tank), data = data)
PC2.lmm.cue <- lmer(PC2 ~ cue + (1|tank), data = data)

#compare models
AIC.table.PC2 <- model.sel(PC2.lmm, PC2.lmm.nowt, PC2.lmm.cue)
(AIC.table.PC2 <- AIC.table.PC2[, c("df", "logLik", "AICc", "delta")])

#model validation
#homogeneity of variance
par(mar = c(4, 4, 0.5, 0.5))
plot(resid(PC2.lmm) ~ fitted(PC2.lmm), xlab = "Predicted values", ylab = "Normalized residuals")
abline(h = 0, lty = 2)

#independence of model residuals
par(mfrow = c(1, 3), mar = c(4, 4, 0.5, 0.5))
plot(resid(PC2.lmm) ~ data$weight_g, xlab = "Weight (g)", ylab = "Normalized residuals")
abline(h = 0, lty = 2)

boxplot(resid(PC2.lmm) ~ data$cue, xlab = "Cue",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)

boxplot(resid(PC2.lmm) ~ data$sex, xlab = "Sex",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)

#normality of residuals
hist(resid(PC2.lmm))

# Residual diagnostics
simulationOutput <- simulateResiduals(fittedModel = PC2.lmm)
plot(simulationOutput)

#interpret results 
summary(PC2.lmm)

#for PC3
PC3.lmm <- lmer(PC3 ~ cue * sex * weight_g + (1|tank), data = data)
PC3.lmm.nowt <- lmer(PC3 ~ cue * sex + (1|tank), data = data)
PC3.lmm.cue <- lmer(PC3 ~ cue + (1|tank), data = data)

#compare models
AIC.table.PC3 <- model.sel(PC3.lmm, PC3.lmm.nowt, PC3.lmm.cue)
(AIC.table.PC3 <- AIC.table.PC3[, c("df", "logLik", "AICc", "delta")])

#model validation
#homogeneity of variance
par(mar = c(4, 4, 0.5, 0.5))
plot(resid(PC3.lmm) ~ fitted(PC3.lmm), xlab = "Predicted values", ylab = "Normalized residuals")
abline(h = 0, lty = 2)

#independence of model residuals
par(mfrow = c(1, 3), mar = c(4, 4, 0.5, 0.5))
plot(resid(PC3.lmm) ~ data$weight_g, xlab = "Weight (g)", ylab = "Normalized residuals")
abline(h = 0, lty = 2)

boxplot(resid(PC3.lmm) ~ data$cue, xlab = "Cue",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)

boxplot(resid(PC3.lmm) ~ data$sex, xlab = "Sex",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)

#normality of residuals
hist(resid(PC3.lmm))

# Residual diagnostics
simulationOutput <- simulateResiduals(fittedModel = PC3.lmm)
plot(simulationOutput)

#interpret results 
summary(PC3.lmm)

# test doing just distance
dist.lmm <- lmer(dist_cm ~ cue * sex * weight_g + (1|tank), data = data)
dist.lmm.nowt <- lmer(dist_cm ~ cue * sex + (1|tank), data = data)
dist.lmm.cue <- lmer(dist_cm ~ cue + (1|tank), data = data)

AIC.table.dist <- model.sel(dist.lmm, dist.lmm.nowt, dist.lmm.cue)
(AIC.table.dist <- AIC.table.dist[, c("df", "logLik", "AICc", "delta")])

simulationOutput <- simulateResiduals(fittedModel = dist.lmm)
plot(simulationOutput)

summary(dist.lmm)

froz.lmm <- lmer(frozen_s ~ cue * sex * weight_g + (1|tank), data = data)
froz.lmm.nowt <- lmer(frozen_s ~ cue * sex + (1|tank), data = data)
froz.lmm.cue <- lmer(frozen_s ~ cue + (1|tank), data = data)

AIC.table.froz <- model.sel(froz.lmm, froz.lmm.nowt, froz.lmm.cue)
(AIC.table.froz <- AIC.table.froz[, c("df", "logLik", "AICc", "delta")])

simulationOutput <- simulateResiduals(fittedModel = froz.lmm)
plot(simulationOutput)

summary(froz.lmm)

shelt.lmm <- lmer(time_shelter_s ~ cue * sex * weight_g + (1|tank), data = data)
shelt.lmm.nowt <- lmer(time_shelter_s ~ cue * sex + (1|tank), data = data)
shelt.lmm.cue <- lmer(time_shelter_s ~ cue + (1|tank), data = data)

AIC.table.shelt <- model.sel(shelt.lmm, shelt.lmm.nowt, shelt.lmm.cue)
(AIC.table.shelt <- AIC.table.shelt[, c("df", "logLik", "AICc", "delta")])

simulationOutput <- simulateResiduals(fittedModel = shelt.lmm)
plot(simulationOutput)

summary(shelt.lmm)

out.lmm <- lmer(time_outer_s ~ cue * sex * weight_g + (1|tank), data = data)
out.lmm.nowt <- lmer(time_outer_s ~ cue * sex + (1|tank), data = data)
out.lmm.cue <- lmer(time_outer_s ~ cue + (1|tank), data = data)

AIC.table.out <- model.sel(out.lmm, out.lmm.nowt, out.lmm.cue)
(AIC.table.out <- AIC.table.out[, c("df", "logLik", "AICc", "delta")])

simulationOutput <- simulateResiduals(fittedModel = out.lmm)
plot(simulationOutput)

summary(out.lmm)

## Export plots and metadata ##
ggsave(here("gup_cue_exp", "plots", "openfield", "of_distplot.jpg"), plot = dist.plot)
ggsave(here("gup_cue_exp","plots", "openfield", "of_sheltplot.jpg"), plot = shelt.plot)
ggsave(here("gup_cue_exp","plots", "openfield", "of_edgeplot.jpg"), plot = edge.plot)
ggsave(here("gup_cue_exp","plots", "openfield", "of_frozenplot.jpg"), plot = frozen.plot)
ggsave(here("gup_cue_exp","plots", "openfield", "of_cue_pca_plot.jpg"), plot = cue.pca.plot)
ggsave(here("gup_cue_exp","plots", "openfield", "of_sex_pca_plot.jpg"), plot = sex.pca.plot)
ggsave(here("gup_cue_exp","plots", "openfield", "of_pca_var_plot.jpg"), plot = of.pca.var.plot)

write_csv(data, file(here("gup_cue_exp", "data", "metadata", "of_metadata.csv")))

