############################################
### Goal: Plot and analyze behavioural data
### Author: Janay Fox
### R script
############################################

## Set up ##
#install packages
#install.packages("ggplot2")
#install.packages("dplyr")
#install.packages("forcats")
#install.packages("reshape2")
#install.packages("factoextra")
#install.packages("tidyverse")
#install.packages("lme4")
#install.packages("ggpubr")
#install.packages("MuMIn")
#install.packages("car")
#install.packages("lmerTest")
#install.packages("r2glmm")

#load packaes
library(ggplot2)
library(dplyr)
library(forcats)
library(reshape2)
library(factoextra)
library(tidyverse)
library(lme4)
library(ggpubr)
library(MuMIn)
library(car)
library(lmerTest)
library(r2glmm)

#read data files
of.data <- read_csv("./dev_exp/data/clean/clean_openfield.csv")
sh.data <- read_csv("./dev_exp/data/clean/clean_shoaling.csv")
wt.data <- read_csv("./dev_exp/data/clean/dev_size.csv") 

#merge datasets
of.data <- merge(of.data, wt.data[,1:4], by = "ID") 
sh.data <- merge(sh.data, wt.data[,1:4], by = "ID") 

#convert character columns to factors
of.data[sapply(of.data, is.character)] <- lapply(of.data[sapply(of.data, is.character)], as.factor)
sh.data[sapply(sh.data, is.character)] <- lapply(sh.data[sapply(sh.data, is.character)], as.factor)

## Check for impact of cue on weight ##
t.test(weight_g ~ cue, data = of.data[of.data$sex == "f",])
t.test(weight_g ~ cue, data = of.data[of.data$sex == "m",])

size.plot <- of.data %>% ggplot(aes(x=sex, y=weight_g, fill=cue)) + 
  geom_boxplot() + labs(y="Weight (g)", x="Sex") + 
  scale_x_discrete(labels=c("f" = "Females", "m" = "Males")) + geom_hline(yintercept=0, linetype='dashed', size = 1) +
  theme(axis.text = element_text(size=12), legend.text = element_text(size=12), axis.title = element_text(size=15)) + 
  stat_compare_means(method="t.test", label = "p.signif") + 
  scale_fill_discrete(name = "Cue", labels = c("ac" = "Alarm", "c" = "Control"))
size.plot

## Analysis of Open Field Data ##
# Exploratory plots #
#Distance travelled plot
dist.plot <-of.data %>% ggplot(aes(x=cue, y=dist_cm, fill=sex)) + theme_bw() +
  geom_boxplot() + labs(y="Distance travelled (cm)", x="Cue") + 
  scale_x_discrete(labels=c("ac" = "Alarm cue", "c" = "Control cue")) + 
  theme(axis.text = element_text(size=12), axis.title = element_text(size=15), legend.text = element_text(size=11)) +
  scale_fill_manual(labels = c("Females", "Males"), name = "Sex", values = c("#FFE17B", "#9d4edd")) 
dist.plot

#Time in shelter plot
shelt.plot <-of.data %>% ggplot(aes(x=cue, y=time_shelter_s, fill=sex)) + theme_bw() +
  geom_boxplot() + labs(y="Time in shelter (s)", x="Cue") + 
  scale_x_discrete(labels=c("ac" = "Alarm cue", "c" = "Control cue")) + 
  theme(axis.text = element_text(size=12), axis.title = element_text(size=15), legend.text = element_text(size=11)) +
  scale_fill_manual(labels = c("Females", "Males"), name = "Sex", values = c("#FFE17B", "#9d4edd")) 
shelt.plot

#Time in outer edge
edge.plot <-of.data %>% ggplot(aes(x=cue, y=time_outer_s, fill=sex)) + theme_bw() +
  geom_boxplot() + labs(y="Time in outer edge (s)", x="Cue") + 
  scale_x_discrete(labels=c("ac" = "Alarm cue", "c" = "Control cue")) + 
  theme(axis.text = element_text(size=12), axis.title = element_text(size=15), legend.text = element_text(size=11)) +
  scale_fill_manual(labels = c("Females", "Males"), name = "Sex", values = c("#FFE17B", "#9d4edd")) 
edge.plot

#Time frozen
frozen.plot <-of.data %>% ggplot(aes(x=cue, y=frozen_s, fill=sex)) + theme_bw() +
  geom_boxplot() + labs(y="Time frozen (s)", x="Cue") + 
  scale_x_discrete(labels=c("ac" = "Alarm cue", "c" = "Control cue")) + 
  theme(axis.text = element_text(size=12), axis.title = element_text(size=15), legend.text = element_text(size=11)) +
  scale_fill_manual(labels = c("Females", "Males"), name = "Sex", values = c("#FFE17B", "#9d4edd")) 
frozen.plot

#scatterplots between variables 
dist.shelt.scatterplot <-of.data %>% ggplot(aes(x=dist_cm, y=time_shelter_s)) + geom_point() + 
  labs(y="Distance Travelled (cm)", x="Time Spent in Shelter (s)") + 
  theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=15)) 
dist.shelt.scatterplot

dist.frozen.scatterplot <-of.data %>% ggplot(aes(x=dist_cm, y=frozen_s)) + geom_point() + 
  labs(y="Distance Travelled (cm)", x="Time Spent in Frozen (s)") + 
  theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=15)) 
dist.frozen.scatterplot

dist.edge.scatterplot <-of.data %>% ggplot(aes(x=dist_cm, y=time_outer_s)) + geom_point() + 
  labs(y="Distance Travelled (cm)", x="Time Spent in Outer Edge (s)") + 
  theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=15)) 
dist.edge.scatterplot

shelt.frozen.scatterplot <-of.data %>% ggplot(aes(x=time_shelter_s, y=frozen_s)) + geom_point() + 
  labs(y="Time spent in shelter (s)", x="Time spent frozen (s)") + 
  theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=15)) 
shelt.frozen.scatterplot

shelt.edge.scatterplot <-of.data %>% ggplot(aes(x=time_shelter_s, y=time_outer_s)) + geom_point() + 
  labs(y="Time spent in shelter (s)", x="Time spent in outer edge (s)") + 
  theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=15)) 
shelt.edge.scatterplot

frozen.edge.scatterplot <-of.data %>% ggplot(aes(x=frozen_s, y=time_outer_s)) + geom_point() + 
  labs(y="Time spent frozen (s)", x="Time spent in outer edge (s)") + 
  theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=15)) 
frozen.edge.scatterplot

#check for correlation between variables 
cor.test(of.data$dist_cm, of.data$time_shelter_s)
cor.test(of.data$dist_cm, of.data$frozen_s)
cor.test(of.data$dist_cm, of.data$time_outer_s)
cor.test(of.data$time_shelter_s, of.data$frozen_s)
cor.test(of.data$time_shelter_s, of.data$time_outer_s)
cor.test(of.data$frozen_s, of.data$time_outer_s)

#run PCA
#for all samples
of.pca <- prcomp(of.data[,c(5,8,11,14)], scale = TRUE)
of.eig <- get_eigenvalue(of.pca) #get eignenvalues
of.var <- get_pca_var(of.pca) #get variance
of.var$contrib

#extract first 3 PCs
of.data <- cbind(of.data, of.pca$x[,1:3])

#plot PCA 
plot.pca <- ggplot(data = of.data, aes(x = PC1, y = PC2, color = cue, shape = sex)) + theme_bw() + geom_point() +
    geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") + 
    stat_ellipse(aes(group = cue)) + xlab("PC1") + ylab("PC2") +
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 11), legend.text = element_text(size=12), 
          legend.title = element_text(size = 13))
plot.pca

# Run LMMs #
#check for need to include tank in models
lm.test <- lm(dist_cm ~ cue * sex, data = of.data)
lm.test.resid <- rstandard(lm.test)
plot(lm.test.resid ~ of.data$tank, xlab = "Tank", ylab = "Standardized residuals")

#check for need to include time of day in models
lm.test <- lm(dist_cm ~ cue * sex, data = of.data)
lm.test.resid <- rstandard(lm.test)
plot(lm.test.resid ~ of.data$time, xlab = "Tank", ylab = "Standardized residuals")

# distance travelled #
#run LMM
dist.lmm <- lmer(dist_cm ~ cue * sex * weight_g + (1|tank), data = of.data, REML = TRUE)

#model validation
#homogeneity of variance
plot(resid(dist.lmm) ~ fitted(dist.lmm), xlab = "Predicted values", ylab = "Normalized residuals")
abline(h = 0, lty = 2)

#independence of model residuals
par(mfrow = c(1, 3), mar = c(4, 4, 0.5, 0.5))
plot(resid(dist.lmm) ~ of.data$weight_g, xlab = "Weight (g)", ylab = "Normalized residuals")
abline(h = 0, lty = 2)

boxplot(resid(dist.lmm) ~ of.data$cue, xlab = "Cue",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)

boxplot(resid(dist.lmm) ~ of.data$sex, xlab = "Sex",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)

#normality of residuals
hist(resid(dist.lmm))

#interpret results 
summary(dist.lmm)

#calculate stats with car package 
car::Anova(dist.lmm, type = 3)

#calulate signifcance of random effects 
ranova(dist.lmm)

#calculate partial R2 
r2beta(dist.lmm, partial = TRUE,  method = 'sgv')

# time frozen #
#run LMM
frozen.lmm <- lmer(frozen_s ~ cue * sex * weight_g + (1|tank), data = of.data, REML = TRUE)

#model validation
#homogeneity of variance
par(mar = c(4, 4, 0.5, 0.5))
plot(resid(frozen.lmm) ~ fitted(frozen.lmm), xlab = "Predicted values", ylab = "Normalized residuals")
abline(h = 0, lty = 2)

#independence of model residuals
par(mfrow = c(1, 3), mar = c(4, 4, 0.5, 0.5))
plot(resid(frozen.lmm) ~ of.data$weight_g, xlab = "Weight (g)", ylab = "Normalized residuals")
abline(h = 0, lty = 2)

boxplot(resid(frozen.lmm) ~ of.data$cue, xlab = "Cue",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)

boxplot(resid(frozen.lmm) ~ of.data$sex, xlab = "Sex",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)

#normality of residuals
hist(resid(frozen.lmm))

#interpret results 
summary(frozen.lmm)

#calculate stats with car package 
car::Anova(frozen.lmm, type = 3)
car::Anova(frozen.lmm, type = 2)

#calulate signifcance of random effects 
ranova(frozen.lmm)

#calculate partial R2 
r2beta(frozen.lmm, partial = TRUE,  method = 'sgv')

#for time in outer edge
outer.lmm <- lmer(time_outer_s ~ cue * sex * weight_g + (1|tank), data = of.data)

#model validation
#homogeneity of variance
par(mar = c(4, 4, 0.5, 0.5))
plot(resid(outer.lmm) ~ fitted(outer.lmm), xlab = "Predicted values", ylab = "Normalized residuals")
abline(h = 0, lty = 2)

#independence of model residuals
par(mfrow = c(1, 3), mar = c(4, 4, 0.5, 0.5))
plot(resid(outer.lmm) ~ of.data$weight_g, xlab = "Weight (g)", ylab = "Normalized residuals")
abline(h = 0, lty = 2)

boxplot(resid(outer.lmm) ~ of.data$cue, xlab = "Cue",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)

boxplot(resid(outer.lmm) ~ of.data$sex, xlab = "Sex",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)

#normality of residuals
hist(resid(outer.lmm))

#interpret results 
summary(outer.lmm)

#calculate stats with car package 
car::Anova(outer.lmm, type = 3)
car::Anova(outer.lmm, type = 2)

#calulate signifcance of random effects 
ranova(outer.lmm)

#calculate partial R2 
r2beta(outer.lmm, partial = TRUE,  method = 'sgv')

#for time in shelter edge
shelter.lmm <- lmer(time_shelter_s ~ cue * sex * weight_g + (1|tank), data = of.data)

#model validation
#homogeneity of variance
par(mar = c(4, 4, 0.5, 0.5))
plot(resid(shelter.lmm) ~ fitted(shelter.lmm), xlab = "Predicted values", ylab = "Normalized residuals")
abline(h = 0, lty = 2)

#independence of model residuals
par(mfrow = c(1, 3), mar = c(4, 4, 0.5, 0.5))
plot(resid(shelter.lmm) ~ of.data$weight_g, xlab = "Weight (g)", ylab = "Normalized residuals")
abline(h = 0, lty = 2)

boxplot(resid(shelter.lmm) ~ of.data$cue, xlab = "Cue",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)

boxplot(resid(shelter.lmm) ~ of.data$sex, xlab = "Sex",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)

#normality of residuals
hist(resid(shelter.lmm))

#interpret results 
summary(shelter.lmm)

#calculate stats with car package 
car::Anova(shelter.lmm, type = 3)
car::Anova(shelter.lmm, type = 2)

#calulate signifcance of random effects 
ranova(shelter.lmm)

#calculate partial R2 
r2beta(shelter.lmm, partial = TRUE,  method = 'sgv')

## Analysis of shoaling test data ##
#calculate time differences for shoal vs empty container
sh.data$tight_diff_s <- sh.data$tight_shoal_s - sh.data$tight_emp_s
sh.data$lse_diff_s <- sh.data$lse_shoal_s - sh.data$lse_emp_s

# Exploratory plots #
#Time loose shoaling
lse.plot <- sh.data %>% ggplot(aes(x=cue, y=lse_diff_s, fill=sex)) + theme_bw() + 
  geom_boxplot() + labs(y="Preference for loose shoaling (s)", x="Cue") + 
  scale_x_discrete(labels=c("ac" = "Alarm cue", "c" = "Control cue")) + 
  theme(axis.text = element_text(size=12), axis.title = element_text(size=15), legend.text = element_text(size=11)) +
  scale_fill_manual(labels = c("Females", "Males"), name = "Sex", values = c("#FFE17B", "#9d4edd")) 
lse.plot

tight.plot <- sh.data %>% ggplot(aes(x=cue, y=tight_diff_s, fill=sex)) + theme_bw() + 
  geom_boxplot() + labs(y="Preference for tight shoaling (s)", x="Cue") + 
  scale_x_discrete(labels=c("ac" = "Alarm cue", "c" = "Control cue")) + 
  theme(axis.text = element_text(size=12), axis.title = element_text(size=15), legend.text = element_text(size=11)) + 
  scale_fill_manual(labels = c("Females", "Males"), name = "Sex", values = c("#FFE17B", "#9d4edd")) 
tight.plot

#check for correlation between variables 
cor.test(sh.data$lse_diff_s, sh.data$tight_diff_s)

# Run LMMs #
#check for need to include tank in models
lm.test <- lm(tight_diff_s ~ cue * sex, data = sh.data)
lm.test.resid <- rstandard(lm.test)
plot(lm.test.resid ~ sh.data$tank, xlab = "Tank", ylab = "Standardized residuals")

#check for need to include time of day in models
plot(lm.test.resid ~ sh.data$time.x, xlab = "Time of day", ylab = "Standardized residuals")

# tight shoaling #
#run LMM
tight.lmm <- lmer(tight_diff_s ~ cue * sex + (1|tank), data = sh.data, REML = TRUE)

#model validation
#homogeneity of variance
plot(resid(tight.lmm) ~ fitted(tight.lmm), xlab = "Predicted values", ylab = "Normalized residuals")
abline(h = 0, lty = 2)

#independence of model residuals
boxplot(resid(tight.lmm) ~ sh.data$cue, xlab = "Cue",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)

boxplot(resid(tight.lmm) ~ sh.data$sex, xlab = "Sex",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)

#normality of residuals
hist(resid(tight.lmm))

#interpret results 
summary(tight.lmm)

#calculate stats with car package 
car::Anova(tight.lmm, type = 3)
car::Anova(tight.lmm, type = 2)

#calulate signifcance of random effects 
ranova(tight.lmm)

#calculate partial R2 
r2beta(tight.lmm, partial = TRUE,  method = 'sgv')

# loose shoaling #
#run lmm
loose.lmm <- lmer(lse_diff_s ~ cue * sex + (1|tank), data = sh.data, REML = TRUE)

#model validation
#homogeneity of variance
plot(resid(loose.lmm) ~ fitted(loose.lmm), xlab = "Predicted values", ylab = "Normalized residuals")
abline(h = 0, lty = 2)

#independence of model residuals
boxplot(resid(loose.lmm) ~ sh.data$cue, xlab = "Cue",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)

boxplot(resid(loose.lmm) ~ sh.data$sex, xlab = "Sex",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)

#normality of residuals
hist(resid(loose.lmm))

#interpret results 
summary(loose.lmm)

#calculate stats with car package 
car::Anova(loose.lmm, type = 3)
car::Anova(loose.lmm, type = 2)

#calulate signifcance of random effects 
ranova(loose.lmm)

#calculate partial R2 
r2beta(loose.lmm, partial = TRUE,  method = 'sgv')

## Make panels ##
of.panel <- ggarrange(dist.plot, frozen.plot, shelt.plot, edge.plot, labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")
of.panel

shoal.panel <- ggarrange(lse.plot, tight.plot, labels = c("A", "B"), ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom")
shoal.panel

#save panels
tiff("openField_panel.tiff", units="in", width = 7, height = 7, res = 600)
of.panel
dev.off()

tiff("shoaling_panel.tiff", units="in", width = 7, height = 5, res = 600)
shoal.panel
dev.off()

