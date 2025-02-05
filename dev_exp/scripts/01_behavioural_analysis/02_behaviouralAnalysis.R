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
library(DHARMa)
library(ggtext)

setwd("/Users/janayfox/R_projects/guppy/gup_cue_exp")

#read data files
of.data <- read.csv("./dev_exp/data/clean/clean_openfield_gupEpi_2024.csv")
sh.data <- read.csv("./dev_exp/data/clean/clean_shoaling_gupEpi_2024.csv")
wt.data <- read.csv("./dev_exp/data/clean/clean_size_gupEpi_2024.csv") 

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
#create new boldness variable 
of.data$boldness_score <- of.data$time_shelter_s + of.data$time_outer_s 

#create new exploration variable
count_data <- of.data[,17:44]
cumDur_data <- of.data[,73:100]

unique_data <- count_data %>%
  rowwise() %>%
  mutate(unique_squares = sum(c_across(starts_with("S")) > 2)) %>%
  ungroup()

unique_data_cumDur <- cumDur_data %>%
  rowwise() %>%
  mutate(unique_squares = sum(c_across(starts_with("S")) > 3)) %>%
  ungroup()

of.data$uni_sq <- unique_data$unique_squares
of.data$uni_sq_cumDur <- unique_data_cumDur$unique_squares

#check correlations 
cor(of.data[,c(5,7,8,11,14,105:107)])

# Exploratory plots #
#Distance travelled plot
dist.plot <-of.data %>% ggplot(aes(x=cue, y=dist_cm, fill=sex)) + theme_bw() +
  geom_boxplot() + labs(y=paste0("<span style='font-size: 16pt'>Activity</span><br><br><span style='font-size: 13pt'>Distance travelled (cm)</span>"), x="Cue") + 
  scale_x_discrete(labels=c("ac" = "Alarm cue", "c" = "Control cue")) + 
  theme(axis.text = element_text(size=12), axis.title.y = ggtext::element_markdown(), legend.text = element_text(size=13),
        axis.title.x = element_text(size=13),legend.title = element_text(size=14)) + 
  scale_fill_manual(labels = c("Females", "Males"), name = "Sex", values = c("#FFE17B", "skyblue")) 
dist.plot

#Time in shelter of frozen plot
shelt.plot <-of.data %>% ggplot(aes(x=cue, y=boldness_score, fill=sex)) + theme_bw() +
  geom_boxplot() + labs(y=paste0("<span style='font-size: 16pt'>Boldness</span><br><br><span style='font-size: 13pt'>Time in shelter or edge (s)</span>"), x="Cue") + 
  scale_x_discrete(labels=c("ac" = "Alarm cue", "c" = "Control cue")) + 
  theme(axis.text = element_text(size=12), axis.title.y = ggtext::element_markdown(), legend.text = element_text(size=13),
        axis.title.x = element_text(size=13),legend.title = element_text(size=14)) + 
  scale_fill_manual(labels = c("Females", "Males"), name = "Sex", values = c("#FFE17B", "skyblue")) 
shelt.plot

#Number of unique squares visited
sq.plot <-of.data %>% ggplot(aes(x=cue, y=uni_sq_cumDur, fill=sex)) + theme_bw() +
  geom_boxplot() + labs(y=paste0("<span style='font-size: 16pt'>Exploration</span><br><br><span style='font-size: 13pt'>Number of squares explored</span>"), x="Cue") + 
  scale_x_discrete(labels=c("ac" = "Alarm cue", "c" = "Control cue")) + 
  theme(axis.text = element_text(size=12), axis.title.y = ggtext::element_markdown(), legend.text = element_text(size=13),
        axis.title.x = element_text(size=13),legend.title = element_text(size=14)) + 
  scale_fill_manual(labels = c("Females", "Males"), name = "Sex", values = c("#FFE17B", "skyblue")) 
sq.plot

#run PCA
#for all samples
of.pca <- prcomp(of.data[,c(5,105,107)], scale = TRUE)
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

autoplot(of.pca, data = of.data, 
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3)

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

#interpret results 
summary(dist.lmm)

#calculate stats with car package 
car::Anova(dist.lmm, type = 3)

#DHARMa validation 
simOutput_dist <- simulateResiduals(fittedModel = dist.lmm)
plot(simOutput_dist)

#remove interaction for cue 
dist.lmm <- lmer(dist_cm ~ cue + sex * weight_g + (1|tank), data = of.data, REML = TRUE)

#interpret results 
summary(dist.lmm)
car::Anova(dist.lmm, type = 3)

#calculate partial R2 
r2beta(dist.lmm, partial = TRUE,  method = 'nsj')

# boldness #
#run LMM
bold.lmm <- lmer(boldness_score ~ cue * sex * weight_g + (1|tank), data = of.data, REML = TRUE)

#DHARMa validation 
simOutput_bold <- simulateResiduals(fittedModel = bold.lmm)
plot(simOutput_bold)

#interpret results 
summary(bold.lmm)

#calculate stats with car package 
car::Anova(bold.lmm, type = 3)

#rerun without interactions 
bold.lmm.2 <- lmer(boldness_score ~ cue + sex + weight_g + (1|tank), data = of.data, REML = TRUE)

#interpret results 
summary(bold.lmm.2)
car::Anova(bold.lmm.2, type = 2)

#calculate partial R2 
r2beta(bold.lmm.2, partial = TRUE,  method = 'nsj')

# exploration #
#run glm
exp.glm <- glmer(uni_sq_cumDur ~ cue * sex * weight_g + (1|tank), data = of.data, family = poisson)

#DHARMa validation 
simOutput_exp_glm <- simulateResiduals(fittedModel = exp.glm)
plot(simOutput_exp_glm)

#interpret results 
summary(exp.glm)

#calculate stats with car package 
car::Anova(exp.glm, type = 3)

#rerun without interactions 
exp.glm.2 <- glmer(uni_sq_cumDur ~ cue + weight_g + sex + (1|tank), data = of.data, family = poisson)

#interpret results 
summary(exp.glm.2)
car::Anova(exp.glm.2, type = 3)

#calculate partial R2 
r2beta(exp.glm.2, partial = TRUE,  method = 'nsj')

## Analysis of shoaling test data ##
#calculate total time shoaling 
sh.data$total_shoal_s <- sh.data$tight_shoal_s + sh.data$lse_shoal_s
sh.data$total_emp_s <- sh.data$tight_emp_s + sh.data$lse_emp_s

#calculate time differences for shoal vs empty container
sh.data$tight_diff_s <- sh.data$tight_shoal_s - sh.data$tight_emp_s
sh.data$lse_diff_s <- sh.data$lse_shoal_s - sh.data$lse_emp_s
sh.data$total_diff_s <- sh.data$total_shoal_s - sh.data$total_emp_s

# Exploratory plots #
# total.plot <- sh.data %>% ggplot(aes(x=cue, y=total_diff_s, fill=sex)) + theme_bw() + 
#   geom_boxplot() + labs(y=expression(atop("Shoaling", "Preference for shoaling (s)")), x="Cue") + 
#   scale_x_discrete(labels=c("ac" = "Alarm cue", "c" = "Control cue")) + 
#   theme(axis.text = element_text(size=13), axis.title = element_text(size=15), legend.text = element_text(size=13),
#         legend.title = element_text(size=14)) + 
#   scale_fill_manual(labels = c("Females", "Males"), name = "Sex", values = c("#FFE17B", "skyblue")) 
# total.plot

total.plot <- sh.data %>% ggplot(aes(x=cue, y=total_diff_s, fill=sex)) + theme_bw() + 
  geom_boxplot() + labs(y=paste0("<span style='font-size: 16pt'>Shoaling</span><br><br><span style='font-size: 13pt'>Preference for shoal (s)</span>"), x="Cue") + 
  scale_x_discrete(labels=c("ac" = "Alarm cue", "c" = "Control cue")) + 
  theme(axis.text = element_text(size=13), axis.title.y = ggtext::element_markdown(), legend.text = element_text(size=13), 
        axis.title.x = element_text(size=13),legend.title = element_text(size=14)) + 
  scale_fill_manual(labels = c("Females", "Males"), name = "Sex", values = c("#FFE17B", "skyblue")) 
total.plot

#check for correlation between variables 
cor.test(sh.data$lse_diff_s, sh.data$tight_diff_s)

# Run LMMs #
#check for need to include tank in models
lm.test <- lm(tight_diff_s ~ cue * sex, data = sh.data)
lm.test.resid <- rstandard(lm.test)
plot(lm.test.resid ~ sh.data$tank, xlab = "Tank", ylab = "Standardized residuals")

#check for need to include time of day in models
plot(lm.test.resid ~ sh.data$time.x, xlab = "Time of day", ylab = "Standardized residuals")

#run lmm
combined.lmm <- lmer(total_diff_s ~ cue * sex * weight_g + (1|tank), data = sh.data, REML = TRUE)

#DHARMa validation 
simOutput_sh <- simulateResiduals(fittedModel = combined.lmm)
plot(simOutput_sh)

#interpret results 
summary(combined.lmm)

#calculate stats with car package 
car::Anova(combined.lmm, type = 3)

#rerun without interactions and weight
combined.lmm.2 <- lmer(total_diff_s ~ cue + sex + (1|tank), data = sh.data, REML = TRUE)

#interpret results 
summary(combined.lmm.2)
car::Anova(combined.lmm.2, type = 3)

#calculate partial R2 
r2beta(combined.lmm.2, partial = TRUE,  method = 'nsj')

## Make panels ##
of.panel <- ggarrange(dist.plot, frozen.plot, shelt.plot, edge.plot, labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")
of.panel

shoal.panel <- ggarrange(lse.plot, tight.plot, labels = c("A", "B"), ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom")
shoal.panel

mix.panel <- ggarrange(dist.plot, shelt.plot, sq.plot, total.plot,  labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")
mix.panel

#save panels
tiff("mixed_panel_v3.tiff", units="in", width = 7, height = 8, res = 600)
mix.panel
dev.off()

tiff("openField_panel.tiff", units="in", width = 7, height = 7, res = 600)
of.panel
dev.off()

tiff("shoaling_panel.tiff", units="in", width = 7, height = 5, res = 600)
shoal.panel
dev.off()

ggsave(plot = total.plot, filename =  "total_shoaling.tiff", width = 5, height = 5, units = "in", dpi = 600)

## Check for preference of shoal container ## 
sh.data %>% filter(sex == 'f' & cue == 'ac') %>% with(t.test(total_shoal_s, total_emp_s, paired = TRUE))
sh.data %>% filter(sex == 'f' & cue == 'c') %>% with(t.test(total_shoal_s, total_emp_s, paired = TRUE))
sh.data %>% filter(sex == 'm' & cue == 'ac') %>% with(t.test(total_shoal_s, total_emp_s, paired = TRUE))
sh.data %>% filter(sex == 'm' & cue == 'c') %>% with(t.test(total_shoal_s, total_emp_s, paired = TRUE))


