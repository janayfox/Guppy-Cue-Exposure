#######################################################
### Goal: Plot and Analyze cue response data
### Author: Janay Fox
### R script
#######################################################

## Set up ##
#install packages
#install.packages("ggplot2")
#install.packages("dplyr")
#install.packages("here")
#install.packages("forcats")
#install.packages("reshape2")
#install.packages("tidyverse")
#install.packages("lme4")
#install.packages("lmerTest")
#install.packages("ggpubr")
#install.packages("RLRsim")

#load packaes
library(ggplot2)
library(dplyr)
library(here)
library(forcats)
library(reshape2)
library(tidyverse)
library(lme4)
library(lmerTest)
library(ggpubr)
library(RLRsim)

#read data files
cue.data <- read_csv("./shortTerm_exp/data/clean/clean_cue_fivemin.csv")

##Prepare data##
#replace NAs with 0s 
cue.data[is.na(cue.data)] <- 0

#convert character columns to factors
cue.data[sapply(cue.data, is.character)] <- lapply(cue.data[sapply(cue.data, is.character)], as.factor)

#order factors
cue.data$time <- factor(cue.data$time, levels=c('before', 'after'))

#remove tank that doesnt have after data 
cue.data <- subset(cue.data, tank != "ST2AC13")

## Calculate behavioural measurements ##
#calculate time fish spent at bottom of tank (in Q5 or Q6)
cue.data$bottom_s <- (cue.data$Q5_s + cue.data$Q6_s)

#calculate substrate use (time at bottom of tank minus time spent foraging)
cue.data$subuse_s <- (cue.data$bottom_s - cue.data$feeding_s)

#convert substrate use to a proportion
cue.data$subuse_prop <- cue.data$subuse_s/cue.data$length_obs_s

#convert time frozen to proportion
cue.data$frozen_prop <- cue.data$frozen_s/cue.data$length_obs_s

#combine mating data 
cue.data$matbehav_s <- cue.data$sig_s + cue.data$chasing_s

#convert mating data to proportion
cue.data$matbehav_prop <- cue.data$matbehav_s/cue.data$length_obs_s

#calculate change in substrate use
diff.sub <- cue.data %>% dcast(fish_ID ~ time, value.var = "subuse_prop", fill=0) %>%
  mutate(diff_sub_use = after - before) %>% select(fish_ID, diff_sub_use)

#calculate change in mating behaviour
diff.mating <- cue.data %>% dcast(fish_ID ~ time, value.var = "matbehav_prop", fill=0) %>%
  mutate(diff_mating = after - before) %>% select(fish_ID, diff_mating)

#calculate total mating behaviour and total trial time
total.mating <- cue.data %>% dcast(fish_ID ~ time, value.var = "matbehav_s", fill=0) %>%
  mutate(total_mating = after + before) %>% select(fish_ID, total_mating)

total.trial <- cue.data %>% dcast(fish_ID ~ time, value.var = "length_obs_s", fill=0) %>%
  mutate(total_trial = after + before) %>% select(fish_ID, total_trial)
  
#merge mating dataframes
mating.data <- merge(diff.mating, total.mating, by = "fish_ID")
mating.data <- merge(mating.data, total.trial, by = "fish_ID")

#need to add cue and back onto new dataframes
diff.sub <- distinct(merge(diff.sub, cue.data[, c("fish_ID", "tank", "cue", "sex")], by = "fish_ID")) 
mating.data <- distinct(merge(mating.data, cue.data[, c("fish_ID", "tank", "cue", "sex")], by = "fish_ID")) 

#remove females from mating data 
mating.data <- subset(mating.data, sex == "male")

#calculate prop total mating 
mating.data$total_matingProp <- mating.data$total_mating/mating.data$total_trial

## Make plots ##
#boxplot of proportion of substrate use
subuse.prop.plot <- cue.data %>% ggplot(aes(x=time, y=subuse_prop, fill=time)) + 
  geom_boxplot() + labs(y="Proportion of trial spent near substrate", x="Cue") + geom_point() + geom_line(aes(group = fish_ID)) +
  scale_x_discrete(labels=c("alarm" = "Alarm cue", "control" = "Control cue")) + 
  theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=15)) + 
  scale_fill_discrete(labels = c("Before cue", "After cue")) + facet_wrap(~ cue)
#messy so use change in sub use plot 

#boxplot of change in substrate use
#across all tanks
subusediff.plot <- diff.sub %>% ggplot(aes(x=cue, y=diff_sub_use)) + theme_bw() +
  geom_boxplot(fill = c("#E14D2A", "#31C6D4")) + labs(y="Difference in substrate use", x="Cue") +
  scale_x_discrete(labels=c("alarm" = "Alarm cue", "control" = "Control cue")) + geom_hline(yintercept=0, linetype='dashed', size = 1) +
  theme(legend.position = "none", axis.text = element_text(size=12), axis.title = element_text(size=15)) 

subusediff.plot

ggsave("subuseDiff_plot.tiff", plot = subusediff.plot, device = "tiff", width = 4, height = 4, units = "in")

#with p-val
pval.fivemin.subusediff.plot <- fivemin.diff.sub %>% ggplot(aes(x=sex, y=diff_sub_use, fill=cue)) + 
  geom_boxplot() + labs(y="Difference in substrate use", x="Sex") + 
  scale_x_discrete(labels=c("F" = "Females", "M" = "Males")) + geom_hline(yintercept=0, linetype='dashed', size = 1) +
  theme(axis.text = element_text(size=12), legend.text = element_text(size=12), axis.title = element_text(size=15)) + 
  stat_compare_means(method="t.test", label = "p.signif", size = 12) + 
  scale_fill_discrete(name = "Cue", labels = c("Alarm","Control"))

pval.twomin.subusediff.plot <- twomin.diff.sub %>% ggplot(aes(x=sex, y=diff_sub_use, fill=cue)) + 
  geom_boxplot() + labs(y="Difference in substrate use", x="Sex") + 
  scale_x_discrete(labels=c("F" = "Females", "M" = "Males")) + geom_hline(yintercept=0, linetype='dashed', size = 1) +
  theme(axis.text = element_text(size=12), legend.text = element_text(size=12), axis.title = element_text(size=15)) + 
  stat_compare_means(method="t.test", label = "p.signif", size = 12) + 
  scale_fill_discrete(name = "Cue", labels = c("Alarm","Control"))

#per tank
tank.subusediff.plot <- diff.sub %>% ggplot(aes(x=cue, y=diff_sub_use, fill=tank)) + 
  geom_boxplot() + labs(y="Difference in substrate use", x="Cue") + 
  scale_x_discrete(labels=c("alarm" = "Alarm cue", "control" = "Control cue")) + geom_hline(yintercept=0, linetype='dashed', size = 1) +
  theme(legend.position = "none", axis.text = element_text(size=12), axis.title = element_text(size=15))

#per tank look at mating behaviour
tank.matbehav.plot <- mating.data %>% ggplot(aes(x=cue, y=total_matingProp)) + 
  geom_boxplot() + labs(y="Proportion of time with mating attempts", x="Cue") + 
  scale_x_discrete(labels=c("alarm" = "Alarm cue", "control" = "Control cue")) + geom_hline(yintercept=0, linetype='dashed', size = 1) +
  theme(legend.position = "none", axis.text = element_text(size=12), axis.title = element_text(size=15))

## Statistical Analyses ##
#t-test on change in sub use between sexes within cues to see if there is any reason to include sex in analysis 
#test assumption of equal variance 
var.test(diff_sub_use ~ sex, subset(diff.sub, cue == "alarm"))
var.test(diff_sub_use ~ sex, subset(diff.sub, cue == "control"))

alarm.test.t <- t.test(diff_sub_use ~ sex, subset(diff.sub, cue == "alarm"))
alarm.test.t
control.test.t <- t.test(diff_sub_use ~ sex, subset(diff.sub, cue == "control"))
control.test.t
#no differences from sex 

#run lmm
lmm <- lmer(diff_sub_use ~ cue + sex + (1|tank), diff.sub, REML = TRUE)

#model validation
#homogeneity of variance
par(mar = c(4, 4, 0.5, 0.5))
plot(resid(lmm) ~ fitted(lmm), xlab = "Predicted values", ylab = "Normalized residuals")
abline(h = 0, lty = 2)

#independence of model residuals
boxplot(resid(lmm) ~ diff.sub$cue, xlab = "Cue",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)

boxplot(resid(lmm) ~ diff.sub$sex, xlab = "Sex",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)

#normality of residuals
hist(resid(lmm))

#run signficiance tests
ranova(lmm)

summary(lmm)

car::Anova(lmm, type = 2)

#t-test for difference in mating behaviour
mating.test.t <- t.test(total_matingProp ~ cue, mating.data)
mating.test.t

