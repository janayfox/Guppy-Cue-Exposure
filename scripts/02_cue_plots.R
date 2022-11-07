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
fivemin.cue.data <- read_csv(here("gup_cue_exp", "data", "clean", "clean_cue_fivemin.csv"))
twomin.cue.data <- read_csv(here("gup_cue_exp", "data", "clean", "clean_cue_twomin.csv"))

#convert character columns to factors
fivemin.cue.data[sapply(fivemin.cue.data, is.character)] <- lapply(fivemin.cue.data[sapply(fivemin.cue.data, is.character)], as.factor)
twomin.cue.data[sapply(twomin.cue.data, is.character)] <- lapply(twomin.cue.data[sapply(twomin.cue.data, is.character)], as.factor)

#order factors
fivemin.cue.data$time <- factor(fivemin.cue.data$time, levels=c('before', 'after'))
twomin.cue.data$time <- factor(twomin.cue.data$time, levels=c('before', 'after'))

#remove tank that doesnt have after data 
fivemin.cue.data <- subset(fivemin.cue.data, tank != "ST2AC13")
twomin.cue.data <- subset(twomin.cue.data, tank != "ST2AC13")

## Calculate behavioural measurements ##
#calculate time fish spent at bottom of tank (in Q5 or Q6)
fivemin.cue.data$bottom_s <- (fivemin.cue.data$Q5_s + fivemin.cue.data$Q6_s)
twomin.cue.data$bottom_s <- (twomin.cue.data$Q5_s + twomin.cue.data$Q6_s)

#calculate substate use (time at bottom of tank minus time spent foraging)
fivemin.cue.data$subuse_s <- (fivemin.cue.data$bottom_s - fivemin.cue.data$feeding_s)
twomin.cue.data$subuse_s <- (twomin.cue.data$bottom_s - twomin.cue.data$feeding_s)

#convert substrate use to a proportion
fivemin.cue.data$subuse_prop <- fivemin.cue.data$subuse_s/fivemin.cue.data$length_obs_s
twomin.cue.data$subuse_prop <- twomin.cue.data$subuse_s/twomin.cue.data$length_obs_s

#convert time frozen to proportion
fivemin.cue.data$frozen_prop <- fivemin.cue.data$frozen_s/fivemin.cue.data$length_obs_s
twomin.cue.data$frozen_prop <- twomin.cue.data$frozen_s/twomin.cue.data$length_obs_s

#combine mating data 
fivemin.cue.data$matbehav_s <- fivemin.cue.data$sig_s + fivemin.cue.data$chasing_s
twomin.cue.data$matbehav_s <- twomin.cue.data$sig_s + twomin.cue.data$chasing_s

#calculate change in substrate use
fivemin.diff.sub <- fivemin.cue.data %>% dcast(fish_ID ~ time, value.var = "subuse_prop", fill=0) %>%
  mutate(diff_sub_use = after - before) %>% select(fish_ID, diff_sub_use)

twomin.diff.sub <- twomin.cue.data %>% dcast(fish_ID ~ time, value.var = "subuse_prop", fill=0) %>%
  mutate(diff_sub_use = after - before) %>% select(fish_ID, diff_sub_use)

#need to add cue and back onto diff subuse dataframe 
fivemin.diff.sub <- distinct(merge(fivemin.diff.sub, fivemin.cue.data[, c("fish_ID", "tank", "cue", "sex")], by = "fish_ID")) 
twomin.diff.sub <- distinct(merge(twomin.diff.sub, fivemin.cue.data[, c("fish_ID", "tank", "cue", "sex")], by = "fish_ID")) 

## Make plots ##
#boxplot of proportion of substrate use
##if I want to use these plots I will have to find a way to put lines between points overtop of the boxplots
fivemin.prop.plot <- fivemin.cue.data %>% ggplot(aes(x=cue, y=subuse_prop, fill=time)) + 
  geom_boxplot() + labs(y="Proportion of trial spent near substrate", x="Cue") + 
  scale_x_discrete(labels=c("alarm" = "Alarm cue", "control" = "Control cue")) + 
  theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=15)) + 
  scale_fill_discrete(labels = c("Before cue", "After cue")) 

twomin.prop.plot <- twomin.cue.data %>% ggplot(aes(x=cue, y=subuse_prop, fill=time)) + 
  geom_boxplot() + labs(y="Proportion of trial spent near substrate", x="Cue") + 
  scale_x_discrete(labels=c("alarm" = "Alarm cue", "control" = "Control cue")) + 
  theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=15)) + 
  scale_fill_discrete(labels = c("Before cue", "After cue")) 

#make sex specific
fem.fivemin.prop.plot <- subset(fivemin.cue.data, sex == "female") %>% ggplot(aes(x=cue, y=subuse_prop, fill=time)) + 
  geom_boxplot() + labs(y="Proportion of trial spent near substrate", x="Cue") + 
  scale_x_discrete(labels=c("alarm" = "Alarm cue", "control" = "Control cue")) + 
  theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=15)) + 
  scale_fill_discrete(labels = c("Before cue", "After cue")) 

mal.fivemin.prop.plot <- subset(fivemin.cue.data, sex == "male") %>% ggplot(aes(x=cue, y=subuse_prop, fill=time)) + 
  geom_boxplot() + labs(y="Proportion of trial spent near substrate", x="Cue") + 
  scale_x_discrete(labels=c("alarm" = "Alarm cue", "control" = "Control cue")) + 
  theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=15)) + 
  scale_fill_discrete(labels = c("Before cue", "After cue")) 

fem.twomin.prop.plot <- subset(twomin.cue.data, sex == "female") %>% ggplot(aes(x=cue, y=subuse_prop, fill=time)) + 
  geom_boxplot() + labs(y="Proportion of trial spent near substrate", x="Cue") + 
  scale_x_discrete(labels=c("alarm" = "Alarm cue", "control" = "Control cue")) + 
  theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=15)) + 
  scale_fill_discrete(labels = c("Before cue", "After cue")) 

mal.twomin.prop.plot <- subset(twomin.cue.data, sex == "male") %>% ggplot(aes(x=cue, y=subuse_prop, fill=time)) + 
  geom_boxplot() + labs(y="Proportion of trial spent near substrate", x="Cue") + 
  scale_x_discrete(labels=c("alarm" = "Alarm cue", "control" = "Control cue")) + 
  theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=15)) + 
  scale_fill_discrete(labels = c("Before cue", "After cue")) 

#boxplot of change in substrate use
#across all tanks
fivemin.subusediff.plot <- fivemin.diff.sub %>% ggplot(aes(x=cue, y=diff_sub_use)) + 
  geom_boxplot() + labs(y="Difference in substrate use", x="Cue") + 
  scale_x_discrete(labels=c("alarm" = "Alarm cue", "control" = "Control cue")) + geom_hline(yintercept=0, linetype='dashed', size = 1) +
  theme(legend.position = "none", axis.text = element_text(size=12), axis.title = element_text(size=15)) 

twomin.subusediff.plot <- twomin.diff.sub %>% ggplot(aes(x=cue, y=diff_sub_use)) + 
  geom_boxplot() + labs(y="Difference in substrate use", x="Cue") + 
  scale_x_discrete(labels=c("alarm" = "Alarm cue", "control" = "Control cue")) + geom_hline(yintercept=0, linetype='dashed', size = 1) +
  theme(legend.position = "none", axis.text = element_text(size=12), axis.title = element_text(size=15))

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
per.fivemin.subusediff.plot <- fivemin.diff.sub %>% ggplot(aes(x=cue, y=diff_sub_use, fill=tank)) + 
  geom_boxplot() + labs(y="Difference in substrate use", x="Cue") + 
  scale_x_discrete(labels=c("alarm" = "Alarm cue", "control" = "Control cue")) + geom_hline(yintercept=0, linetype='dashed', size = 1) +
  theme(legend.position = "none", axis.text = element_text(size=12), axis.title = element_text(size=15))

per.twomin.subusediff.plot <- twomin.diff.sub %>% ggplot(aes(x=cue, y=diff_sub_use, fill=tank)) + 
   geom_boxplot() + labs(y="Difference in substrate use", x="Cue") + 
   scale_x_discrete(labels=c("alarm" = "Alarm cue", "control" = "Control cue")) + geom_hline(yintercept=0, linetype='dashed', size = 1) +
   theme(legend.position = "none", axis.text = element_text(size=12), axis.title = element_text(size=15))

#per tank look at mating behaviour
per.fivemin.matbehav.plot <- fivemin.cue.data %>% ggplot(aes(x=cue, y=matbehav_s, fill=tank)) + 
  geom_boxplot() + labs(y="Time male spent engaging in mating behaviour (s)", x="Cue") + 
  scale_x_discrete(labels=c("alarm" = "Alarm cue", "control" = "Control cue")) + geom_hline(yintercept=0, linetype='dashed', size = 1) +
  theme(legend.position = "none", axis.text = element_text(size=12), axis.title = element_text(size=15))

## Export plots and data ##
ggsave(here("gup_cue_exp", "plots", "cue_exp", "fivemin.prop.plot.jpg"), plot = fivemin.prop.plot)
ggsave(here("gup_cue_exp", "plots", "cue_exp", "twomin.prop.plot.jpg"), plot = twomin.prop.plot)
ggsave(here("gup_cue_exp", "plots", "cue_exp", "fem.fivemin.prop.plot.jpg"), plot = fem.fivemin.prop.plot)
ggsave(here("gup_cue_exp", "plots", "cue_exp", "mal.fivemin.prop.plot.jpg"), plot = mal.fivemin.prop.plot)
ggsave(here("gup_cue_exp", "plots", "cue_exp", "fem.twomin.prop.plot.jpg"), plot = fem.twomin.prop.plot)
ggsave(here("gup_cue_exp", "plots", "cue_exp", "mal.twomin.prop.plot.jpg"), plot = fem.twomin.prop.plot)

ggsave(here("gup_cue_exp", "plots", "cue_exp", "all.tank.fivemin.subusediff.plot.jpg"), plot = fivemin.subusediff.plot)
ggsave(here("gup_cue_exp", "plots", "cue_exp", "all.tank.twomin.subusediff.plot.jpg"), plot = twomin.subusediff.plot)
ggsave(here("gup_cue_exp", "plots", "cue_exp", "per.tank.fivemin.subusediff.plot.jpg"), plot = per.fivemin.subusediff.plot)
ggsave(here("gup_cue_exp", "plots", "cue_exp", "per.tank.twomin.subusediff.plot.jpg"), plot = per.twomin.subusediff.plot)

ggsave(here("gup_cue_exp", "plots", "cue_exp", "pval.tank.fivemin.subusediff.plot.jpg"), plot = pval.fivemin.subusediff.plot)
ggsave(here("gup_cue_exp", "plots", "cue_exp", "pval.tank.twomin.subusediff.plot.jpg"), plot = pval.twomin.subusediff.plot)

write_csv(fivemin.diff.sub, file(here("gup_cue_exp", "data", "metadata", "cue_diffsub_fivemin_metadata.csv")))
write_csv(twomin.diff.sub, file(here("gup_cue_exp", "data", "metadata", "cue_diffsub_twomin_metadata.csv")))

write_csv(fivemin.cue.data, file(here("gup_cue_exp", "data", "metadata", "cue_fivemin_metadata.csv")))
write_csv(twomin.cue.data, file(here("gup_cue_exp", "data", "metadata", "cue_twomin_metadata.csv")))

## Statistical Analyses ##
#t-test on change in sub use between sexes within cues to see if there is any reason to include sex in analysis 
#test assumption of equal variance 
var.test(diff_sub_use ~ sex, subset(fivemin.diff.sub, cue == "alarm"))
var.test(diff_sub_use ~ sex, subset(fivemin.diff.sub, cue == "control"))
var.test(diff_sub_use ~ sex, subset(twomin.diff.sub, cue == "alarm"))
var.test(diff_sub_use ~ sex, subset(twomin.diff.sub, cue == "control"))

alarm.test.t <- t.test(diff_sub_use ~ sex, subset(fivemin.diff.sub, cue == "alarm"))
alarm.test.t
control.test.t <- t.test(diff_sub_use ~ sex, subset(fivemin.diff.sub, cue == "control"))
control.test.t
alarm.test.t <- t.test(diff_sub_use ~ sex, subset(twomin.diff.sub, cue == "alarm"))
alarm.test.t
control.test.t <- t.test(diff_sub_use ~ sex, subset(twomin.diff.sub, cue == "control"))
control.test.t
#no differences from sex so do not include in further analysis

#do t.test within sex to test between treatments 
var.test(diff_sub_use ~ cue, subset(fivemin.diff.sub, sex == "female"))
fivemin.test.t <- t.test(diff_sub_use ~ cue, subset(fivemin.diff.sub, sex == "female"))
fivemin.test.t

var.test(diff_sub_use ~ cue, subset(twomin.diff.sub, sex == "female"))
twomin.test.t <- t.test(diff_sub_use ~ cue, subset(twomin.diff.sub, sex == "female"))
twomin.test.t

var.test(diff_sub_use ~ cue, subset(fivemin.diff.sub, sex == "male"))
fivemin.test.t <- t.test(diff_sub_use ~ cue, subset(fivemin.diff.sub, sex == "male"))
fivemin.test.t

var.test(diff_sub_use ~ cue, subset(twomin.diff.sub, sex == "male"))
twomin.test.t <- t.test(diff_sub_use ~ cue, subset(twomin.diff.sub, sex == "male"))
twomin.test.t

#run lmms
five.lmm <- lmer(diff_sub_use ~ cue + (1|tank), fivemin.diff.sub)
two.lmm <- lmer(diff_sub_use ~ cue + (1|tank), twomin.diff.sub)

#run simulation to determine if variance between tanks is significant
exactRLRT(five.lmm)
exactRLRT(two.lmm)

summary(five.lmm)
summary(two.lmm)

#check assumptions
par(mar = c(4, 4, 0.5, 0.5))
#homogenous dispersion of residuals
plot(resid(five.lmm) ~ fitted(five.lmm), xlab = "Predicted values", ylab = "Normalized residuals")
abline(h = 0, lty = 2)
#independence of the model residuals with each co-variate
par(mfrow = c(1, 2), mar = c(4, 4, 0.5, 0.5))
boxplot(resid(five.lmm) ~ cue, data = fivemin.diff.sub, xlab = "Cue",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)
boxplot(resid(five.lmm) ~ tank, data = fivemin.diff.sub, xlab = "Tank", ylab = "Normalized residuals")
abline(h = 0, lty = 2)
#normality of residuals
hist(resid(five.lmm))

