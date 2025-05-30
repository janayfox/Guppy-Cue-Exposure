#######################################################
### Goal: Plot and Analyze cue response data
### Author: Janay Fox
### R script
#######################################################

## Set up ##
#install packages
#install.packages("ggplot2")
#install.packages("dplyr")
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
library(forcats)
library(reshape2)
library(tidyverse)
library(lme4)
library(lmerTest)
library(ggpubr)
library(RLRsim)
library(ggpubr)

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

cue.data$sig_prop <- cue.data$sig_s/cue.data$length_obs_s

cue.data$chase_prop <- cue.data$chasing_s/cue.data$length_obs_s

#calculate change in substrate use
diff.sub <- cue.data %>% dcast(fish_ID ~ time, value.var = "subuse_prop", fill=0) %>%
  mutate(diff_sub_use = after - before) %>% select(fish_ID, diff_sub_use)

#calculate averages
mean(filter(cue.data, time == "before" & cue == "alarm")$subuse_prop)
mean(filter(cue.data, time == "before" & cue == "control")$subuse_prop)
mean(filter(cue.data, time == "after" & cue == "alarm")$subuse_prop)
mean(filter(cue.data, time == "after" & cue == "control")$subuse_prop)

#calculate change in mating behaviour
diff.mating <- cue.data %>% dcast(fish_ID ~ time, value.var = "matbehav_prop", fill=0) %>%
  mutate(diff_mating = after - before) %>% select(fish_ID, diff_mating)

diff.sig<- cue.data %>% dcast(fish_ID ~ time, value.var = "sig_prop", fill=0) %>%
  mutate(diff_sig = after - before) %>% select(fish_ID, diff_sig)

diff.chasing <- cue.data %>% dcast(fish_ID ~ time, value.var = "chase_prop", fill=0) %>%
  mutate(diff_chase = after - before) %>% select(fish_ID, diff_chase)

#calculate averages
mean(filter(cue.data, time == "before" & cue == "alarm" & sex == "male")$matbehav_prop)
mean(filter(cue.data, time == "before" & cue == "control" & sex == "male")$matbehav_prop)
mean(filter(cue.data, time == "after" & cue == "alarm" & sex == "male")$matbehav_prop)
mean(filter(cue.data, time == "after" & cue == "control" & sex == "male")$matbehav_prop)

#calculate total mating behaviour and total trial time
total.mating <- cue.data %>% dcast(fish_ID ~ time, value.var = "matbehav_s", fill=0) %>%
  mutate(total_mating = after + before) %>% select(fish_ID, total_mating)

total.trial <- cue.data %>% dcast(fish_ID ~ time, value.var = "length_obs_s", fill=0) %>%
  mutate(total_trial = after + before) %>% select(fish_ID, total_trial)

#merge mating dataframes
mating.data <- merge(diff.mating, total.mating, by = "fish_ID")
mating.data <- merge(mating.data, total.trial, by = "fish_ID")
mating.data <- merge(mating.data, diff.sig, by = "fish_ID")
mating.data <- merge(mating.data, diff.chasing, by = "fish_ID")

#need to add cue and tank back onto new dataframes
diff.sub <- distinct(merge(diff.sub, cue.data[, c("fish_ID", "tank", "cue", "sex")], by = "fish_ID")) 
all.diff.sub <- distinct(merge(all.diff.sub, all.diff.mating, by = c("fish_ID", "bin"))) 
all.diff.sub <- distinct(merge(all.diff.sub, all.diff.sig, by = c("fish_ID", "bin"))) 

mating.data <- distinct(merge(mating.data, cue.data[, c("fish_ID", "tank", "cue", "sex")], by = "fish_ID")) 

#remove females from mating data 
mating.data <- subset(mating.data, sex == "male")

#calculate prop total mating 
mating.data$total_matingProp <- mating.data$total_mating/mating.data$total_trial

## Make plots ##
#boxplot of change in substrate use
#across all tanks
subusediff.plot <- diff.sub %>% ggplot(aes(x=cue, y=diff_sub_use)) + theme_bw() +
  geom_boxplot(fill = c("#E14D2A", "#31C6D4")) + labs(y="Diff. prop. substrate use", x="Cue") +
  scale_x_discrete(labels=c("alarm" = "Alarm cue", "control" = "Control cue")) + geom_hline(yintercept=0, linetype='dashed', size = 1) +
  theme(legend.position = "none", axis.text = element_text(size=12), axis.title = element_text(size=15)) 

subusediff.plot

ggsave("subuseDiff_plot.tiff", plot = subusediff.plot, device = "tiff", width = 4, height = 4, units = "in")

## Statistical Analyses ##
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
#difference in change in proportion
mating.test.t <- t.test(diff_mating ~ cue, mating.data)
mating.test.t

#difference in change in sig
mating.test.t <- t.test(diff_sig ~ cue, mating.data)
mating.test.t

#difference in change in chase
mating.test.t <- t.test(diff_chase ~ cue, mating.data)
mating.test.t

#calculate average proportion of time spent foraging 
cue.data$prop_foraging <- cue.data$feeding_s/cue.data$length_obs_s

mean(filter(cue.data, time == "before" & cue == "alarm")$prop_foraging)
mean(filter(cue.data, time == "before" & cue == "control")$prop_foraging)
mean(filter(cue.data, time == "after" & cue == "alarm")$prop_foraging)
mean(filter(cue.data, time == "after" & cue == "control")$prop_foraging)
