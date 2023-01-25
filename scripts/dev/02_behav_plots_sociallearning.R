###########
### Goal: Plot and analyze feeding data
### Author: Janay Fox
###########

### Set up ###
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

#read data files
feeding.data <- read_csv(here("gup_cue_exp", "data", "clean", "clean_feeding.csv"))
full.data <- read_csv(here("gup_cue_exp", "data", "clean", "clean_full_sl.csv"))
half.data <- read_csv(here("gup_cue_exp", "data", "clean", "clean_5min_sl.csv"))

#make list of fish that fed during feeding trial
peck_ID <- subset(feeding.data, peck_nb > 0, select=c(ID)) 

#subset only fish in peck list
peck.full.data <- full.data[full.data$ID %in% peck_ID$ID,]
peck.half.data <- half.data[half.data$ID %in% peck_ID$ID,]

#remove males from peck subset since very few fed 
peck.full.data <- subset(peck.full.data, sex == "f")
peck.half.data <- subset(peck.half.data, sex == "f")

#convert character columns to factors
full.data[sapply(full.data, is.character)] <- lapply(full.data[sapply(full.data, is.character)], as.factor)
half.data[sapply(half.data, is.character)] <- lapply(half.data[sapply(half.data, is.character)], as.factor)
peck.full.data[sapply(peck.full.data, is.character)] <- lapply(peck.full.data[sapply(peck.full.data, is.character)], as.factor)
peck.half.data[sapply(peck.half.data, is.character)] <- lapply(peck.half.data[sapply(peck.half.data, is.character)], as.factor)

#calculate preference 
full.data$in_diff_s <- full.data$in_dem_s - full.data$in_nondem_s
full.data$nr_diff_s <- full.data$nr_dem_s - full.data$nr_nondem_s
full.data$far_diff_s <- full.data$far_dem_s - full.data$far_nondem_s
full.data$in_lat_diff_s <- full.data$in_dem_lat_s - full.data$in_nondem_lat_s

half.data$in_diff_s <- half.data$in_dem_s - half.data$in_nondem_s
half.data$nr_diff_s <- half.data$nr_dem_s - half.data$nr_nondem_s
half.data$far_diff_s <- half.data$far_dem_s - half.data$far_nondem_s
half.data$in_lat_diff_s <- half.data$in_dem_lat_s - half.data$in_nondem_lat_s

peck.full.data$in_diff_s <- peck.full.data$in_dem_s - peck.full.data$in_nondem_s
peck.full.data$nr_diff_s <- peck.full.data$nr_dem_s - peck.full.data$nr_nondem_s
peck.full.data$far_diff_s <- peck.full.data$far_dem_s - peck.full.data$far_nondem_s
peck.full.data$in_lat_diff_s <- peck.full.data$in_dem_lat_s - peck.full.data$in_nondem_lat_s

peck.half.data$in_diff_s <- peck.half.data$in_dem_s - peck.half.data$in_nondem_s
peck.half.data$nr_diff_s <- peck.half.data$nr_dem_s - peck.half.data$nr_nondem_s
peck.half.data$far_diff_s <- peck.half.data$far_dem_s - peck.half.data$far_nondem_s
peck.half.data$in_lat_diff_s <- peck.half.data$in_dem_lat_s - peck.half.data$in_nondem_lat_s

### Exploratory plots ###
#pref for inside demonstrated feeder
full.in.pref.plot <-full.data %>% ggplot(aes(x=cue, y=in_diff_s, fill=sex)) + 
  geom_boxplot() + labs(y="Difference in time inside feeder (s)", x="Cue") + 
  scale_x_discrete(labels=c("ac" = "Alarm cue", "c" = "Control cue")) + 
  theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=15)) + 
  scale_fill_discrete(labels = c("Females", "Males"))
full.in.pref.plot

half.in.pref.plot <-half.data %>% ggplot(aes(x=cue, y=in_diff_s, fill=sex)) + 
  geom_boxplot() + labs(y="Difference in time inside feeder (s)", x="Cue") + 
  scale_x_discrete(labels=c("ac" = "Alarm cue", "c" = "Control cue")) + 
  theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=15)) + 
  scale_fill_discrete(labels = c("Females", "Males"))
half.in.pref.plot

peck.full.in.pref.plot <-peck.full.data %>% ggplot(aes(x=cue, y=in_diff_s, fill=sex)) + 
  geom_boxplot() + labs(y="Difference in time inside feeder (s)", x="Cue") + 
  scale_x_discrete(labels=c("ac" = "Alarm cue", "c" = "Control cue")) + 
  theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=15)) + 
  scale_fill_discrete(labels = c("Females", "Males"))
peck.full.in.pref.plot

peck.half.in.pref.plot <-peck.half.data %>% ggplot(aes(x=cue, y=in_diff_s, fill=sex)) + 
  geom_boxplot() + labs(y="Difference in time inside feeder (s)", x="Cue") + 
  scale_x_discrete(labels=c("ac" = "Alarm cue", "c" = "Control cue")) + 
  theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=15)) + 
  scale_fill_discrete(labels = c("Females", "Males"))
peck.half.in.pref.plot

#pref for near demonstrated feeder
full.nr.pref.plot <-full.data %>% ggplot(aes(x=cue, y=nr_diff_s, fill=sex)) + 
  geom_boxplot() + labs(y="Difference in time near feeder (s)", x="Cue") + 
  scale_x_discrete(labels=c("ac" = "Alarm cue", "c" = "Control cue")) + 
  theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=15)) + 
  scale_fill_discrete(labels = c("Females", "Males"))
full.nr.pref.plot

half.nr.pref.plot <-half.data %>% ggplot(aes(x=cue, y=nr_diff_s, fill=sex)) + 
  geom_boxplot() + labs(y="Difference in time near feeder (s)", x="Cue") + 
  scale_x_discrete(labels=c("ac" = "Alarm cue", "c" = "Control cue")) + 
  theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=15)) + 
  scale_fill_discrete(labels = c("Females", "Males"))
half.nr.pref.plot

peck.full.nr.pref.plot <-peck.full.data %>% ggplot(aes(x=cue, y=nr_diff_s, fill=sex)) + 
  geom_boxplot() + labs(y="Difference in time near feeder (s)", x="Cue") + 
  scale_x_discrete(labels=c("ac" = "Alarm cue", "c" = "Control cue")) + 
  theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=15)) + 
  scale_fill_discrete(labels = c("Females", "Males"))
peck.full.nr.pref.plot

peck.half.nr.pref.plot <-peck.half.data %>% ggplot(aes(x=cue, y=nr_diff_s, fill=sex)) + 
  geom_boxplot() + labs(y="Difference in time near feeder (s)", x="Cue") + 
  scale_x_discrete(labels=c("ac" = "Alarm cue", "c" = "Control cue")) + 
  theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=15)) + 
  scale_fill_discrete(labels = c("Females", "Males"))
peck.half.nr.pref.plot

#pref for far demonstrated feeder
full.far.pref.plot <-full.data %>% ggplot(aes(x=cue, y=far_diff_s, fill=sex)) + 
  geom_boxplot() + labs(y="Difference in time near feeder (s)", x="Cue") + 
  scale_x_discrete(labels=c("ac" = "Alarm cue", "c" = "Control cue")) + 
  theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=15)) + 
  scale_fill_discrete(labels = c("Females", "Males"))
full.far.pref.plot

half.far.pref.plot <-half.data %>% ggplot(aes(x=cue, y=far_diff_s, fill=sex)) + 
  geom_boxplot() + labs(y="Difference in time near feeder (s)", x="Cue") + 
  scale_x_discrete(labels=c("ac" = "Alarm cue", "c" = "Control cue")) + 
  theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=15)) + 
  scale_fill_discrete(labels = c("Females", "Males"))
half.far.pref.plot

peck.full.far.pref.plot <-peck.full.data %>% ggplot(aes(x=cue, y=far_diff_s, fill=sex)) + 
  geom_boxplot() + labs(y="Difference in time near feeder (s)", x="Cue") + 
  scale_x_discrete(labels=c("ac" = "Alarm cue", "c" = "Control cue")) + 
  theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=15)) + 
  scale_fill_discrete(labels = c("Females", "Males"))
peck.full.far.pref.plot

peck.half.far.pref.plot <-peck.half.data %>% ggplot(aes(x=cue, y=far_diff_s, fill=sex)) + 
  geom_boxplot() + labs(y="Difference in time near feeder (s)", x="Cue") + 
  scale_x_discrete(labels=c("ac" = "Alarm cue", "c" = "Control cue")) + 
  theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=15)) + 
  scale_fill_discrete(labels = c("Females", "Males"))
peck.half.far.pref.plot

#diff in latency to enter feeder
full.in.lat.diff.plot <-full.data %>% ggplot(aes(x=cue, y=in_lat_diff_s, fill=sex)) + 
  geom_boxplot() + labs(y="Difference in latency to enter feeder (s)", x="Cue") + 
  scale_x_discrete(labels=c("ac" = "Alarm cue", "c" = "Control cue")) + 
  theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=15)) + 
  scale_fill_discrete(labels = c("Females", "Males"))
full.in.lat.diff.plot

half.in.lat.diff.plot <-half.data %>% ggplot(aes(x=cue, y=in_lat_diff_s, fill=sex)) + 
  geom_boxplot() + labs(y="Difference in latency to enter feeder (s)", x="Cue") + 
  scale_x_discrete(labels=c("ac" = "Alarm cue", "c" = "Control cue")) + 
  theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=15)) + 
  scale_fill_discrete(labels = c("Females", "Males"))
half.in.lat.diff.plot

peck.full.in.lat.diff.plot <-peck.full.data %>% ggplot(aes(x=cue, y=in_lat_diff_s, fill=sex)) + 
  geom_boxplot() + labs(y="Difference in latency to enter feeder (s)", x="Cue") + 
  scale_x_discrete(labels=c("ac" = "Alarm cue", "c" = "Control cue")) + 
  theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=15)) + 
  scale_fill_discrete(labels = c("Females", "Males"))
peck.full.in.lat.diff.plot

peck.half.in.lat.diff.plot <-peck.half.data %>% ggplot(aes(x=cue, y=in_lat_diff_s, fill=sex)) + 
  geom_boxplot() + labs(y="Difference in latency to enter feeder (s)", x="Cue") + 
  scale_x_discrete(labels=c("ac" = "Alarm cue", "c" = "Control cue")) + 
  theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=15)) + 
  scale_fill_discrete(labels = c("Females", "Males"))
peck.half.in.lat.diff.plot

### Run models ###
#check for need to include tank as a random effect
lm.test <- lm(in_diff_s ~ cue * sex, data = full.data)
lm.test.resid <- rstandard(lm.test)
plot(lm.test.resid ~ full.data$tank, xlab = "Tank", ylab = "Standardized residuals")

lm.test <- lm(in_diff_s ~ cue * sex, data = half.data)
lm.test.resid <- rstandard(lm.test)
plot(lm.test.resid ~ half.data$tank, xlab = "Tank", ylab = "Standardized residuals")

lm.test <- lm(in_diff_s ~ cue, data = peck.full.data)
lm.test.resid <- rstandard(lm.test)
plot(lm.test.resid ~ peck.full.data$tank, xlab = "Tank", ylab = "Standardized residuals")

lm.test <- lm(in_diff_s ~ cue, data = peck.half.data)
lm.test.resid <- rstandard(lm.test)
plot(lm.test.resid ~ peck.half.data$tank, xlab = "Tank", ylab = "Standardized residuals")
#residuals vary so should be included 

#run LMMs
full.in.diff.lmm.sc <- lmer(in_diff_s ~ cue * sex + (1|tank), data = full.data)
full.in.diff.lmm.c <- lmer(in_diff_s ~ cue + (1|tank), data = full.data)
full.in.diff.lmm.s <- lmer(in_diff_s ~ sex + (1|tank), data = full.data)

full.nr.diff.lmm.sc <- lmer(nr_diff_s ~ cue * sex + (1|tank), data = full.data)
full.nr.diff.lmm.c <- lmer(nr_diff_s ~ cue + (1|tank), data = full.data)
full.nr.diff.lmm.s <- lmer(nr_diff_s ~ sex + (1|tank), data = full.data)

full.far.diff.lmm.sc <- lmer(far_diff_s ~ cue * sex + (1|tank), data = full.data)
full.far.diff.lmm.c <- lmer(far_diff_s ~ cue + (1|tank), data = full.data)
full.far.diff.lmm.s <- lmer(far_diff_s ~ sex + (1|tank), data = full.data)

full.in.lat.diff.lmm.sc <- lmer(in_lat_diff_s ~ cue * sex + (1|tank), data = full.data)
full.in.lat.diff.lmm.c <- lmer(in_lat_diff_s ~ cue + (1|tank), data = full.data)
full.in.lat.diff.lmm.s <- lmer(in_lat_diff_s ~ sex + (1|tank), data = full.data)

half.in.diff.lmm.sc <- lmer(in_diff_s ~ cue * sex + (1|tank), data = half.data)
half.in.diff.lmm.c <- lmer(in_diff_s ~ cue + (1|tank), data = half.data)
half.in.diff.lmm.s <- lmer(in_diff_s ~ sex + (1|tank), data = half.data)

half.nr.diff.lmm.sc <- lmer(nr_diff_s ~ cue * sex + (1|tank), data = half.data)
half.nr.diff.lmm.c <- lmer(nr_diff_s ~ cue + (1|tank), data = half.data)
half.nr.diff.lmm.s <- lmer(nr_diff_s ~ sex + (1|tank), data = half.data)

half.far.diff.lmm.sc <- lmer(far_diff_s ~ cue * sex + (1|tank), data = half.data)
half.far.diff.lmm.c <- lmer(far_diff_s ~ cue + (1|tank), data = half.data)
half.far.diff.lmm.s <- lmer(far_diff_s ~ sex + (1|tank), data = half.data)

half.in.lat.diff.lmm.sc <- lmer(in_lat_diff_s ~ cue * sex + (1|tank), data = half.data)
half.in.lat.diff.lmm.c <- lmer(in_lat_diff_s ~ cue + (1|tank), data = half.data)
half.in.lat.diff.lmm.s <- lmer(in_lat_diff_s ~ sex + (1|tank), data = half.data)

peck.full.in.diff.lmm.c <- lmer(in_diff_s ~ cue + (1|tank), data = peck.full.data)
peck.full.nr.diff.lmm.c <- lmer(nr_diff_s ~ cue + (1|tank), data = peck.full.data)
peck.full.far.diff.lmm.c <- lmer(far_diff_s ~ cue + (1|tank), data = peck.full.data)
peck.full.in.lat.diff.lmm.c <- lmer(in_lat_diff_s ~ cue + (1|tank), data = peck.full.data)

peck.half.in.diff.lmm.c <- lmer(in_diff_s ~ cue + (1|tank), data = peck.half.data)
peck.half.nr.diff.lmm.c <- lmer(nr_diff_s ~ cue + (1|tank), data = peck.half.data)
peck.half.far.diff.lmm.c <- lmer(far_diff_s ~ cue + (1|tank), data = peck.half.data)
peck.half.in.lat.diff.lmm.c <- lmer(in_lat_diff_s ~ cue + (1|tank), data = peck.half.data)

#compare models
AIC.table.in.full <- model.sel(full.in.diff.lmm.sc, full.in.diff.lmm.c, full.in.diff.lmm.s)
(AIC.table.in.full <- AIC.table.in.full[, c("df", "logLik", "AICc", "delta")])
AIC.table.nr.full <- model.sel(full.nr.diff.lmm.sc, full.nr.diff.lmm.c, full.nr.diff.lmm.s)
(AIC.table.nr.full <- AIC.table.nr.full[, c("df", "logLik", "AICc", "delta")])
AIC.table.far.full <- model.sel(full.far.diff.lmm.sc, full.far.diff.lmm.c, full.far.diff.lmm.s)
(AIC.table.far.full <- AIC.table.far.full[, c("df", "logLik", "AICc", "delta")])
AIC.table.in.lat.full <- model.sel(full.in.lat.diff.lmm.sc, full.in.lat.diff.lmm.c, full.in.lat.diff.lmm.s)
(AIC.table.in.lat.full <- AIC.table.in.lat.full[, c("df", "logLik", "AICc", "delta")])

AIC.table.in.half <- model.sel(half.in.diff.lmm.sc, half.in.diff.lmm.c, half.in.diff.lmm.s)
(AIC.table.in.half <- AIC.table.in.half[, c("df", "logLik", "AICc", "delta")])
AIC.table.nr.half <- model.sel(half.nr.diff.lmm.sc, half.nr.diff.lmm.c, half.nr.diff.lmm.s)
(AIC.table.nr.half <- AIC.table.nr.half[, c("df", "logLik", "AICc", "delta")])
AIC.table.far.half <- model.sel(half.far.diff.lmm.sc, half.far.diff.lmm.c, half.far.diff.lmm.s)
(AIC.table.far.half <- AIC.table.far.half[, c("df", "logLik", "AICc", "delta")])
AIC.table.in.lat.half <- model.sel(half.in.lat.diff.lmm.sc, half.in.lat.diff.lmm.c, half.in.lat.diff.lmm.s)
(AIC.table.in.lat.half <- AIC.table.in.lat.half[, c("df", "logLik", "AICc", "delta")])

#model validation
#homogeneity of variance
par(mar = c(4, 4, 0.5, 0.5))
plot(resid(full.in.diff.lmm.sc) ~ fitted(full.in.diff.lmm.sc), xlab = "Predicted values", ylab = "Normalized residuals")
abline(h = 0, lty = 2)

par(mar = c(4, 4, 0.5, 0.5))
plot(resid(full.nr.diff.lmm.sc) ~ fitted(full.nr.diff.lmm.sc), xlab = "Predicted values", ylab = "Normalized residuals")
abline(h = 0, lty = 2)

par(mar = c(4, 4, 0.5, 0.5))
plot(resid(full.far.diff.lmm.sc) ~ fitted(full.far.diff.lmm.sc), xlab = "Predicted values", ylab = "Normalized residuals")
abline(h = 0, lty = 2)

par(mar = c(4, 4, 0.5, 0.5))
plot(resid(full.in.lat.diff.lmm.sc) ~ fitted(full.in.lat.diff.lmm.sc), xlab = "Predicted values", ylab = "Normalized residuals")
abline(h = 0, lty = 2)

par(mar = c(4, 4, 0.5, 0.5))
plot(resid(half.in.diff.lmm.sc) ~ fitted(half.in.diff.lmm.sc), xlab = "Predicted values", ylab = "Normalized residuals")
abline(h = 0, lty = 2)

par(mar = c(4, 4, 0.5, 0.5))
plot(resid(half.nr.diff.lmm.sc) ~ fitted(half.nr.diff.lmm.sc), xlab = "Predicted values", ylab = "Normalized residuals")
abline(h = 0, lty = 2)

par(mar = c(4, 4, 0.5, 0.5))
plot(resid(half.far.diff.lmm.sc) ~ fitted(half.far.diff.lmm.sc), xlab = "Predicted values", ylab = "Normalized residuals")
abline(h = 0, lty = 2)

par(mar = c(4, 4, 0.5, 0.5))
plot(resid(half.in.lat.diff.lmm.sc) ~ fitted(half.in.lat.diff.lmm.sc), xlab = "Predicted values", ylab = "Normalized residuals")
abline(h = 0, lty = 2)

par(mar = c(4, 4, 0.5, 0.5))
plot(resid(peck.full.in.diff.lmm.c) ~ fitted(peck.full.in.diff.lmm.c), xlab = "Predicted values", ylab = "Normalized residuals")
abline(h = 0, lty = 2)

par(mar = c(4, 4, 0.5, 0.5))
plot(resid(peck.full.nr.diff.lmm.c) ~ fitted(peck.full.nr.diff.lmm.c), xlab = "Predicted values", ylab = "Normalized residuals")
abline(h = 0, lty = 2)

par(mar = c(4, 4, 0.5, 0.5))
plot(resid(peck.full.far.diff.lmm.c) ~ fitted(peck.full.far.diff.lmm.c), xlab = "Predicted values", ylab = "Normalized residuals")
abline(h = 0, lty = 2)

par(mar = c(4, 4, 0.5, 0.5))
plot(resid(peck.full.in.lat.diff.lmm.c) ~ fitted(peck.full.in.lat.diff.lmm.c), xlab = "Predicted values", ylab = "Normalized residuals")
abline(h = 0, lty = 2)

par(mar = c(4, 4, 0.5, 0.5))
plot(resid(peck.half.in.diff.lmm.c) ~ fitted(peck.half.in.diff.lmm.c), xlab = "Predicted values", ylab = "Normalized residuals")
abline(h = 0, lty = 2)

par(mar = c(4, 4, 0.5, 0.5))
plot(resid(peck.half.nr.diff.lmm.c) ~ fitted(peck.half.nr.diff.lmm.c), xlab = "Predicted values", ylab = "Normalized residuals")
abline(h = 0, lty = 2)

par(mar = c(4, 4, 0.5, 0.5))
plot(resid(peck.half.far.diff.lmm.c) ~ fitted(peck.half.far.diff.lmm.c), xlab = "Predicted values", ylab = "Normalized residuals")
abline(h = 0, lty = 2)

par(mar = c(4, 4, 0.5, 0.5))
plot(resid(peck.half.in.lat.diff.lmm.c) ~ fitted(peck.half.in.lat.diff.lmm.c), xlab = "Predicted values", ylab = "Normalized residuals")
abline(h = 0, lty = 2)

#independence of model residuals
boxplot(resid(full.in.diff.lmm.sc) ~ full.data$cue, xlab = "Cue",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)
boxplot(resid(full.in.diff.lmm.sc) ~ full.data$sex, xlab = "Sex",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)

boxplot(resid(full.nr.diff.lmm.sc) ~ full.data$cue, xlab = "Cue",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)
boxplot(resid(full.nr.diff.lmm.sc) ~ full.data$sex, xlab = "Sex",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)

boxplot(resid(full.far.diff.lmm.sc) ~ full.data$cue, xlab = "Cue",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)
boxplot(resid(full.far.diff.lmm.sc) ~ full.data$sex, xlab = "Sex",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)

boxplot(resid(full.in.lat.diff.lmm.sc) ~ full.data$cue, xlab = "Cue",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)
boxplot(resid(full.in.lat.diff.lmm.sc) ~ full.data$sex, xlab = "Sex",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)

boxplot(resid(half.in.diff.lmm.sc) ~ half.data$cue, xlab = "Cue",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)
boxplot(resid(half.in.diff.lmm.sc) ~ half.data$sex, xlab = "Sex",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)

boxplot(resid(half.nr.diff.lmm.sc) ~ half.data$cue, xlab = "Cue",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)
boxplot(resid(half.nr.diff.lmm.sc) ~ half.data$sex, xlab = "Sex",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)

boxplot(resid(half.far.diff.lmm.sc) ~ half.data$cue, xlab = "Cue",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)
boxplot(resid(half.far.diff.lmm.sc) ~ half.data$sex, xlab = "Sex",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)

boxplot(resid(half.in.lat.diff.lmm.sc) ~ half.data$cue, xlab = "Cue",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)
boxplot(resid(half.in.lat.diff.lmm.sc) ~ half.data$sex, xlab = "Sex",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)

boxplot(resid(peck.full.in.diff.lmm.c) ~ peck.full.data$cue, xlab = "Cue",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)

boxplot(resid(peck.full.nr.diff.lmm.c) ~ peck.full.data$cue, xlab = "Cue",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)

boxplot(resid(peck.full.far.diff.lmm.c) ~ peck.full.data$cue, xlab = "Cue",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)

boxplot(resid(peck.full.in.lat.diff.lmm.c) ~ peck.full.data$cue, xlab = "Cue",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)

boxplot(resid(peck.half.in.diff.lmm.c) ~ peck.full.data$cue, xlab = "Cue",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)

boxplot(resid(peck.half.nr.diff.lmm.c) ~ peck.half.data$cue, xlab = "Cue",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)

boxplot(resid(peck.half.far.diff.lmm.c) ~ peck.half.data$cue, xlab = "Cue",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)

boxplot(resid(peck.half.in.lat.diff.lmm.c) ~ peck.half.data$cue, xlab = "Cue",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)

#normality of residuals
hist(resid(full.in.diff.lmm.sc))
hist(resid(full.nr.diff.lmm.sc))
hist(resid(full.far.diff.lmm.sc))
hist(resid(full.in.lat.diff.lmm.sc))

hist(resid(half.in.diff.lmm.sc))
hist(resid(half.nr.diff.lmm.sc))
hist(resid(half.far.diff.lmm.sc))
hist(resid(half.in.lat.diff.lmm.sc))

hist(resid(peck.full.in.diff.lmm.c))
hist(resid(peck.full.nr.diff.lmm.c))
hist(resid(peck.full.far.diff.lmm.c))
hist(resid(peck.full.in.lat.diff.lmm.c))

hist(resid(peck.half.in.diff.lmm.c))
hist(resid(peck.half.nr.diff.lmm.c))
hist(resid(peck.half.far.diff.lmm.c))
hist(resid(peck.half.in.lat.diff.lmm.c))

#residual diagnostics
simulationOutput <- simulateResiduals(fittedModel = full.in.diff.lmm.sc)
plot(simulationOutput)
simulationOutput <- simulateResiduals(fittedModel = full.nr.diff.lmm.sc)
plot(simulationOutput)
simulationOutput <- simulateResiduals(fittedModel = full.far.diff.lmm.sc)
plot(simulationOutput)
simulationOutput <- simulateResiduals(fittedModel = full.in.lat.diff.lmm.sc)
plot(simulationOutput)

simulationOutput <- simulateResiduals(fittedModel = half.in.diff.lmm.sc)
plot(simulationOutput)
simulationOutput <- simulateResiduals(fittedModel = half.nr.diff.lmm.sc)
plot(simulationOutput)
simulationOutput <- simulateResiduals(fittedModel = half.far.diff.lmm.sc)
plot(simulationOutput)
simulationOutput <- simulateResiduals(fittedModel = half.in.lat.diff.lmm.sc)
plot(simulationOutput)

simulationOutput <- simulateResiduals(fittedModel = peck.full.in.diff.lmm.c)
plot(simulationOutput)
simulationOutput <- simulateResiduals(fittedModel = peck.full.nr.diff.lmm.c)
plot(simulationOutput)
simulationOutput <- simulateResiduals(fittedModel = peck.full.far.diff.lmm.c)
plot(simulationOutput)
simulationOutput <- simulateResiduals(fittedModel = peck.full.in.lat.diff.lmm.c)
plot(simulationOutput)

simulationOutput <- simulateResiduals(fittedModel = peck.half.in.diff.lmm.c)
plot(simulationOutput)
simulationOutput <- simulateResiduals(fittedModel = peck.half.nr.diff.lmm.c)
plot(simulationOutput)
simulationOutput <- simulateResiduals(fittedModel = peck.half.far.diff.lmm.c)
plot(simulationOutput)
simulationOutput <- simulateResiduals(fittedModel = peck.half.in.lat.diff.lmm.c)
plot(simulationOutput)

#interpret results
summary(full.in.diff.lmm.sc)
summary(full.nr.diff.lmm.sc)
summary(full.far.diff.lmm.sc)
summary(full.in.lat.diff.lmm.sc)

summary(half.in.diff.lmm.sc)
summary(half.nr.diff.lmm.sc)
summary(half.far.diff.lmm.sc)
summary(half.in.lat.diff.lmm.sc)

summary(peck.full.in.diff.lmm.c)
summary(peck.full.nr.diff.lmm.c)
summary(peck.full.far.diff.lmm.c)
summary(peck.full.in.lat.diff.lmm.c)

summary(peck.half.in.diff.lmm.c)
summary(peck.half.nr.diff.lmm.c)
summary(peck.half.far.diff.lmm.c)
summary(peck.half.in.lat.diff.lmm.c)

### Check for preference of left or right feeder ### 
LR.data <- full.data[,c(1,2,4,22:24)] #subset just columns needed
LR.data <- gather(LR.data, side, in_s, in_L_s:in_R_s, factor_key=TRUE) #convert to long format
LR.glmm <- lmer(in_s ~ side * sex * cue + (1|tank) + (1|ID), data = LR.data)
summary(LR.glmm)

full.peck.LR.data <- peck.full.data[,c(1,2,4,22:24)] #subset just columns needed
full.peck.LR.data <- gather(full.peck.LR.data, side, in_s, in_L_s:in_R_s, factor_key=TRUE) #convert to long format
LR.glmm <- lmer(in_s ~ side * cue + (1|tank) + (1|ID), data = full.peck.LR.data)
summary(LR.glmm)

full.peck.LR.data <- peck.full.data[,c(1,2,4,22:24)] #subset just columns needed
full.peck.LR.data <- gather(full.peck.LR.data, side, in_s, in_L_s:in_R_s, factor_key=TRUE) #convert to long format
LR.glmm <- lmer(in_s ~ side *  cue + (1|tank) + (1|ID), data = full.peck.LR.data)
summary(LR.glmm)

### Export plots and metadata ####
ggsave(here("gup_cue_exp", "plots", "social_learning", "sl_full_diff_in_plot.jpg"), plot = full.in.pref.plot)
ggsave(here("gup_cue_exp", "plots", "social_learning", "sl_full_diff_nr_plot.jpg"), plot = full.nr.pref.plot)
ggsave(here("gup_cue_exp", "plots", "social_learning", "sl_full_diff_far_plot.jpg"), plot = full.far.pref.plot)
ggsave(here("gup_cue_exp", "plots", "social_learning", "sl_full_diff_in_lat_plot.jpg"), plot = full.in.lat.diff.plot)

ggsave(here("gup_cue_exp", "plots", "social_learning", "sl_half_diff_in_plot.jpg"), plot = half.in.pref.plot)
ggsave(here("gup_cue_exp", "plots", "social_learning", "sl_half_diff_nr_plot.jpg"), plot = half.nr.pref.plot)
ggsave(here("gup_cue_exp", "plots", "social_learning", "sl_half_diff_far_plot.jpg"), plot = half.far.pref.plot)
ggsave(here("gup_cue_exp", "plots", "social_learning", "sl_half_diff_in_lat_plot.jpg"), plot = half.in.lat.diff.plot)

ggsave(here("gup_cue_exp", "plots", "social_learning", "sl_peck_full_diff_in_plot.jpg"), plot = peck.full.in.pref.plot)
ggsave(here("gup_cue_exp", "plots", "social_learning", "sl_peck_full_diff_nr_plot.jpg"), plot = peck.full.nr.pref.plot)
ggsave(here("gup_cue_exp", "plots", "social_learning", "sl_peck_full_diff_far_plot.jpg"), plot = peck.full.far.pref.plot)
ggsave(here("gup_cue_exp", "plots", "social_learning", "sl_peck_full_diff_in_lat_plot.jpg"), plot = peck.full.in.lat.diff.plot)

ggsave(here("gup_cue_exp", "plots", "social_learning", "sl_peck_half_diff_in_plot.jpg"), plot = peck.half.in.pref.plot)
ggsave(here("gup_cue_exp", "plots", "social_learning", "sl_peck_half_diff_nr_plot.jpg"), plot = peck.half.nr.pref.plot)
ggsave(here("gup_cue_exp", "plots", "social_learning", "sl_peck_half_diff_far_plot.jpg"), plot = peck.half.far.pref.plot)
ggsave(here("gup_cue_exp", "plots", "social_learning", "sl_peck_half_diff_in_lat_plot.jpg"), plot = peck.half.in.lat.diff.plot)

write_csv(peck.half.data, file(here("gup_cue_exp", "data", "metadata", "sl_halfpeck_metadata.csv")))
write_csv(peck.full.data, file(here("gup_cue_exp", "data", "metadata", "sl_fullpeck_metadata.csv")))
