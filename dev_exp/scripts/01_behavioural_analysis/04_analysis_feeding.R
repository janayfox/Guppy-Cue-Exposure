############################################
### Goal: Plot and analyze feeding data
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
data <- read_csv(here("gup_cue_exp", "data", "clean", "clean_feeding.csv"))

#convert character columns to factors
data[sapply(data, is.character)] <- lapply(data[sapply(data, is.character)], as.factor)

## Exploratory plots ##
#number of pecks
peck.plot <-data %>% ggplot(aes(x=cue, y=peck_nb, fill=sex)) + 
  geom_boxplot() + labs(y="# of pecks", x="Cue") + 
  scale_x_discrete(labels=c("ac" = "Alarm cue", "c" = "Control cue")) + 
  theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=15)) + 
  scale_fill_discrete(labels = c("Females", "Males"))
peck.plot

#time within
within.plot <-data %>% ggplot(aes(x=cue, y=within_s, fill=sex)) + 
  geom_boxplot() + labs(y="Time within feeder (s)", x="Cue") + 
  scale_x_discrete(labels=c("ac" = "Alarm cue", "c" = "Control cue")) + 
  theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=15)) + 
  scale_fill_discrete(labels = c("Females", "Males"))
within.plot

#time near
near.plot <-data %>% ggplot(aes(x=cue, y=near_s, fill=sex)) + 
  geom_boxplot() + labs(y="Time near feeder (s)", x="Cue") + 
  scale_x_discrete(labels=c("ac" = "Alarm cue", "c" = "Control cue")) + 
  theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=15)) + 
  scale_fill_discrete(labels = c("Females", "Males"))
near.plot

## Run models ##
#check for need to include tank as a random effect
lm.test <- lm(within_s ~ cue * sex, data = data)
lm.test.resid <- rstandard(lm.test)
plot(lm.test.resid ~ data$tank, xlab = "Tank", ylab = "Standardized residuals")
#residuals vary so should be included 

#run LMMs
in.lmm.sc <- lmer(within_s ~ cue * sex + (1|tank), data = data)
in.lmm.c <- lmer(within_s ~ cue + (1|tank), data = data)
in.lmm.s <- lmer(within_s ~ sex + (1|tank), data = data)

#compare models
AIC.table.in <- model.sel(in.lmm.sc, in.lmm.c, in.lmm.s)
(AIC.table.in <- AIC.table.in[, c("df", "logLik", "AICc", "delta")])

#model validation
#homogeneity of variance
par(mar = c(4, 4, 0.5, 0.5))
plot(resid(in.lmm.sc) ~ fitted(in.lmm.sc), xlab = "Predicted values", ylab = "Normalized residuals")
abline(h = 0, lty = 2)

#independence of model residuals
boxplot(resid(in.lmm.sc) ~ data$cue, xlab = "Cue",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)

boxplot(resid(in.lmm.sc) ~ data$sex, xlab = "Sex",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)

#normality of residuals
hist(resid(in.lmm.sc))

#residual diagnostics
simulationOutput <- simulateResiduals(fittedModel = in.lmm.sc)
plot(simulationOutput)

#interpret results
summary(in.lmm.sc)

## Export plots and metadata ##
ggsave(here("gup_cue_exp", "plots", "feeding", "feeding_within_plot.jpg"), plot = within.plot)
ggsave(here("gup_cue_exp", "plots", "feeding", "feeding_peck_plot.jpg"), plot = peck.plot)

write_csv(data, file(here("gup_cue_exp", "data", "metadata", "feed_metadata.csv")))

