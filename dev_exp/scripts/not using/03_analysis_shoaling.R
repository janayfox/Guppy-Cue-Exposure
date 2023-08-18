############################################
### Goal: Plot and analyze shoaling data
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
data <- read_csv(here("gup_cue_exp", "data", "clean", "clean_shoaling.csv"))

#convert character columns to factors
data[sapply(data, is.character)] <- lapply(data[sapply(data, is.character)], as.factor)

#calculate time differences for shoal vs empty container
data$int_diff_s <- data$int_shoal_s - data$int_emp_s
data$int_lat_diff_s <- data$int_emp_lat_s - data$int_shoal_lat_s 
data$tight_diff_s <- data$tight_shoal_s - data$tight_emp_s
data$tight_lat_diff_s <- data$tight_emp_lat_s - data$tight_shoal_lat_s
data$lse_diff_s <- data$lse_shoal_s - data$lse_emp_s
data$lse_lat_diff_s <- data$lse_emp_lat_s - data$lse_shoal_lat_s 

#calculate percentage of time shoaling
data$per_int <- data$int_shoal_s / 300 * 100
data$per_tight <- data$tight_shoal_s / 300 * 100
data$per_lse <- data$lse_shoal_s / 300 * 100


## Exploratory plots ##
#Time loose shoaling
lse.plot <-data %>% ggplot(aes(x=cue, y=lse_shoal_s, fill=sex)) + 
  geom_boxplot() + labs(y="Time spent loose shoaling (s)", x="Cue") + 
  scale_x_discrete(labels=c("ac" = "Alarm cue", "c" = "Control cue")) + 
  theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=15)) + 
  scale_fill_discrete(labels = c("Females", "Males"))
lse.plot

#Percentage time loose shoaling
per.lse.plot <-data %>% ggplot(aes(x=cue, y=per_lse, fill=sex)) + 
  geom_boxplot() + labs(y="Time spent loose shoaling (s)", x="Cue") + 
  scale_x_discrete(labels=c("ac" = "Alarm cue", "c" = "Control cue")) + 
  theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=15)) + 
  scale_fill_discrete(labels = c("Females", "Males"))
per.lse.plot

#Latency to loose shoal
lat.lse.plot <-data %>% ggplot(aes(x=cue, y=lse_shoal_lat_s, fill=sex)) + 
  geom_boxplot() + labs(y="Latency to loose shoal (s)", x="Cue") + 
  scale_x_discrete(labels=c("ac" = "Alarm cue", "c" = "Control cue")) + 
  theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=15)) + 
  scale_fill_discrete(labels = c("Females", "Males"))
lat.lse.plot

#Time tight shoaling 
tight.plot <-data %>% ggplot(aes(x=cue, y=tight_shoal_s, fill=sex)) + 
  geom_boxplot() + labs(y="Time spent tight shoaling (s)", x="Cue") + 
  scale_x_discrete(labels=c("ac" = "Alarm cue", "c" = "Control cue")) + 
  theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=15)) + 
  scale_fill_discrete(labels = c("Females", "Males"))
tight.plot

#Percentage time tight shoaling
per.tight.plot <-data %>% ggplot(aes(x=cue, y=per_tight, fill=sex)) + 
  geom_boxplot() + labs(y="Time spent loose shoaling (s)", x="Cue") + 
  scale_x_discrete(labels=c("ac" = "Alarm cue", "c" = "Control cue")) + 
  theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=15)) + 
  scale_fill_discrete(labels = c("Females", "Males"))
per.tight.plot

#Latency to tight shoal
lat.tight.plot <-data %>% ggplot(aes(x=cue, y=tight_shoal_lat_s, fill=sex)) + 
  geom_boxplot() + labs(y="Latency to tight shoal (s)", x="Cue") + 
  scale_x_discrete(labels=c("ac" = "Alarm cue", "c" = "Control cue")) + 
  theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=15)) + 
  scale_fill_discrete(labels = c("Females", "Males"))
lat.tight.plot

#Time spent interacting 
int.plot <-data %>% ggplot(aes(x=cue, y=int_shoal_s, fill=sex)) + 
  geom_boxplot() + labs(y="Time spent interacting (s)", x="Cue") + 
  scale_x_discrete(labels=c("ac" = "Alarm cue", "c" = "Control cue")) + 
  theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=15)) + 
  scale_fill_discrete(labels = c("Females", "Males"))
int.plot

#Percentage time interacting
per.int.plot <-data %>% ggplot(aes(x=cue, y=per_int, fill=sex)) + 
  geom_boxplot() + labs(y="Time spent loose shoaling (s)", x="Cue") + 
  scale_x_discrete(labels=c("ac" = "Alarm cue", "c" = "Control cue")) + 
  theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=15)) + 
  scale_fill_discrete(labels = c("Females", "Males"))
per.int.plot

#Latency to interact 
lat.int.plot <-data %>% ggplot(aes(x=cue, y=int_shoal_lat_s, fill=sex)) + 
  geom_boxplot() + labs(y="Latency to interact (s)", x="Cue") + 
  scale_x_discrete(labels=c("ac" = "Alarm cue", "c" = "Control cue")) + 
  theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=15)) + 
  scale_fill_discrete(labels = c("Females", "Males"))
lat.int.plot

#Time difference loose 
diff.lse.plot <-data %>% ggplot(aes(x=cue, y=lse_diff_s, fill=sex)) + 
  geom_boxplot() + labs(y="Time difference loose shoal (s)", x="Cue") + 
  scale_x_discrete(labels=c("ac" = "Alarm cue", "c" = "Control cue")) + 
  theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=15)) + 
  scale_fill_discrete(labels = c("Females", "Males"))
diff.lse.plot

#Time difference loose latency
diff.lse.lat.plot <-data %>% ggplot(aes(x=cue, y=lse_lat_diff_s, fill=sex)) + 
  geom_boxplot() + labs(y="Time difference latency to loose shoal (s)", x="Cue") + 
  scale_x_discrete(labels=c("ac" = "Alarm cue", "c" = "Control cue")) + 
  theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=15)) + 
  scale_fill_discrete(labels = c("Females", "Males"))
diff.lse.lat.plot

#Time difference tight 
diff.tight.plot <-data %>% ggplot(aes(x=cue, y=tight_diff_s, fill=sex)) + 
  geom_boxplot() + labs(y="Time difference tight shoal (s)", x="Cue") + 
  scale_x_discrete(labels=c("ac" = "Alarm cue", "c" = "Control cue")) + 
  theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=15)) + 
  scale_fill_discrete(labels = c("Females", "Males"))
diff.tight.plot

#Time difference tight latency
diff.tight.lat.plot <-data %>% ggplot(aes(x=cue, y=tight_lat_diff_s, fill=sex)) + 
  geom_boxplot() + labs(y="Time difference latency to tight shoal (s)", x="Cue") + 
  scale_x_discrete(labels=c("ac" = "Alarm cue", "c" = "Control cue")) + 
  theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=15)) + 
  scale_fill_discrete(labels = c("Females", "Males"))
diff.tight.lat.plot

#Time difference interacting 
diff.int.plot <-data %>% ggplot(aes(x=cue, y=int_diff_s, fill=sex)) + 
  geom_boxplot() + labs(y="Time difference interacting (s)", x="Cue") + 
  scale_x_discrete(labels=c("ac" = "Alarm cue", "c" = "Control cue")) + 
  theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=15)) + 
  scale_fill_discrete(labels = c("Females", "Males"))
diff.int.plot

#Time difference interact latency
diff.int.lat.plot <-data %>% ggplot(aes(x=cue, y=int_lat_diff_s, fill=sex)) + 
  geom_boxplot() + labs(y="Time difference latency to interact (s)", x="Cue") + 
  scale_x_discrete(labels=c("ac" = "Alarm cue", "c" = "Control cue")) + 
  theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=15)) + 
  scale_fill_discrete(labels = c("Females", "Males"))
diff.int.lat.plot

## Run models ##
#Check for correlations between measurements
res.diff <- cor(data[,c(32,34,36)])
round(res.diff, 2)

res.per <- cor(data[,c(38,39,40)])
round(res.per, 2)

#check for need to include tank as a random effect
lm.test <- lm(int_diff_s ~ cue * sex, data = data)
lm.test.resid <- rstandard(lm.test)
plot(lm.test.resid ~ data$tank, xlab = "Tank", ylab = "Standardized residuals")
#residuals vary so should be included 

#run LMMs
int.lmm.sc <- lmer(int_diff_s ~ cue * sex + (1|tank), data = data)
int.lmm.c <- lmer(int_diff_s ~ cue + (1|tank), data = data)
int.lmm.s <- lmer(int_diff_s ~ sex + (1|tank), data = data)

int.lat.lmm.sc <- lmer(int_lat_diff_s ~ cue * sex + (1|tank), data = data)
int.lat.lmm.c <- lmer(int_lat_diff_s ~ cue + (1|tank), data = data)
int.lat.lmm.s <- lmer(int_lat_diff_s ~ sex + (1|tank), data = data)

tight.lmm.sc <- lmer(tight_diff_s ~ cue * sex + (1|tank), data = data)
tight.lmm.c <- lmer(tight_diff_s ~ cue + (1|tank), data = data)
tight.lmm.s <- lmer(tight_diff_s ~ sex + (1|tank), data = data)

tight.lat.lmm.sc <- lmer(tight_lat_diff_s ~ cue * sex + (1|tank), data = data)
tight.lat.lmm.c <- lmer(tight_lat_diff_s ~ cue + (1|tank), data = data)
tight.lat.lmm.s <- lmer(tight_lat_diff_s ~ sex + (1|tank), data = data)

lse.lmm.sc <- lmer(lse_diff_s ~ cue * sex + (1|tank), data = data)
lse.lmm.c <- lmer(lse_diff_s ~ cue + (1|tank), data = data)
lse.lmm.s <- lmer(lse_diff_s ~ sex + (1|tank), data = data)

lse.lat.lmm.sc <- lmer(lse_lat_diff_s ~ cue * sex + (1|tank), data = data)
lse.lat.lmm.c <- lmer(lse_lat_diff_s ~ cue + (1|tank), data = data)
lse.lat.lmm.s <- lmer(lse_lat_diff_s ~ sex + (1|tank), data = data)

#compare models
AIC.table.int <- model.sel(int.lmm.sc, int.lmm.c, int.lmm.s)
(AIC.table.int <- AIC.table.int[, c("df", "logLik", "AICc", "delta")])

AIC.table.lat.int <- model.sel(int.lat.lmm.sc, int.lat.lmm.c, int.lat.lmm.s)
(AIC.table.lat.int <- AIC.table.lat.int[, c("df", "logLik", "AICc", "delta")])

AIC.table.tight <- model.sel(tight.lmm.sc, tight.lmm.c, tight.lmm.s)
(AIC.table.tight <- AIC.table.tight[, c("df", "logLik", "AICc", "delta")])

AIC.table.lat.tight <- model.sel(tight.lat.lmm.sc, tight.lat.lmm.c, tight.lat.lmm.s)
(AIC.table.lat.tight <- AIC.table.lat.tight[, c("df", "logLik", "AICc", "delta")])

AIC.table.lse <- model.sel(lse.lmm.sc, lse.lmm.c, lse.lmm.s)
(AIC.table.lse <- AIC.table.lse[, c("df", "logLik", "AICc", "delta")])

AIC.table.lat.lse <- model.sel(lse.lat.lmm.sc, lse.lat.lmm.c, lse.lat.lmm.s)
(AIC.table.lat.lse <- AIC.table.lat.lse[, c("df", "logLik", "AICc", "delta")])

#model validation
#homogeneity of variance
par(mar = c(4, 4, 0.5, 0.5))
plot(resid(int.lmm.sc) ~ fitted(int.lmm.sc), xlab = "Predicted values", ylab = "Normalized residuals")
abline(h = 0, lty = 2)

par(mar = c(4, 4, 0.5, 0.5))
plot(resid(int.lat.lmm.sc) ~ fitted(int.lat.lmm.sc), xlab = "Predicted values", ylab = "Normalized residuals")
abline(h = 0, lty = 2)

par(mar = c(4, 4, 0.5, 0.5))
plot(resid(tight.lmm.sc) ~ fitted(tight.lmm.sc), xlab = "Predicted values", ylab = "Normalized residuals")
abline(h = 0, lty = 2)

par(mar = c(4, 4, 0.5, 0.5))
plot(resid(tight.lat.lmm.sc) ~ fitted(tight.lat.lmm.sc), xlab = "Predicted values", ylab = "Normalized residuals")
abline(h = 0, lty = 2)

par(mar = c(4, 4, 0.5, 0.5))
plot(resid(lse.lmm.sc) ~ fitted(lse.lmm.sc), xlab = "Predicted values", ylab = "Normalized residuals")
abline(h = 0, lty = 2)

par(mar = c(4, 4, 0.5, 0.5))
plot(resid(lse.lat.lmm.sc) ~ fitted(lse.lat.lmm.sc), xlab = "Predicted values", ylab = "Normalized residuals")
abline(h = 0, lty = 2)

#independence of model residuals
boxplot(resid(int.lmm.sc) ~ data$cue, xlab = "Cue",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)
boxplot(resid(int.lmm.sc) ~ data$sex, xlab = "Sex",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)

boxplot(resid(int.lat.lmm.sc) ~ data$cue, xlab = "Cue",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)
boxplot(resid(int.lat.lmm.sc) ~ data$sex, xlab = "Sex",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)

boxplot(resid(tight.lmm.sc) ~ data$cue, xlab = "Cue",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)
boxplot(resid(tight.lmm.sc) ~ data$sex, xlab = "Sex",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)

boxplot(resid(tight.lat.lmm.sc) ~ data$cue, xlab = "Cue",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)
boxplot(resid(tight.lat.lmm.sc) ~ data$sex, xlab = "Sex",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)

boxplot(resid(lse.lmm.sc) ~ data$cue, xlab = "Cue",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)
boxplot(resid(lse.lmm.sc) ~ data$sex, xlab = "Sex",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)

boxplot(resid(lse.lat.lmm.sc) ~ data$cue, xlab = "Cue",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)
boxplot(resid(lse.lat.lmm.sc) ~ data$sex, xlab = "Sex",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)

#normality of residuals
hist(resid(int.lmm.sc))
hist(resid(int.lat.lmm.sc))

hist(resid(tight.lmm.sc))
hist(resid(tight.lat.lmm.sc))

hist(resid(lse.lmm.sc))
hist(resid(lse.lat.lmm.sc))

#residual diagnostics
simulationOutput <- simulateResiduals(fittedModel = int.lmm.sc)
plot(simulationOutput)

simulationOutput <- simulateResiduals(fittedModel = int.lat.lmm.sc)
plot(simulationOutput)

simulationOutput <- simulateResiduals(fittedModel = lse.lmm.sc)
plot(simulationOutput)

simulationOutput <- simulateResiduals(fittedModel = lse.lat.lmm.sc)
plot(simulationOutput)

simulationOutput <- simulateResiduals(fittedModel = tight.lmm.sc)
plot(simulationOutput)

simulationOutput <- simulateResiduals(fittedModel = tight.lat.lmm.sc)
plot(simulationOutput)

#interpret results 
summary(int.lmm.sc)
summary(int.lat.lmm.sc)

summary(tight.lmm.sc)
summary(tight.lat.lmm.sc)

summary(lse.lmm.sc)
summary(lse.lat.lmm.sc)

## Export plots and metadata ##
ggsave(here("gup_cue_exp", "plots", "shoaling", "sh_diff_lse_plot.jpg"), plot = diff.lse.plot)
ggsave(here("gup_cue_exp", "plots", "shoaling", "sh_diff_tight_plot.jpg"), plot = diff.tight.plot)
ggsave(here("gup_cue_exp", "plots", "shoaling", "sh_diff_int_plot.jpg"), plot = diff.int.plot)

ggsave(here("gup_cue_exp", "plots", "shoaling", "sh_diff_lat_lse_plot.jpg"), plot = diff.lse.lat.plot)
ggsave(here("gup_cue_exp", "plots", "shoaling", "sh_diff_lat_tight_plot.jpg"), plot = diff.tight.lat.plot)
ggsave(here("gup_cue_exp", "plots", "shoaling", "sh_diff_lat_int_plot.jpg"), plot = diff.int.lat.plot)

write_csv(data, file(here("gup_cue_exp", "data", "metadata", "sh_metadata.csv")))

