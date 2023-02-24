###########
### Goal: Figure out sample numbers submitted
### Author: Janay Fox
###########


sub.data <- read_csv(here("data", "raw", "sub.csv"))

sub.data$study <- str_sub(sub.data$ID, 1, 1)
sub.data$sex <- str_sub(sub.data$ID, -2)
sub.data$sex <- gsub('[0-9]*', '', sub.data$sex)

dev.data <- subset(sub.data, study == "D")
dev.data$sex[dev.data$ID == "DAC5M1-52"] <- "M"

st.data <- subset(sub.data, study == "S")

mal.dev.data <- subset(dev.data, sex == "M")
fem.dev.data <- subset(dev.data, sex == "F")
