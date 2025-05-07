library(ggplot2)
library(dplyr)
setwd("C:/Users/annik/OneDrive/Desktop/BU Fall 2024/MA675/Genetics Project")

HKSM <- read.csv("EG_HKSM_20E_col_03_11.csv")

TwentyE <- read.csv("EG_20E_col_03_11.csv")

control <- read.csv("EG_Con_col_03_11.csv")


control <- control %>% mutate(treatment = "H2O")
HKSM <- HKSM %>% mutate(treatment = "HKSM")
TwentyE <- TwentyE %>% mutate(treatment = "20E")
cdata <- bind_rows(control, HKSM, TwentyE)
cdata <- cdata[, -32]
cdata <- cdata %>%
  mutate(treatment = factor(treatment, levels = c("h20", "HKSM", "20E")))


############################################################
# Some EDA Looking at the HKSM data

ggplot(cdata, aes(x=new_act_score)) + geom_histogram(color="black", fill="white")+
  ggtitle("Histogram of activity score for all treatments")

### Histograms By Each Treatment
ggplot(HKSM, aes(x=new_act_score)) + geom_histogram(color="black", fill="white")+
  ggtitle("Histogram of activity score for hksm")


ggplot(TwentyE, aes(x=new_act_score)) + geom_histogram(color="black", fill="white")+
  ggtitle("Histogram of activity score for 20e")


ggplot(control, aes(x=new_act_score)) + geom_histogram(color="black", fill="white")+
  ggtitle("Histogram of activity score for control")























