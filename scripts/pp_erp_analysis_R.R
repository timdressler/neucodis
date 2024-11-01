#Practical Project Tim Dreßler
#supervised by Prof. Dr. Herrmann & Prof. Dr. Hildebrandt

#Neuronal Markers of Corollary Discharge

#Oct24

#-------------------------------------Set up------------------------------------
#load packages
library(tidyr) 
library(RVAideMemoire) 
library(readxl)
library(car) 
library(corrplot) 
library(dplyr) 
library(rstudioapi)
library(ez)
library(ggplot2)
library(ggstatsplot)
library(DescTools)
library(ggpubr)
library(sjmisc)
library(sjlabelled)
library(sjPlot)
library(sjstats)
library(ggeffects)
library(jtools)
library(rstatix)
library(cowplot) 
library(tidyverse)
library(devtools)
library(smplot2)
library(glmmTMB)
library(rmcorr)
library(psych)
library(MoMAColors)

#clean WS
rm(list=ls())

options(scipen = 999)

#set costum colors
colors <- list()
colors$main_blue <- "#004F9F"
colors$main_red <- "#D53D0E"
colors$main_green <- "#00786B"
colors$light_blue <- "#5BC5F2"
colors$main_yellow <- "#FDC300"
colors$UI <- "grey"

#setup paths
SCRIPTPATH <- dirname(rstudioapi::getSourceEditorContext()$path)
if (grepl("neucodis/scripts", SCRIPTPATH)) {
  cat("Path OK\n") 
} else {
  stop("Path not OK")  
}
MAINPATH <- gsub("neucodis/scripts", "", SCRIPTPATH)
INPATH <- file.path(MAINPATH, "data", "analysis_data")
OUTPATH <- file.path(MAINPATH, "data", "analysis_data")

setwd(INPATH)

#---------------------------------N100 Analysis---------------------------------

#load data
df_erp <- read.csv(file.path(INPATH, "erp_analysis.csv"))

#ERP1.0
#T Test (N100 Amplitude ~ Condition)
#check if condition (talk/no talk) has an effect on the N100 amplitude
t.test(data = df_erp, erp_amp ~ cond, paired = TRUE) #
df_erp %>% 
  cohens_d(erp_amp ~ cond, paired = TRUE) 

psych::describeBy(df_erp$erp_amp,
                  group = df_erp$cond)

#PLOT: N100 amplitude by condition
P1 <- df_erp %>%
  ggplot(aes(x = cond, y = erp_amp, fill = cond)) +
  sm_raincloud(boxplot.params = list(fill = "white", outlier.shape = NA), violin.params = list(alpha = 1)
               , point.params = list(), legends = F) +
  scale_fill_manual(values = c(colors$main_blue, colors$main_red)) +
  scale_x_discrete(labels = c('talk' = 'talk', 'listen' = 'listen')) +
  scale_y_continuous(n.breaks = 15) +
  labs(x = NULL, y = "Amplitude") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 12)) +
  theme(legend.position="none") +
  geom_signif(comparisons=list(c("talk", "listen")), annotations="WATCHOUT",
              y_position = 1, tip_length = 0.02,  vjust=0.4) 
P1

#ERP1.1
#T Test (N100 Latency ~ Condition)
#check if condition (talk/no talk) has an effect on the N100 amplitude
t.test(data = df_erp, erp_lat ~ cond, paired = TRUE) #
df_erp %>% 
  cohens_d(erp_lat ~ cond, paired = TRUE) 

psych::describeBy(df_erp$erp_lat,
                  group = df_erp$cond)

#PLOT: N100 latency by condition
P2 <- df_erp %>%
  ggplot(aes(x = cond, y = erp_lat, fill = cond)) +
  sm_raincloud(boxplot.params = list(fill = "white", outlier.shape = NA), violin.params = list(alpha = 1)
               , point.params = list(), legends = F) +
  scale_fill_manual(values = c(colors$main_blue, colors$main_red)) +
  scale_x_discrete(labels = c('talk' = 'talk', 'listen' = 'listen')) +
  scale_y_continuous(n.breaks = 10, limits = c(50,150)) +
  labs(x = NULL, y = "Latency [ms]") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 12)) +
  theme(legend.position="none") +
  geom_signif(comparisons=list(c("talk", "listen")), annotations="WATCHOUT",
              y_position = 130, tip_length = 0.02,  vjust=0.4) 
P2

#--------------------------------Sample_analysis--------------------------------








