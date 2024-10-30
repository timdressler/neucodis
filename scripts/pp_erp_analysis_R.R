#Bachelor Thesis Analysis Tim Dreßler
#supervised by Kaja Loock, M.Sc. & Prof. Dr. Thomas Bäumer

#Time-dependency of mechanisms involved in episodic memory modulation induced by prediction errors

#Feb23

#TO DO:

#-------------------------------------Set up------------------------------------
#load packages
library(tidyr) 
library(RVAideMemoire) #for Shapiro-Tests
library(readxl)
library(car) 
library(corrplot) 
library(dplyr) 
library(rstudioapi)
library(ez)
library(ggplot2)
library(ggstatsplot)
library(lme4)
library(lmerTest)
library(flexplot)
##library(apaTables) 
library(DescTools)
library(emmeans)
##library(nlme)
library(ggpubr)
library(sjmisc)
library(sjlabelled)
library(sjPlot)
library(sjstats)
library(ggeffects)
library(jtools)
##library(interactions) #very bad performance
##library(reghelper)
##library(rockchalk) #does not work for GLMMs
##library(arsenal) #for comparing df's, not used
library(rstatix)
##library(PairedData) #for paired plots
library(cowplot) #for merging plots
library(tidyverse)
library(devtools)
library(smplot2)
library(glmmTMB)
##library(lsr)
library(rmcorr)
library(psych)
library(simr) #for power analysis
library(MoMAColors)

#clean WS
rm(list=ls())

options(scipen = 999)

#set WS
setwd("C:/Users/timdr/OneDrive/Dokumente/_Microsoft Office & Dokumente/_Dateien/HFT Stuttgart/6. Semester/Bachelor Thesis/Data")