###R script creating box plots for TZA duration of immunity paper
#Naomi and Alex
#October 2025

#Set up: Libaries; working directory; theme:
### scatter plot ####
setwd ("C:/Users/NaomiMiall/Documents/Polio work/Tanzania")

## library
library(haven)
library(ggplot2)
library(ggpubr)
library(ggtext)
library(scales)
library(dplyr)
library(tidyverse)
library(tidyr)
library(readr)
library(tibble)
library(lme4)
library(lmerTest)
library(emmeans)
library(pbkrtest)

## mytheme
mytheme2 <- theme_bw() + 
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.border = element_blank(),
        legend.background = element_rect(fill=NA),
        legend.position = "bottom")

longitudinal <- read_dta("Both_datasets_long.dta")
#vis <- subset (longitudinal, visit=="bl"|visit=="w2"|visit=="w4"|visit=="w6"|visit=="w8")
#vis_ser <- subset (longitudinal, vis=="bl"|vis=="w4"|vis=="w8")

longitudinal$visit <- factor(longitudinal$visit, levels = c("1", "2", "3", "4", "5", "6", "7", "8"),
                             labels = c("Baseline", "One month", "3 months", "6 months", "1 year", "2 years", "3 years", "8-9 years"))

#Explore age outliers
agepv1_visit <- ggplot (longitudinal, aes(x=factor(visit), y=logPV1_stool))+
  geom_jitter(aes(color = factor(flag_ageoutlier)),alpha=0.6, position=position_jitter(width = 0.2, height = 0), size=2)+
  scale_color_manual(name = "Age outlier", values = c("0" = "grey", "1" = "red"),
                     labels = c("0" = "0", "1" = "1")) +
  mytheme2+
  scale_y_continuous(breaks = c(2, 4, 6, 8, 10), limits=c(1,10)) +
  labs( x = "visit", y = ("PV1 log<sub>2</sub> neutralization titer"))+
  scale_x_discrete(labels= c("Baseline", "1 month", "3 months", "6 months", "1 year", "2 years", "3 years", "8-9 years"))+
  theme(plot.caption = element_text(hjust = 0.5, size = 14),legend.position = "none", axis.text.x = element_text(size = 18, angle=45, hjust=0.75), axis.text.y = element_text(size = 14), axis.title = element_markdown(size = 18), axis.title.x = element_blank())
ggsave("Figures/Age_outliers_PV1.png", width = 6, height = 6, dpi = 300)
agepv1_visit

agepv1 <- ggplot (longitudinal, aes(x=age_bopv14w, y=logPV1_stool))+
  geom_jitter(aes(color = factor(flag_ageoutlier)),alpha=0.6, position=position_jitter(width = 0.2, height = 0), size=2)+
  scale_color_manual(name = "Age outlier", values = c("0" = "grey", "1" = "red"),
                     labels = c("0" = "0", "1" = "1")) +
  mytheme2+
  scale_y_continuous(breaks = c(2, 4, 6, 8, 10), limits=c(1,10)) +
  labs( x = "Age at vaccination", y = ("PV1 log<sub>2</sub> neutralization titer"))+
   theme(plot.caption = element_text(hjust = 0.5, size = 14),legend.position = "none", axis.text.x = element_text(size = 18, angle=45, hjust=0.75), axis.text.y = element_text(size = 14), axis.title = element_markdown(size = 18), axis.title.x = element_blank())
#ggsave("Figures/Age_outliers_PV1.png", width = 6, height = 6, dpi = 300)
agepv1

agepv3_visit <- ggplot (longitudinal, aes(x=factor(visit), y=logPV3_stool))+
  geom_jitter(aes(color = factor(flag_ageoutlier)),alpha=0.6, position=position_jitter(width = 0.2, height = 0), size=2)+
  scale_color_manual(name = "Age outlier", values = c("0" = "grey", "1" = "red"),
                     labels = c("0" = "0", "1" = "1")) +
  mytheme2+
  scale_y_continuous(breaks = c(2, 4, 6, 8, 10), limits=c(1,10)) +
  labs( x = "visit", y = ("PV3 log<sub>2</sub> neutralization titer"))+
  scale_x_discrete(labels= c("Baseline", "1 month", "3 months", "6 months", "1 year", "2 years", "3 years", "8-9 years"))+
  theme(plot.caption = element_text(hjust = 0.5, size = 14),legend.position = "none", axis.text.x = element_text(size = 18, angle=45, hjust=0.75), axis.text.y = element_text(size = 14), axis.title = element_markdown(size = 18), axis.title.x = element_blank())
ggsave("Figures/Age_outliers_PV3.png", width = 6, height = 6, dpi = 300)
agepv3_visit

agepv3 <- ggplot (longitudinal, aes(x=age_bopv14w, y=logPV3_stool))+
  geom_jitter(aes(color = factor(flag_ageoutlier)),alpha=0.6, position=position_jitter(width = 0.2, height = 0), size=2)+
  scale_color_manual(name = "Age outlier", values = c("0" = "grey", "1" = "red"),
                     labels = c("0" = "0", "1" = "1")) +
  mytheme2+
  scale_y_continuous(breaks = c(2, 4, 6, 8, 10), limits=c(1,10)) +
  labs( x = "Age at vaccination", y = ("PV3 log<sub>2</sub> neutralization titer"))+
  theme(plot.caption = element_text(hjust = 0.5, size = 14),legend.position = "none", axis.text.x = element_text(size = 18, angle=45, hjust=0.75), axis.text.y = element_text(size = 14), axis.title = element_markdown(size = 18), axis.title.x = element_blank())
#ggsave("Figures/Box_Stool_PV1_neut.png", width = 6, height = 6, dpi = 300)
agepv3

agepv2_visit <- ggplot (longitudinal, aes(x=factor(visit), y=logPV2_stool))+
  geom_jitter(aes(color = factor(flag_ageoutlier)),alpha=0.6, position=position_jitter(width = 0.2, height = 0), size=2)+
  scale_color_manual(name = "Age outlier", values = c("0" = "grey", "1" = "red"),
                     labels = c("0" = "0", "1" = "1")) +
  mytheme2+
  scale_y_continuous(breaks = c(2, 4, 6, 8, 10), limits=c(1,10)) +
  labs( x = "visit", y = ("PV2 log<sub>2</sub> neutralization titer"))+
  scale_x_discrete(labels= c("Baseline", "1 month", "3 months", "6 months", "1 year", "2 years", "3 years", "8-9 years"))+
  theme(plot.caption = element_text(hjust = 0.5, size = 14),legend.position = "none", axis.text.x = element_text(size = 18, angle=45, hjust=0.75), axis.text.y = element_text(size = 14), axis.title = element_markdown(size = 18), axis.title.x = element_blank())
ggsave("Figures/Age_outliers_PV2.png", width = 6, height = 6, dpi = 300)
agepv2_visit


#Remove the age outlier participants
longitudinal <- filter(longitudinal, flag_ageoutlier==0)


#################################################################################
## stool neutr plot: Box plots
#PV1 unstratified in longitudinal sample //3 are missing stool nuet data
m1nbox <- ggplot (longitudinal, aes(x=factor(visit), y=logPV1_stool))+
  geom_jitter(aes(color = factor(sampletype)),alpha=0.6, position=position_jitter(width = 0.2, height = 0), size=2)+
  geom_boxplot(fatten=3, outlier.shape = NA, color="black", fill=NA,  alpha=0.3)+
  scale_color_manual(name = "Sample", values = c("1" = "#1f78b4", "2" = "#E66218"),
    labels = c("1" = "Cross sectional", "2" = "Longitudinal")) +
   mytheme2+
  scale_y_continuous(breaks = c(2, 4, 6, 8, 10), limits=c(1,10)) +
  labs( x = "visit", y = ("PV1 log<sub>2</sub> neutralization titer"))+
  scale_x_discrete(labels= c("Baseline", "1 month", "3 months", "6 months", "1 year", "2 years", "3 years", "8-9 years"))+
 # theme(plot.caption = element_text(hjust = 0.5, size = 14), legend.text = element_text(size=14), legend.title = element_text(size=14), axis.text.x = element_text(size = 18, angle=45, hjust=0.75), axis.text.y = element_text(size = 14), axis.title = element_markdown(size = 18), axis.title.x = element_blank())
  theme(plot.caption = element_text(hjust = 0.5, size = 14),legend.position = "none", axis.text.x = element_text(size = 18, angle=45, hjust=0.75), axis.text.y = element_text(size = 14), axis.title = element_markdown(size = 18), axis.title.x = element_blank())
ggsave("Figures/Box_Stool_PV1_neut.png", width = 6, height = 6, dpi = 300)
m1nbox

#PV3 unstratified in longitdudinal sample
m3nbox <- ggplot (longitudinal, aes(x=factor(visit), y=logPV3_stool))+
  geom_jitter(aes(color = factor(sampletype)),alpha=0.6, position=position_jitter(width = 0.2, height = 0), size=2)+
  geom_boxplot(fatten=3, outlier.shape = NA, color="black", fill=NA,  alpha=0.3)+
  scale_color_manual(name = "Sample", values = c("1" = "#1f78b4", "2" = "#E66218"),
                     labels = c("1" = "Cross sectional", "2" = "Longitudinal")) +
  mytheme2+  
  labs( x = "visit", y = expression("PV3 log<sub>2</sub> neutralization titer"))+
  scale_x_discrete(labels= c("Baseline", "1 month", "3 months", "6 months", "1 year", "2 years", "3 years", "8-9 years"))+
  scale_y_continuous(breaks = c(0,2,4,6,8,10), limits=c(1,10)) +
  theme(plot.caption = element_text(hjust = 0.5, size = 14),legend.position = "none", axis.text.x = element_text(size = 18, angle=45, hjust=0.75), axis.text.y = element_text(size = 14), axis.title = element_markdown(size = 18), axis.title.x = element_blank()); m3nbox
ggsave("Figures/Box_Stool_PV3_neut.png", width = 6, height = 6, dpi = 300)

#PV2 unstratified in longitdudinal sample
m2nbox <- ggplot (longitudinal, aes(x=factor(visit), y=logPV2_stool))+
  geom_jitter(aes(color = factor(sampletype)),alpha=0.6, position=position_jitter(width = 0.2, height = 0), size=2)+
  geom_boxplot(fatten=3, outlier.shape = NA, color="black", fill=NA,  alpha=0.3)+
  scale_color_manual(name = "Sample", values = c("1" = "#1f78b4", "2" = "#E66218"),
                     labels = c("1" = "Cross sectional", "2" = "Longitudinal")) +
  mytheme2+
  labs( x = "visit", y = expression("PV2 log<sub>2</sub> neutralization titer"))+
  scale_x_discrete(labels= c("Baseline", "1 month", "3 months", "6 months", "1 year", "2 years", "3 years", "8-9 years"))+
  scale_y_continuous(breaks = c(0,2,4,6,8,10), limits=c(1,10)) +
  theme(plot.caption = element_text(hjust = 0.5, size = 14),legend.position = "none", axis.text.x = element_text(size = 18, angle=45, hjust=0.75), axis.text.y = element_text(size = 14), axis.title = element_markdown(size = 18), axis.title.x = element_blank())
ggsave("Figures/Box_Stool_PV2_neut.png", width = 6, height = 6, dpi = 300)
m2nbox

# #PV1 amongst strata with high pv1 baseline in longitudinal sample
# b<- ggplot (highpv1baseline, aes(x=factor(visit), y=logPV1_stool))+
#   geom_jitter(alpha=0.6, position=position_jitter(width = 0.2, height = 0), size=2)+
#   geom_boxplot(fatten=3, outlier.shape = NA, color="black", fill=NA,  alpha=0.3)+
#   mytheme2+
#   labs( x = "visit", y = expression("PV1 neutralization titer " [2]))+
#   scale_x_discrete(labels= c("Baseline", "One month", "3 months", "6 months"))+
#   ylim (1, 10)+
#   theme(plot.caption = element_text(hjust = 0.5, size = 18), legend.text = element_text(size=16), legend.title = element_text(size=16), axis.text = element_text(size = 16), axis.title = element_text(size = 16), axis.title.x = element_blank())
# b

# #PV1 amongst strata with low pv1 baseline in longitudinal sample
# c<- ggplot (lowpv1baseline, aes(x=factor(visit), y=logPV1_stool))+
#   geom_jitter(alpha=0.6, position=position_jitter(width = 0.2, height = 0), size=2)+
#   geom_boxplot(fatten=3, outlier.shape = NA, color="black", fill=NA,  alpha=0.3)+
#   mytheme2+
#   labs( x = "visit", y = expression("PV1 neutralization titer " [2]))+
#   scale_x_discrete(labels= c("Baseline", "One month", "3 months", "6 months"))+
#   ylim (1, 10)+
#   theme(plot.caption = element_text(hjust = 0.5, size = 18), legend.text = element_text(size=16), legend.title = element_text(size=16), axis.text = element_text(size = 16), axis.title = element_text(size = 16), axis.title.x = element_blank())
# c

#################################################################################
## Serum neutr plot: Box plots
#PV1 unstratified in longitudinal sample - 25 are missing serum values
s1nbox <- ggplot (longitudinal, aes(x=factor(visit), y=logPV1_serum))+
  geom_jitter(aes(color = factor(sampletype)),alpha=0.6, position=position_jitter(width = 0.2, height = 0), size=2)+
  geom_boxplot(fatten=3, outlier.shape = NA, color="black", fill=NA,  alpha=0.3)+
  scale_color_manual(name = "Sample", values = c("1" = "#1f78b4", "2" = "#E66218"),
                     labels = c("1" = "Cross sectional", "2" = "Longitudinal")) +
  mytheme2+
  labs( x = "visit", y = expression("PV1 log<sub>2</sub> neutralization titer"))+
  scale_x_discrete(labels= c("Baseline", "1 month", "3 months", "6 months", "1 year", "2 years", "3 years", "8-9 years"))+
  scale_y_continuous(breaks = c(0,2,4,6,8,10,12,14,16), limits=c(4,16)) +
  theme(plot.caption = element_text(hjust = 0.5, size = 14),legend.position = "none", axis.text.x = element_text(size = 18, angle=45, hjust=0.75), axis.text.y = element_text(size = 14), axis.title = element_markdown(size = 18), axis.title.x = element_blank())
ggsave("Figures/Box_Serum_PV1_neut.png", width = 6, height = 6, dpi = 300)
s1nbox

#PV3 unstratified in longitdudinal sample
s3nbox <- ggplot (longitudinal, aes(x=factor(visit), y=logPV3_serum))+
  geom_jitter(aes(color = factor(sampletype)),alpha=0.6, position=position_jitter(width = 0.2, height = 0), size=2)+
  geom_boxplot(fatten=3, outlier.shape = NA, color="black", fill=NA,  alpha=0.3)+
  scale_color_manual(name = "Sample", values = c("1" = "#1f78b4", "2" = "#E66218"),
                     labels = c("1" = "Cross sectional", "2" = "Longitudinal")) +
  mytheme2+
  labs( x = "visit", y = expression("PV3 log<sub>2</sub> neutralization titer"))+
  scale_x_discrete(labels= c("Baseline", "1 month", "3 months", "6 months", "1 year", "2 years", "3 years", "8-9 years"))+
  scale_y_continuous(breaks = c(0,2,4,6,8,10,12,14,16), limits=c(4,16)) +
  theme(plot.caption = element_text(hjust = 0.5, size = 14),legend.position = "none", axis.text.x = element_text(size = 18, angle=45, hjust=0.75), axis.text.y = element_text(size = 14), axis.title = element_markdown(size = 18), axis.title.x = element_blank())
ggsave("Figures/Box_Serum_PV3_neut.png", width = 6, height = 6, dpi = 300)
s3nbox

#PV2 unstratified in longitdudinal sample
s2nbox <- ggplot (longitudinal, aes(x=factor(visit), y=logPV2_serum))+
  geom_jitter(aes(color = factor(sampletype)),alpha=0.6, position=position_jitter(width = 0.2, height = 0), size=2)+
  geom_boxplot(fatten=3, outlier.shape = NA, color="black", fill=NA,  alpha=0.3)+
  scale_color_manual(name = "Sample", values = c("1" = "#1f78b4", "2" = "#E66218"),
                     labels = c("1" = "Cross sectional", "2" = "Longitudinal")) +
  mytheme2+
  labs( x = "visit", y = expression("PV2 log<sub>2</sub> neutralization titer"))+
  scale_x_discrete(labels= c("Baseline", "1 month", "3 months", "6 months", "1 year", "2 years", "3 years", "8-9 years"))+
  scale_y_continuous(breaks = c(0,2,4,6,8,10,12,14,16), limits=c(4,16)) +
  theme(plot.caption = element_text(hjust = 0.5, size = 14),legend.position = "none", axis.text.x = element_text(size = 18, angle=45, hjust=0.75), axis.text.y = element_text(size = 14), axis.title = element_markdown(size = 18), axis.title.x = element_blank())
ggsave("Figures/Box_Serum_PV2_neut.png", width = 6, height = 6, dpi = 300)
s2nbox
# 
# #PV1 amongst strata with high pv1 baseline in longitudinal sample
# b<- ggplot (highpv1baseline, aes(x=factor(visit), y=logPV1_serum))+
#   geom_jitter(alpha=0.6, position=position_jitter(width = 0.2, height = 0), size=2)+
#   geom_boxplot(fatten=3, outlier.shape = NA, color="black", fill=NA,  alpha=0.3)+
#   mytheme2+
#   labs( x = "visit", y = expression("PV1 neutralization titer " [2]))+
#   scale_x_discrete(labels= c("Baseline", "One month", "3 months", "6 months"))+
#   ylim (1, 10)+
#   theme(plot.caption = element_text(hjust = 0.5, size = 18), legend.text = element_text(size=16), legend.title = element_text(size=16), axis.text = element_text(size = 16), axis.title = element_text(size = 16), axis.title.x = element_blank())
# b
# 
# #PV1 amongst strata with low pv1 baseline in longitudinal sample
# c<- ggplot (lowpv1baseline, aes(x=factor(visit), y=logPV1_serum))+
#   geom_jitter(alpha=0.6, position=position_jitter(width = 0.2, height = 0), size=2)+
#   geom_boxplot(fatten=3, outlier.shape = NA, color="black", fill=NA,  alpha=0.3)+
#   mytheme2+
#   labs( x = "visit", y = expression("PV1 neutralization titer " [2]))+
#   scale_x_discrete(labels= c("Baseline", "One month", "3 months", "6 months"))+
#   ylim (1, 10)+
#   theme(plot.caption = element_text(hjust = 0.5, size = 18), legend.text = element_text(size=16), legend.title = element_text(size=16), axis.text = element_text(size = 16), axis.title = element_text(size = 16), axis.title.x = element_blank())
# c
