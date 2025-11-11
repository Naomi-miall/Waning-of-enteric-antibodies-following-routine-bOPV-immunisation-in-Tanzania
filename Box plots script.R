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
library(cowplot)
library(patchwork)

## mytheme
mytheme2 <- theme_bw() + 
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.border = element_blank(),
        legend.background = element_rect(fill=NA),
        legend.position = "bottom")

#################################################################################
#Neutralisation in stool
longitudinal <- read_dta("Box_plots.dta")
longitudinal <- subset (longitudinal, assay<3)
longitudinal$visit <- factor(longitudinal$visit, levels = c("1", "2", "3", "4", "5", "6", "7", "8"),
                             labels = c("Baseline", "One month", "3 months", "6 months", "1 year", "2 years", "3 years", "8-9 years"))
longitudinal$assay <- factor(longitudinal$assay, levels = c("1", "2"),
                             labels = c("Poliovirus type-1 in stool", "Poliovirus type-3 in stool"))

box_neut_stool <- ggplot(longitudinal, aes(x = factor(visit), y = logPV)) +
  geom_jitter(aes(color = factor(sampletype)),
              alpha = 0.6,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7),
              size = 2) +
  geom_boxplot(aes(group = interaction(sampletype, visit)),
               position = position_dodge2(width = 0.6, preserve = "single", padding=0.1),
               fatten = 3,
               width = 0.7, 
               outlier.shape = NA,
               color = "black",
               fill = NA,
               alpha = 0.3) +
  scale_color_manual(values = c("1" = "#1f78b4", "2" = "#E66218"),
                     name = "Cohort",
                     labels = c("1" = "Cross sectional", "2" = "Longitudinal")) +
  mytheme2 +
  scale_y_continuous(breaks = c(2, 6, 10, 14), limits = c(1, 14)) +
  labs(x = "Time since fourth dose of bOPV", y = "Log<sub>2</sub> neutralization titer") +
  scale_x_discrete(labels = c("Baseline", "1 month", "3 months", "6 months", "1 year", "2 years", "3 years", "8–9 years")) +
  facet_wrap(~assay) +
  theme(plot.caption = element_text(hjust = 0.5, size = 14),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_markdown(size = 16, face="bold"),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_line(color = "black"),
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 14),
        strip.background = element_rect(fill="white"),
        strip.text = element_text(face="bold", size=20),
        strip.placement = "outside",
        panel.spacing = unit(2, "lines"),
        plot.tag = element_text(face = "bold", size = 40))
box_neut_stool
###########################################################################################
#Neutralisation in serum
longitudinal <- read_dta("Box_plots.dta")
longitudinal <- subset (longitudinal, assay==3 | assay==4)
longitudinal$visit <- factor(longitudinal$visit, levels = c("1", "2", "3", "4", "5", "6", "7", "8"),
                             labels = c("Baseline", "One month", "3 months", "6 months", "1 year", "2 years", "3 years", "8-9 years"))
longitudinal$assay <- factor(longitudinal$assay, levels = c("3", "4"),
                             labels = c("Poliovirus type-1 in serum", "Poliovirus type-3 in serum"))
box_neut_serum <- ggplot(longitudinal, aes(x = factor(visit), y = logPV)) +
  geom_jitter(aes(color = factor(sampletype)),
              alpha = 0.6,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7),
              size = 2) +
  geom_boxplot(aes(group = interaction(sampletype, visit)),
               position = position_dodge2(width = 0.6, preserve = "single", padding=0.1),
               fatten = 3,
               width = 0.7, 
               outlier.shape = NA,
               color = "black",
               fill = NA,
               alpha = 0.3) +
  scale_color_manual(values = c("1" = "#1f78b4", "2" = "#E66218"),
                     name = "Cohort",
                     labels = c("1" = "Cross sectional", "2" = "Longitudinal")) +
  mytheme2 +
  scale_y_continuous(breaks = c(2, 6, 10, 14), limits = c(1, 14)) +
  labs(x = "Time since fourth dose of bOPV", y = "Log<sub>2</sub> neutralization titer") +
  scale_x_discrete(labels = c("Baseline", "1 month", "3 months", "6 months", "1 year", "2 years", "3 years", "8–9 years")) +
  facet_grid(~assay, scales="free", space="free") +
  theme(plot.caption = element_text(hjust = 0.5, size = 14),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_markdown(size = 16, face="bold"),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_line(color = "black"),
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 14),
        strip.background = element_rect(fill="white"),
        strip.text = element_text(face="bold", size=20),
        strip.placement = "outside",
        panel.spacing = unit(2, "lines"),
        plot.tag = element_text(face = "bold", size = 40))
box_neut_serum


###################################################################################################
#IGA stool
longitudinal <- read_dta("Box_plots.dta")
longitudinal <- subset (longitudinal, assay==5 | assay==6)
longitudinal$visit <- factor(longitudinal$visit, levels = c("1", "2", "3", "4", "5", "6", "7", "8"),
                             labels = c("Baseline", "One month", "3 months", "6 months", "1 year", "2 years", "3 years", "8-9 years"))
longitudinal$assay <- factor(longitudinal$assay, levels = c("5", "6"),
                             labels = c("Poliovirus type-1 IgA in stool", "Poliovirus type-3 IgA in stool"))
box_IgA_stool <- ggplot(longitudinal, aes(x = factor(visit), y = logPV)) +
  geom_jitter(aes(color = factor(sampletype)),
              alpha = 0.6,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7),
              size = 2) +
  geom_boxplot(aes(group = interaction(sampletype, visit)),
               position = position_dodge2(width = 0.6, preserve = "single", padding=0.1),
               fatten = 3,
               width = 0.7, 
               outlier.shape = NA,
               color = "black",
               fill = NA,
               alpha = 0.3) +
  scale_color_manual(values = c("1" = "#1f78b4", "2" = "#E66218"),
                     name = "Cohort",
                     labels = c("1" = "Cross sectional", "2" = "Longitudinal")) +
  mytheme2 +
  scale_y_continuous(breaks = c(1:4), limits = c(1, 4)) +
  labs(x = "Time since fourth dose of bOPV", y = "Log<sub>10</sub> IgA MFI") +
  scale_x_discrete(labels = c("Baseline", "1 month", "3 months", "6 months", "1 year", "2 years", "3 years", "8–9 years")) +
  facet_wrap(~assay,ncol=2) +
  theme(plot.caption = element_text(hjust = 0.5, size = 14),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_markdown(size = 16, face="bold"),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_line(color = "black"),
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 14),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.spacing = unit(2, "lines"),
        plot.tag = element_text(face = "bold", size = 40))
box_IgA_stool


###################################################################################################
#IGA Serum
longitudinal <- read_dta("Box_plots.dta")
longitudinal <- subset (longitudinal, assay==7 | assay==8)
longitudinal$visit <- factor(longitudinal$visit, levels = c("1", "2", "3", "4", "5", "6", "7", "8"),
                             labels = c("Baseline", "One month", "3 months", "6 months", "1 year", "2 years", "3 years", "8-9 years"))
longitudinal$assay <- factor(longitudinal$assay, levels = c("7", "8"),
                             labels = c("Poliovirus type-1 IgA in serum", "Poliovirus type-3 IgA in serum"))
box_IgA_serum <- ggplot(longitudinal, aes(x = factor(visit), y = logPV)) +
  geom_jitter(aes(color = factor(sampletype)),
              alpha = 0.6,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7),
              size = 2) +
  geom_boxplot(aes(group = interaction(sampletype, visit)),
               position = position_dodge2(width = 0.6, preserve = "single", padding=0.1),
               fatten = 3,
               width = 0.7, 
               outlier.shape = NA,
               color = "black",
               fill = NA,
               alpha = 0.3) +
  scale_color_manual(values = c("1" = "#1f78b4", "2" = "#E66218"),
                     name = "Cohort",
                     labels = c("1" = "Cross sectional", "2" = "Longitudinal")) +
  mytheme2 +
  scale_y_continuous(breaks = c(1:4), limits = c(1, 4)) +
  labs(x = "Time since fourth dose of bOPV", y = "Log<sub>10</sub> IgA MFI") +
  scale_x_discrete(labels = c("Baseline", "1 month", "3 months", "6 months", "1 year", "2 years", "3 years", "8–9 years")) +
  facet_wrap(~assay,ncol=2) +
  theme(plot.caption = element_text(hjust = 0.5, size = 14),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_markdown(size = 16, face="bold"),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_line(color = "black"),
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 14),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.spacing = unit(2, "lines"),
        plot.tag = element_text(face = "bold", size = 40))
box_IgA_serum


###################################################################################################
#IGG Serum
longitudinal <- read_dta("Box_plots.dta")
longitudinal <- subset (longitudinal, assay==9 | assay==10)
longitudinal$visit <- factor(longitudinal$visit, levels = c("1", "2", "3", "4", "5", "6", "7", "8"),
                             labels = c("Baseline", "One month", "3 months", "6 months", "1 year", "2 years", "3 years", "8-9 years"))
longitudinal$assay <- factor(longitudinal$assay, levels = c("9", "10"),
                             labels = c("Poliovirus type-1 IgG in serum", "Poliovirus type-3 IgG in serum"))
box_IgG_serum <- ggplot(longitudinal, aes(x = factor(visit), y = logPV)) +
  geom_jitter(aes(color = factor(sampletype)),
              alpha = 0.6,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7),
              size = 2) +
  geom_boxplot(aes(group = interaction(sampletype, visit)),
               position = position_dodge2(width = 0.6, preserve = "single", padding=0.1),
               fatten = 3,
               width = 0.7, 
               outlier.shape = NA,
               color = "black",
               fill = NA,
               alpha = 0.3) +
  scale_color_manual(values = c("1" = "#1f78b4", "2" = "#E66218"),
                     name = "Cohort",
                     labels = c("1" = "Cross sectional", "2" = "Longitudinal")) +
  mytheme2 +
  scale_y_continuous(breaks = c(1:4), limits = c(1, 4)) +
  labs(x = "Time since fourth dose of bOPV", y = "Log<sub>10</sub> IgG MFI") +
  scale_x_discrete(labels = c("Baseline", "1 month", "3 months", "6 months", "1 year", "2 years", "3 years", "8–9 years")) +
  facet_wrap(~assay,ncol=2) +
  theme(plot.caption = element_text(hjust = 0.5, size = 14),
        legend.position = "bottom",
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_markdown(size = 16, face="bold"),
        axis.title.x = element_text(size = 16, face="bold"),
        axis.line.x = element_line(color = "black"),
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 14),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.spacing = unit(2, "lines"),
        plot.tag = element_text(face = "bold", size = 40))
box_IgG_serum



box <- (box_neut_stool /  plot_spacer() / box_IgA_stool / plot_spacer() / box_neut_serum / plot_spacer() / box_IgA_serum / plot_spacer() / box_IgG_serum)+
  plot_annotation(tag_levels="a") +
  plot_layout(heights = c(1, 0.05, 1, 0.05, 1, 0.05, 1, 0.05, 1))
box
ggsave("Figures/Figure 1.png", width =12, height = 20, dpi = 300)









#PV2 unstratified in longitdudinal sample
longitudinal <- read_dta("Box_plots_PV2.dta")
longitudinal <- subset (longitudinal, assay<3 & visit!=NA)
longitudinal$visit <- factor(longitudinal$visit, levels = c("1", "2", "3", "4", "5", "6", "7", "8"),
                             labels = c("Baseline", "One month", "3 months", "6 months", "1 year", "2 years", "3 years", "8-9 years"))
longitudinal$assay <- factor(longitudinal$assay, levels = c("1", "2"),
                             labels = c("Poliovirus type 2 in stool", "Poliovirus type 2 in serum"))

box2 <- ggplot(longitudinal, aes(x = factor(visit), y = logPV)) +
  geom_jitter(aes(color = factor(sampletype)),
              alpha = 0.6,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7),
              size = 2) +
  geom_boxplot(aes(group = interaction(sampletype, visit)),
               position = position_dodge2(width = 0.6, preserve = "single", padding=0.1),
               fatten = 3,
               width = 0.7, 
               outlier.shape = NA,
               color = "black",
               fill = NA,
               alpha = 0.3) +
  scale_color_manual(values = c("1" = "#1f78b4", "2" = "#E66218"),
                     name = "Cohort",
                     labels = c("1" = "Cross sectional", "2" = "Longitudinal")) +
  mytheme2 +
  scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12, 14), limits = c(1, 14)) +
  labs(x = "Time since fourth dose of bOPV", y = "Log<sub>2</sub> neutralization titer") +
  scale_x_discrete(labels = c("Baseline", "1 month", "3 months", "6 months", "1 year", "2 years", "3 years", "8–9 years")) +
  facet_wrap(~assay) +
  theme(plot.caption = element_text(hjust = 0.5, size = 14),
        legend.position = "bottom",
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_markdown(size = 16, face="bold"),
        axis.title.x = element_text(size = 16, face="bold"),
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 14),
        strip.background = element_rect(fill="white"),
        strip.text = element_text(face="bold", size=16),
        panel.spacing = unit(2, "lines"))
box2
ggsave("Figures/Figure 1_PV2.png", width = 10, height = 5, dpi = 300)

#IgA box plot
longitudinal <- read_dta("Box_plots.dta")
longitudinal <- subset (longitudinal, assay>4)
longitudinal$visit <- factor(longitudinal$visit, levels = c("1", "2", "3", "4", "5", "6", "7", "8"),
                             labels = c("Baseline", "One month", "3 months", "6 months", "1 year", "2 years", "3 years", "8-9 years"))
longitudinal$assay <- factor(longitudinal$assay, levels = c("5", "6", "7", "8", "9", "10"),
                             labels = c("Poliovirus type-1 IgA in stool", "Poliovirus type-3 IgA in stool", "Poliovirus type-1 IgA in serum", "Poliovirus type-3 IgA in serum", "Poliovirus type-1 IgG in serum", "Poliovirus type-3 IgG in serum"))
boxIg <- ggplot(longitudinal, aes(x = factor(visit), y = logPV)) +
  geom_jitter(aes(color = factor(sampletype)),
              alpha = 0.6,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7),
              size = 2) +
  geom_boxplot(aes(group = interaction(sampletype, visit)),
               position = position_dodge2(width = 0.6, preserve = "single", padding=0.1),
               fatten = 3,
               width = 0.7, 
               outlier.shape = NA,
               color = "black",
               fill = NA,
               alpha = 0.3) +
  scale_color_manual(values = c("1" = "#1f78b4", "2" = "#E66218"),
                     name = "Cohort",
                     labels = c("1" = "Cross sectional", "2" = "Longitudinal")) +
  mytheme2 +
  scale_y_continuous(breaks = c(1:4), limits = c(1, 4)) +
  labs(x = "Time since fourth dose of bOPV", y = "Log<sub>10</sub> MFI") +
  scale_x_discrete(labels = c("Baseline", "1 month", "3 months", "6 months", "1 year", "2 years", "3 years", "8–9 years")) +
  facet_wrap(~assay,ncol=2) +
  theme(plot.caption = element_text(hjust = 0.5, size = 14),
        legend.position = "right",
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_markdown(size = 16, face="bold"),
        axis.title.x = element_text(size = 16, face="bold"),
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 14),
        strip.background = element_rect(fill="white"),
        strip.text = element_text(face="bold", size=16),
        panel.spacing = unit(2, "lines"))
boxIg
ggsave("Figures/Figure 1 part 2.png", width = 10, height = 10, dpi = 300)


###################################################################################################################################
# infection <- read_dta("infections_dataset.dta")
# infection$assay <- factor(infection$assay, levels = c("1", "2", "3", "4"),
#                              labels = c("PV1 stool", "PV3 stool", "PV1 serum", "PV3 serum"))
# 
# 
# s3nbox <- ggplot(infection, aes(x = factor(assay), y = logPV_)) +
#   geom_hline(yintercept = 3, linetype = "dashed") +
#   geom_jitter(aes(color = factor(any_infection)),
#               alpha = 0.6,
#               position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6),
#               size = 2) +
#   geom_boxplot(aes(group = interaction(any_infection, assay)),
#                position = position_dodge2(width = 0.6, preserve = "single", padding=0.3),
#                fatten = 2,
#                width = 0.4, 
#                outlier.shape = NA,
#                color = "black",
#                fill = NA,
#                alpha = 0.3) +
#   scale_color_manual(name = "Infections status", values = c("0" = "#00798C", "1" = "#D1495B"),
#                        labels = c("0" = "None", "1" = "At least 1")) +
#   mytheme2 +
#   scale_y_continuous(breaks = c(0,2,4,6,8,10,12,14), limits=c(0,14)) +
#   labs(x = "Sample", y = "Log<sub>2</sub> neutralization titer") +
#   scale_x_discrete(labels = c("PV1 in stool", "PV3 in stool", "PV1 in serum", "PV3 in serum")) +
#   theme(plot.caption = element_text(hjust = 0.5, size = 14),
#         legend.position = "right",
#         axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
#         axis.text.y = element_text(size = 14),
#         axis.title = element_markdown(size = 18),
#         axis.title.x = element_blank())
# s3nbox
# ggsave("Figures/Box_infections_neut.png", width = 6, height = 6, dpi = 300)

########################################################################################
#With IgA
longitudinal <- read_dta("Box_plots.dta")

longitudinal$visit <- factor(longitudinal$visit, levels = c("1", "2", "3", "4", "5", "6", "7", "8"),
                             labels = c("Baseline", "One month", "3 months", "6 months", "1 year", "2 years", "3 years", "8-9 years"))

longitudinal$assay <- factor(longitudinal$assay, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"),
                             labels = c("Poliovirus type 1 log2 neutralising titre in stool", "Poliovirus type 3 log2 neutralising titre in stool", "Poliovirus type 1 log2 neutralising titre in serum", "Poliovirus type 3 log2 neutralising titre in serum","Poliovirus type 1 log10 IgA MFI in serum", "Poliovirus type 3 log10 IgA MFI in serum","Poliovirus type 1 log10 IgA MFI in stool", "Poliovirus type 3 log10 IgA MFI in stool", "Poliovirus type 1 log10 IgG MFI in serum", "Poliovirus type 3 log10 IgG MFI in serum"))

## stool neutr plot: Box plots
box <- ggplot(longitudinal, aes(x = factor(visit), y = logPV)) +
  geom_jitter(aes(color = factor(sampletype)),
              alpha = 0.6,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7),
              size = 2) +
  geom_boxplot(aes(group = interaction(sampletype, visit)),
               position = position_dodge2(width = 0.6, preserve = "single", padding=0.1),
               fatten = 3,
               width = 0.7, 
               outlier.shape = NA,
               color = "black",
               fill = NA,
               alpha = 0.3) +
  scale_color_manual(values = c("1" = "#1f78b4", "2" = "#E66218"),
                     name = "Cohort",
                     labels = c("1" = "Cross sectional", "2" = "Longitudinal")) +
  mytheme2 +
  scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12, 14), limits = c(1, 14)) +
  labs(x = "Time since fourth dose of bOPV", y = "Log<sub>2</sub> neutralization titer") +
  scale_x_discrete(labels = c("Baseline", "1 month", "3 months", "6 months", "1 year", "2 years", "3 years", "8–9 years")) +
  facet_wrap(~assay) +
  theme(plot.caption = element_text(hjust = 0.5, size = 14),
        legend.position = "right",
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_markdown(size = 16, face="bold"),
        axis.title.x = element_text(size = 16, face="bold"),
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 14),
        strip.background = element_rect(fill="white"),
        strip.text = element_text(face="bold", size=16),
        panel.spacing = unit(2, "lines"))
box
ggsave("Figures/Figure 1.png", width = 10, height = 10, dpi = 300)
