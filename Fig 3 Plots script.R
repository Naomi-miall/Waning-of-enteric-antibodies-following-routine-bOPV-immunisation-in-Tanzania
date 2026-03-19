###R script creating scatter plots for TZA duration of immunity paper
#Naomi and Alex
#October 2025

#Set up: Libaries; working directory; theme:
### scatter plot ####
setwd ("C:/Users/_________________________________________")

## library
library(haven)
library(ggplot2)
library(ggpubr)
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

longitudinal <- read_dta("merged_dataset_long.dta")

longitudinal$visit <- factor(longitudinal$visit, levels = c("1", "2", "3", "4"),
                             labels = c("Baseline", "One month", "3 months", "6 months"))

#Create data subsets for stratification by baseline; by complete data or not
#PV1 baseline stool - high means >=8 at baseline
highpv1baseline <- subset (longitudinal, high_neut_baseline_pv1==1)
lowpv1baseline <- subset (longitudinal, high_neut_baseline_pv1==0)

#PV3 baseline stool - high means >=8 at baseline
highpv3baseline <- subset (longitudinal, high_neut_baseline_pv3==1)
lowpv3baseline <- subset (longitudinal, high_neut_baseline_pv3==0)

#Complete data
completedata <- subset (longitudinal, no_obs==4 & flag_missing_stool_neuts==0 & flag_missing_serum_neuts==0)
completedata_highpv1baseline <- subset (completedata, high_neut_baseline_pv1==1)
completedata_lowpv1baseline <- subset (completedata, high_neut_baseline_pv1==0)
completedata_highpv3baseline <- subset (completedata, high_neut_baseline_pv3==1)
completedata_lowpv3baseline <- subset (completedata, high_neut_baseline_pv3==0)


visit_labels <- c("Baseline", "1 month", "3 months", "6 months")

##################################################################################
#################################################################################
## stool neutr plot: Spagetti plots

# #Including those with missing data:
# #Unstratified
mean_longitudinal <- longitudinal %>%
  group_by(summ_visit_month) %>%
  summarise(mean_logpv1_stool = mean(logPV1_stool, na.rm = TRUE),
            mean_logpv3_stool = mean(logPV3_stool, na.rm = TRUE),
            mean_logpv2_stool = mean(logPV2_stool, na.rm=TRUE),
            mean_logpv1_serum = mean(logPV1_serum, na.rm = TRUE),
            mean_logpv3_serum = mean(logPV3_serum, na.rm = TRUE),
            mean_logpv2_serum = mean(logPV2_serum, na.rm=TRUE),
            median_logpv2_stool = median(logPV2_stool, na.rm = TRUE),
            median_logpv2_serum = median(logPV2_serum, na.rm = TRUE),
            lq_logpv2_stool = quantile(logPV2_stool, 0.25, na.rm = TRUE),
            lq_logpv2_serum = quantile(logPV2_serum, 0.25, na.rm = TRUE),
            uq_logpv2_stool = quantile(logPV2_stool, 0.75, na.rm = TRUE),
            uq_logpv2_serum = quantile(logPV2_serum, 0.75, na.rm = TRUE))


m2n <- ggplot(longitudinal, aes(x =follow_up_months, y = logPV2_stool, group = pid)) +
  geom_line(alpha = 0.5, size = 0.2) +
  geom_point(alpha = 0.1) +
  geom_hline(yintercept = 3, linetype = "dashed")+
  geom_ribbon(data = mean_longitudinal, aes(x = summ_visit_month, ymin = lq_logpv2_stool, ymax = uq_logpv2_stool), fill = "#725AC1", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = mean_longitudinal, aes(x = summ_visit_month, y = median_logpv2_stool), color = "#725AC1", inherit.aes = FALSE, size=1) +
  geom_point(data = mean_longitudinal, aes(x = summ_visit_month, y = median_logpv2_stool), shape = 17, size = 2, color = "#725AC1", inherit.aes = FALSE) +
  scale_y_continuous(breaks = c(1, 3, 5, 7, 9, 11, 13), limits = c(1,14)) +
  scale_x_continuous(breaks = c(0, 1, 3, 6), limits = c(1,7.4), labels=visit_labels) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 18),
        axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
ggsave("Figures/Stool_PV2_all_median.png", width = 6, height = 4, dpi = 300)
m2n

#High pv1 baseline full sample
mean_highpv1baseline <- highpv1baseline %>%
  group_by(summ_visit_month) %>%
  summarise(mean_logpv1_stool = mean(logPV1_stool, na.rm = TRUE),
            mean_logpv1_serum = mean(logPV1_serum, na.rm = TRUE),
            median_logpv1_stool = median(logPV1_stool, na.rm = TRUE),
            median_logpv1_serum = median(logPV1_serum, na.rm = TRUE),
            lq_logpv1_stool = quantile(logPV1_stool, 0.25, na.rm = TRUE),
            lq_logpv1_serum = quantile(logPV1_serum, 0.25, na.rm = TRUE),
            uq_logpv1_stool = quantile(logPV1_stool, 0.75, na.rm = TRUE),
            uq_logpv1_serum = quantile(logPV1_serum, 0.75, na.rm = TRUE),
            median_logpv1_stool_IgA = median(logPV1_IgA_stool, na.rm=TRUE),
            median_logpv1_serum_IgA = median(logPV1_IgA_serum, na.rm=TRUE),
            median_logpv1_serum_IgG = median(logPV1_IgG_serum, na.rm=TRUE),
            lq_logpv1_stool_IgA = quantile(logPV1_IgA_stool, 0.25, na.rm = TRUE),
            lq_logpv1_serum_IgA = quantile(logPV1_IgA_serum, 0.25, na.rm = TRUE),
            lq_logpv1_serum_IgG = quantile(logPV1_IgG_serum, 0.25, na.rm = TRUE),
            uq_logpv1_stool_IgA = quantile(logPV1_IgA_stool, 0.75, na.rm = TRUE),
            uq_logpv1_serum_IgA = quantile(logPV1_IgA_serum, 0.75, na.rm = TRUE),
            uq_logpv1_serum_IgG = quantile(logPV1_IgG_serum, 0.75, na.rm = TRUE),
            )

#Stool PV1 neuts high baseline
m1nhb <- ggplot(highpv1baseline, aes(x = follow_up_months, y = logPV1_stool, group = pid)) +
  geom_line(alpha = 0.5, size = 0.2) + geom_point(alpha = 0.1) +
  geom_hline(yintercept = 3, linetype = "dashed") +
  geom_ribbon(data = mean_highpv1baseline, aes(x = summ_visit_month, ymin = lq_logpv1_stool, ymax = uq_logpv1_stool), fill = "#725AC1", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = mean_highpv1baseline, aes(x = summ_visit_month, y = median_logpv1_stool), color = "#725AC1", inherit.aes = FALSE, size=1) +
  geom_point(data = mean_highpv1baseline, aes(x = summ_visit_month, y = median_logpv1_stool), shape = 17, size = 2, color = "#725AC1", inherit.aes = FALSE) +
  scale_y_continuous(breaks = c(1, 3, 5, 7, 9, 11, 13), limits = c(1,14)) +
  scale_x_continuous(breaks = c(0, 1, 3, 6), limits = c(0, 7.4), labels=visit_labels) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 18),
        #axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
        axis.text.x = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
m1nhb
ggsave("Figures/Stool_PV1_Neuts_highbaseline_all_medians.png", width = 6, height = 4, dpi = 300)

#Stool PV1 IgA high baseline
m1igahb <- ggplot(highpv1baseline, aes(x = follow_up_months, y = logPV1_IgA_stool, group = pid)) +
  geom_line(alpha = 0.5, size = 0.2) + geom_point(alpha = 0.1) +
  geom_ribbon(data = mean_highpv1baseline, aes(x = summ_visit_month, ymin = lq_logpv1_stool_IgA, ymax = uq_logpv1_stool_IgA), fill = "#725AC1", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = mean_highpv1baseline, aes(x = summ_visit_month, y = median_logpv1_stool_IgA), color = "#725AC1", inherit.aes = FALSE, size=1) +
  geom_point(data = mean_highpv1baseline, aes(x = summ_visit_month, y = median_logpv1_stool_IgA), shape = 17, size = 2, color = "#725AC1", inherit.aes = FALSE) +
  scale_y_continuous(breaks = c(0:4), limits = c(-1,4)) +
  scale_x_continuous(breaks = c(0, 1, 3, 6), limits = c(0,7.4), labels=visit_labels) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
m1igahb
ggsave("Figures/Stool_PV1_IgA_highbaseline_all_medians.png", width = 6, height = 8, dpi = 300)


#Low PV1 baseline full sample
mean_lowpv1baseline <- lowpv1baseline %>%
  group_by(summ_visit_month) %>%
  summarise(mean_logpv1_stool = mean(logPV1_stool, na.rm = TRUE),
            mean_logpv1_serum = mean(logPV1_serum, na.rm = TRUE),
            median_logpv1_stool = median(logPV1_stool, na.rm = TRUE),
            median_logpv1_serum = median(logPV1_serum, na.rm = TRUE),
            lq_logpv1_stool = quantile(logPV1_stool, 0.25, na.rm = TRUE),
            lq_logpv1_serum = quantile(logPV1_serum, 0.25, na.rm = TRUE),
            uq_logpv1_stool = quantile(logPV1_stool, 0.75, na.rm = TRUE),
            uq_logpv1_serum = quantile(logPV1_serum, 0.75, na.rm = TRUE),
            median_logpv1_stool_IgA = median(logPV1_IgA_stool, na.rm=TRUE),
            median_logpv1_serum_IgA = median(logPV1_IgA_serum, na.rm=TRUE),
            median_logpv1_serum_IgG = median(logPV1_IgG_serum, na.rm=TRUE),
            lq_logpv1_stool_IgA = quantile(logPV1_IgA_stool, 0.25, na.rm = TRUE),
            lq_logpv1_serum_IgA = quantile(logPV1_IgA_serum, 0.25, na.rm = TRUE),
            lq_logpv1_serum_IgG = quantile(logPV1_IgG_serum, 0.25, na.rm = TRUE),
            uq_logpv1_stool_IgA = quantile(logPV1_IgA_stool, 0.75, na.rm = TRUE),
            uq_logpv1_serum_IgA = quantile(logPV1_IgA_serum, 0.75, na.rm = TRUE),
            uq_logpv1_serum_IgG = quantile(logPV1_IgG_serum, 0.75, na.rm = TRUE),
  )

#Stool PV1 neuts low baseline
m1nlb <- ggplot(lowpv1baseline, aes(x = follow_up_months, y = logPV1_stool, group = pid)) +
  geom_line(alpha = 0.5, size = 0.2) + geom_point(alpha = 0.1) +
  geom_hline(yintercept = 3, linetype = "dashed") +
  geom_ribbon(data = mean_lowpv1baseline, aes(x = summ_visit_month, ymin = lq_logpv1_stool, ymax = uq_logpv1_stool), fill = "#725AC1", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = mean_lowpv1baseline, aes(x = summ_visit_month, y = median_logpv1_stool), color = "#725AC1", inherit.aes = FALSE, size=1) +
  geom_point(data = mean_lowpv1baseline, aes(x = summ_visit_month, y = median_logpv1_stool), shape = 17, size = 2, color = "#725AC1", inherit.aes = FALSE) +
  scale_y_continuous(breaks = c(1, 3, 5, 7, 9, 11, 13), limits = c(1,14)) +
  scale_x_continuous(breaks = c(0, 1, 3, 6), limits = c(0,7.4), labels=visit_labels) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 18),
        #axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
        axis.text.x = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
m1nlb
ggsave("Figures/Stool_PV1_Neuts_lowbaseline_all_medians.png", width = 6, height = 4, dpi = 300)

#Stool PV1 IgA low baseline
m1igalb <- ggplot(lowpv1baseline, aes(x = follow_up_months, y = logPV1_IgA_stool, group = pid)) +
  geom_line(alpha = 0.5, size = 0.2) + geom_point(alpha = 0.1) +
  geom_ribbon(data = mean_lowpv1baseline, aes(x = summ_visit_month, ymin = lq_logpv1_stool_IgA, ymax = uq_logpv1_stool_IgA), fill = "#725AC1", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = mean_lowpv1baseline, aes(x = summ_visit_month, y = median_logpv1_stool_IgA), color = "#725AC1", inherit.aes = FALSE, size=1) +
  geom_point(data = mean_lowpv1baseline, aes(x = summ_visit_month, y = median_logpv1_stool_IgA), shape = 17, size = 2, color = "#725AC1", inherit.aes = FALSE) +
  scale_y_continuous(breaks = c(0:4), limits = c(-1,4)) +
  scale_x_continuous(breaks = c(0, 1, 3, 6), limits = c(0,7.4), labels=visit_labels) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
m1igalb
ggsave("Figures/Stool_PV1_IgA_lowbaseline_all_medians.png", width = 6, height = 8, dpi = 300)



#High pv3 baseline full sample
mean_highpv3baseline <- highpv3baseline %>%
  group_by(summ_visit_month) %>%
  summarise(mean_logpv3_stool = mean(logPV3_stool, na.rm = TRUE),
            mean_logpv3_serum = mean(logPV3_serum, na.rm = TRUE),
            median_logpv3_stool = median(logPV3_stool, na.rm = TRUE),
            median_logpv3_serum = median(logPV3_serum, na.rm = TRUE),
            lq_logpv3_stool = quantile(logPV3_stool, 0.25, na.rm = TRUE),
            lq_logpv3_serum = quantile(logPV3_serum, 0.25, na.rm = TRUE),
            uq_logpv3_stool = quantile(logPV3_stool, 0.75, na.rm = TRUE),
            uq_logpv3_serum = quantile(logPV3_serum, 0.75, na.rm = TRUE),
            median_logpv3_stool_IgA = median(logPV3_IgA_stool, na.rm=TRUE),
            median_logpv3_serum_IgA = median(logPV3_IgA_serum, na.rm=TRUE),
            median_logpv3_serum_IgG = median(logPV3_IgG_serum, na.rm=TRUE),
            lq_logpv3_stool_IgA = quantile(logPV3_IgA_stool, 0.25, na.rm = TRUE),
            lq_logpv3_serum_IgA = quantile(logPV3_IgA_serum, 0.25, na.rm = TRUE),
            lq_logpv3_serum_IgG = quantile(logPV3_IgG_serum, 0.25, na.rm = TRUE),
            uq_logpv3_stool_IgA = quantile(logPV3_IgA_stool, 0.75, na.rm = TRUE),
            uq_logpv3_serum_IgA = quantile(logPV3_IgA_serum, 0.75, na.rm = TRUE),
            uq_logpv3_serum_IgG = quantile(logPV3_IgG_serum, 0.75, na.rm = TRUE),
  )


#Stool PV3 neuts high baseline
m3nhb <- ggplot(highpv3baseline, aes(x = follow_up_months, y = logPV3_stool, group = pid)) +
  geom_line(alpha = 0.5, size = 0.2) + geom_point(alpha = 0.1) +
  geom_hline(yintercept = 3, linetype = "dashed") +
  geom_ribbon(data = mean_highpv3baseline, aes(x = summ_visit_month, ymin = lq_logpv3_stool, ymax = uq_logpv3_stool), fill = "#725AC1", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = mean_highpv3baseline, aes(x = summ_visit_month, y = median_logpv3_stool), color = "#725AC1", inherit.aes = FALSE, size=1) +
  geom_point(data = mean_highpv3baseline, aes(x = summ_visit_month, y = median_logpv3_stool), shape = 17, size = 2, color = "#725AC1", inherit.aes = FALSE) +
  scale_y_continuous(breaks = c(1, 3, 5, 7, 9, 11, 13), limits = c(1,14)) +
  scale_x_continuous(breaks = c(0, 1, 3, 6), limits = c(0,7.4), labels=visit_labels) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 18),
        # axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
        axis.text.x = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
ggsave("Figures/Stool_PV3_Neuts_highbaseline_all_medians.png", width = 6, height = 4, dpi = 300)
m3nhb 

#Stool PV3 IgA high baseline
m3igahb <- ggplot(highpv3baseline, aes(x = follow_up_months, y  = logPV3_IgA_stool, group = pid)) +
  geom_line(alpha = 0.5, size = 0.2) + geom_point(alpha = 0.1) +
  geom_ribbon(data = mean_highpv3baseline, aes(x = summ_visit_month, ymin = lq_logpv3_stool_IgA, ymax = uq_logpv3_stool_IgA), fill = "#725AC1", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = mean_highpv3baseline, aes(x = summ_visit_month, y = median_logpv3_stool_IgA), color = "#725AC1", inherit.aes = FALSE, size=1) +
  geom_point(data = mean_highpv3baseline, aes(x = summ_visit_month, y = median_logpv3_stool_IgA), shape = 17, size = 2, color = "#725AC1", inherit.aes = FALSE) +
  scale_y_continuous(breaks = c(0:4), limits = c(-1,4)) +
  scale_x_continuous(breaks = c(0, 1, 3, 6), limits = c(0,7.4), labels=visit_labels) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
m3igahb
ggsave("Figures/Stool_PV3_IgA_highbaseline_all_medians.png", width = 6, height = 8, dpi = 300)


#Low PV3 baseline full sample
mean_lowpv3baseline <- lowpv3baseline %>%
  group_by(summ_visit_month) %>%
  summarise(mean_logpv3_stool = mean(logPV3_stool, na.rm = TRUE),
            mean_logpv3_serum = mean(logPV3_serum, na.rm = TRUE),
            median_logpv3_stool = median(logPV3_stool, na.rm = TRUE),
            median_logpv3_serum = median(logPV3_serum, na.rm = TRUE),
            lq_logpv3_stool = quantile(logPV3_stool, 0.25, na.rm = TRUE),
            lq_logpv3_serum = quantile(logPV3_serum, 0.25, na.rm = TRUE),
            uq_logpv3_stool = quantile(logPV3_stool, 0.75, na.rm = TRUE),
            uq_logpv3_serum = quantile(logPV3_serum, 0.75, na.rm = TRUE),
            median_logpv3_stool_IgA = median(logPV3_IgA_stool, na.rm=TRUE),
            median_logpv3_serum_IgA = median(logPV3_IgA_serum, na.rm=TRUE),
            median_logpv3_serum_IgG = median(logPV3_IgG_serum, na.rm=TRUE),
            lq_logpv3_stool_IgA = quantile(logPV3_IgA_stool, 0.25, na.rm = TRUE),
            lq_logpv3_serum_IgA = quantile(logPV3_IgA_serum, 0.25, na.rm = TRUE),
            lq_logpv3_serum_IgG = quantile(logPV3_IgG_serum, 0.25, na.rm = TRUE),
            uq_logpv3_stool_IgA = quantile(logPV3_IgA_stool, 0.75, na.rm = TRUE),
            uq_logpv3_serum_IgA = quantile(logPV3_IgA_serum, 0.75, na.rm = TRUE),
            uq_logpv3_serum_IgG = quantile(logPV3_IgG_serum, 0.75, na.rm = TRUE),
  )

#Stool PV3 neuts low baseline
m3nlb <- ggplot(lowpv3baseline, aes(x =follow_up_months, y = logPV3_stool, group = pid)) +
  geom_line(alpha = 0.5, size = 0.2) +
  geom_point(alpha = 0.1) +
  geom_hline(yintercept = 3, linetype = "dashed")+
  geom_ribbon(data = mean_lowpv3baseline, aes(x = summ_visit_month, ymin = lq_logpv3_stool, ymax = uq_logpv3_stool), fill = "#725AC1", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = mean_lowpv3baseline, aes(x = summ_visit_month, y = median_logpv3_stool), color = "#725AC1", inherit.aes = FALSE, size=1) +
  geom_point(data = mean_lowpv3baseline, aes(x = summ_visit_month, y = median_logpv3_stool), shape = 17, size = 2, color = "#725AC1", inherit.aes = FALSE) +
  scale_y_continuous(breaks = c(1, 3, 5, 7, 9, 11, 13), limits = c(1,14)) +
  scale_x_continuous(breaks = c(0, 1, 3, 6), limits = c(0,7.4), labels=visit_labels) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 18),
        #axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
        axis.text.x = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
m3nlb
ggsave("Figures/Stool_PV3_Neuts_lowbaseline_all_medians.png", width = 6, height = 4, dpi = 300)

#Stool PV3 IgA low baseline
m3igalb <- ggplot(lowpv3baseline, aes(x = follow_up_months, y = logPV3_IgA_stool, group = pid)) +
  geom_line(alpha = 0.5, size = 0.2) + geom_point(alpha = 0.1) +
  geom_ribbon(data = mean_lowpv3baseline, aes(x = summ_visit_month, ymin = lq_logpv3_stool_IgA, ymax = uq_logpv3_stool_IgA), fill = "#725AC1", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = mean_lowpv3baseline, aes(x = summ_visit_month, y = median_logpv3_stool_IgA), color = "#725AC1", inherit.aes = FALSE, size=1) +
  geom_point(data = mean_lowpv3baseline, aes(x = summ_visit_month, y = median_logpv3_stool_IgA), shape = 17, size = 2, color = "#725AC1", inherit.aes = FALSE) +
  scale_y_continuous(breaks = c(0:4), limits = c(-1,4)) +
  scale_x_continuous(breaks = c(0, 1, 3, 6), limits = c(0,7.4), labels=visit_labels) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
m3igalb
ggsave("Figures/Stool_PV3_IgA_lowbaseline_all_medians.png", width = 6, height = 8, dpi = 300)

#################################################################################
#################################################################################
## Serum neutr plot: Spagetti plots

#All participants unstratified for PV2
s2n <- ggplot(longitudinal, aes(x =follow_up_months, y = logPV2_serum, group = pid)) +
  geom_line(alpha = 0.5, size = 0.2) +
  geom_point(alpha = 0.1) +  
  geom_hline(yintercept = 3, linetype = "dashed")+
  geom_ribbon(data = mean_longitudinal, aes(x = summ_visit_month, ymin = lq_logpv2_serum, ymax = uq_logpv2_serum), fill = "#BC3908", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = mean_longitudinal, aes(x = summ_visit_month, y = median_logpv2_serum), color = "#BC3908", inherit.aes = FALSE, size=1) +
  geom_point(data = mean_longitudinal, aes(x = summ_visit_month, y = median_logpv2_serum), shape = 17, size = 2, color = "#BC3908", inherit.aes = FALSE) +
  scale_y_continuous(breaks = c(1, 3, 5, 7, 9, 11, 13), limits = c(1,14)) +
  scale_x_continuous(breaks = c(0, 1, 3, 6), limits = c(0,7.4), labels=visit_labels) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 18),
        axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
ggsave("Figures/Serum_PV2_all_median.png", width = 6, height = 4, dpi = 300)
s2n


#Serum PV1 neuts high baseline
s1nhb<- ggplot(highpv1baseline, aes(x =follow_up_months, y = logPV1_serum, group = pid)) +
  geom_line(alpha = 0.5, size = 0.2) +
  geom_point(alpha = 0.1) +
  geom_hline(yintercept = 3, linetype = "dashed")+
  geom_ribbon(data = mean_highpv1baseline, aes(x = summ_visit_month, ymin = lq_logpv1_serum, ymax = uq_logpv1_serum), fill = "#BC3908", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = mean_highpv1baseline, aes(x = summ_visit_month, y = median_logpv1_serum), color = "#BC3908", inherit.aes = FALSE, size=1) +
  geom_point(data = mean_highpv1baseline, aes(x = summ_visit_month, y = median_logpv1_serum), shape = 17, size = 2, color = "#BC3908", inherit.aes = FALSE) +
  scale_y_continuous(breaks = c(1, 3, 5, 7, 9, 11, 13), limits = c(1,14)) +
  scale_x_continuous(breaks = c(0, 1, 3, 6), limits = c(0,7.4), labels=visit_labels) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 18),
        #axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
        axis.text.x = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
ggsave("Figures/Serum_PV1_Neuts_highbaseline_all_medians.png", width = 6.3, height = 4, dpi = 300)
s1nhb 

#Serum PV1 IgA high baseline
s1igahb<- ggplot(highpv1baseline, aes(x =follow_up_months, y = logPV1_IgA_serum, group = pid)) +
  geom_line(alpha = 0.5, size = 0.2) +
  geom_point(alpha = 0.1) +
  geom_ribbon(data = mean_highpv1baseline, aes(x = summ_visit_month, ymin = lq_logpv1_serum_IgA, ymax = uq_logpv1_serum_IgA), fill = "#BC3908", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = mean_highpv1baseline, aes(x = summ_visit_month, y = median_logpv1_serum_IgA), color = "#BC3908", inherit.aes = FALSE, size=1) +
  geom_point(data = mean_highpv1baseline, aes(x = summ_visit_month, y = median_logpv1_serum_IgA), shape = 17, size = 2, color = "#BC3908", inherit.aes = FALSE) +
  scale_y_continuous(breaks = c(0:4), limits = c(-1,4)) +
  scale_x_continuous(breaks = c(0, 1, 3, 6), limits = c(0,7.4), labels=visit_labels) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
ggsave("Figures/Serum_PV1_IgA_highbaseline_all_medians.png", width = 6.3, height = 8, dpi = 300)
s1igahb 

#Serum PV1 IgG high baseline
s1igghb<- ggplot(highpv1baseline, aes(x =follow_up_months, y = logPV1_IgG_serum, group = pid)) +
  geom_line(alpha = 0.5, size = 0.2) +
  geom_point(alpha = 0.1) +
  geom_ribbon(data = mean_highpv1baseline, aes(x = summ_visit_month, ymin = lq_logpv1_serum_IgG, ymax = uq_logpv1_serum_IgG), fill = "#BC3908", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = mean_highpv1baseline, aes(x = summ_visit_month, y = median_logpv1_serum_IgG), color = "#BC3908", inherit.aes = FALSE, size=1) +
  geom_point(data = mean_highpv1baseline, aes(x = summ_visit_month, y = median_logpv1_serum_IgG), shape = 17, size = 2, color = "#BC3908", inherit.aes = FALSE) +
  scale_y_continuous(breaks = c(0:4), limits = c(-1,4)) +
  scale_x_continuous(breaks = c(0, 1, 3, 6), limits = c(0,7.4), labels=visit_labels) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 18),
        axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
ggsave("Figures/Serum_PV1_IgG_highbaseline_all_medians.png", width = 6.3, height = 9, dpi = 300)
s1igghb

#Serum PV1 neuts low baseline
s1nlb<- ggplot(lowpv1baseline, aes(x =follow_up_months, y = logPV1_serum, group = pid)) +
  geom_line(alpha = 0.5, size = 0.2) +
  geom_point(alpha = 0.1) +
  geom_hline(yintercept = 3, linetype = "dashed")+
  geom_ribbon(data = mean_lowpv1baseline, aes(x = summ_visit_month, ymin = lq_logpv1_serum, ymax = uq_logpv1_serum), fill = "#BC3908", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = mean_lowpv1baseline, aes(x = summ_visit_month, y = median_logpv1_serum), color = "#BC3908",  inherit.aes = FALSE, size=1) +
  geom_point(data = mean_lowpv1baseline, aes(x = summ_visit_month, y = median_logpv1_serum), shape = 17, size = 2, color = "#BC3908", inherit.aes = FALSE) +
  scale_y_continuous(breaks = c(1, 3, 5, 7, 9, 11, 13), limits = c(1,14)) +
  scale_x_continuous(breaks = c(0, 1, 3, 6), limits = c(0,7.4), labels=visit_labels) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 18),
        #axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
        axis.text.x = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
ggsave("Figures/Serum_PV1_Neuts_lowbaseline_all_medians.png", width = 6.3, height = 4, dpi = 300)
s1nlb


#Serum PV1 IgA low baseline
s1igalb<- ggplot(lowpv1baseline, aes(x =follow_up_months, y = logPV1_IgA_serum, group = pid)) +
  geom_line(alpha = 0.5, size = 0.2) +
  geom_point(alpha = 0.1) +
  geom_ribbon(data = mean_lowpv1baseline, aes(x = summ_visit_month, ymin = lq_logpv1_serum_IgA, ymax = uq_logpv1_serum_IgA), fill = "#BC3908", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = mean_lowpv1baseline, aes(x = summ_visit_month, y = median_logpv1_serum_IgA), color = "#BC3908", inherit.aes = FALSE, size=1) +
  geom_point(data = mean_lowpv1baseline, aes(x = summ_visit_month, y = median_logpv1_serum_IgA), shape = 17, size = 2, color = "#BC3908", inherit.aes = FALSE) +
  scale_y_continuous(breaks = c(0:4), limits = c(-1,4)) +
  scale_x_continuous(breaks = c(0, 1, 3, 6), limits = c(0,7.4), labels=visit_labels) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
ggsave("Figures/Serum_PV1_IgA_lowbaseline_all_medians.png", width = 6.3, height = 8, dpi = 300)
s1igalb 



#Serum PV1 IgG low baseline
s1igglb<- ggplot(lowpv1baseline, aes(x =follow_up_months, y = logPV1_IgG_serum, group = pid)) +
  geom_line(alpha = 0.5, size = 0.2) +
  geom_point(alpha = 0.1) +
  geom_ribbon(data = mean_lowpv1baseline, aes(x = summ_visit_month, ymin = lq_logpv1_serum_IgG, ymax = uq_logpv1_serum_IgG), fill = "#BC3908", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = mean_lowpv1baseline, aes(x = summ_visit_month, y = median_logpv1_serum_IgG), color = "#BC3908", inherit.aes = FALSE, size=1) +
  geom_point(data = mean_lowpv1baseline, aes(x = summ_visit_month, y = median_logpv1_serum_IgG), shape = 17, size = 2, color = "#BC3908", inherit.aes = FALSE) +
  scale_y_continuous(breaks = c(0:4), limits = c(-1,4)) +
  scale_x_continuous(breaks = c(0, 1, 3, 6), limits = c(0,7.4), labels=visit_labels) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 18),
        axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
ggsave("Figures/Serum_PV1_IgG_lowbaseline_all_medians.png", width = 6.3, height = 9, dpi = 300)
s1igglb


#Serum PV3 neuts high baseline
s3nhb<- ggplot(highpv3baseline, aes(x =follow_up_months, y = logPV3_serum, group = pid)) +
  geom_line(alpha = 0.5, size = 0.2) +
  geom_point(alpha = 0.1) +
  geom_hline(yintercept = 3, linetype = "dashed")+
  geom_ribbon(data = mean_highpv3baseline, aes(x = summ_visit_month, ymin = lq_logpv3_serum, ymax = uq_logpv3_serum), fill = "#BC3908", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = mean_highpv3baseline, aes(x = summ_visit_month, y = median_logpv3_serum), color = "#BC3908", inherit.aes = FALSE, size=1) +
  geom_point(data = mean_highpv3baseline, aes(x = summ_visit_month, y = median_logpv3_serum), shape = 17, size = 2, color = "#BC3908", inherit.aes = FALSE) +
  scale_y_continuous(breaks = c(1, 3, 5, 7, 9, 11, 13), limits = c(1,14)) +
  scale_x_continuous(breaks = c(0, 1, 3, 6), limits = c(0,7.4), labels=visit_labels) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 18),
        #axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
        axis.text.x = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
ggsave("Figures/Serum_PV3_Neuts_highbaseline_all_medians.png", width = 6.3, height = 4, dpi = 300)
s3nhb

#Serum PV3 IgA high baseline
s3igahb<- ggplot(highpv3baseline, aes(x =follow_up_months, y = logPV3_IgA_serum, group = pid)) +
  geom_line(alpha = 0.5, size = 0.2) +
  geom_point(alpha = 0.1) +
  geom_ribbon(data = mean_highpv3baseline, aes(x = summ_visit_month, ymin = lq_logpv3_serum_IgA, ymax = uq_logpv3_serum_IgA), fill = "#BC3908", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = mean_highpv3baseline, aes(x = summ_visit_month, y = median_logpv3_serum_IgA), color = "#BC3908", inherit.aes = FALSE, size=1) +
  geom_point(data = mean_highpv3baseline, aes(x = summ_visit_month, y = median_logpv3_serum_IgA), shape = 17, size = 2, color = "#BC3908", inherit.aes = FALSE) +
  scale_y_continuous(breaks = c(0:4), limits = c(-1,4)) +
  scale_x_continuous(breaks = c(0, 1, 3, 6), limits = c(0,7.4), labels=visit_labels) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
ggsave("Figures/Serum_PV3_IgA_highbaseline_all_medians.png", width = 6.3, height = 8, dpi = 300)
s3igahb 


#Serum PV3 IgG high baseline
s3igghb<- ggplot(highpv3baseline, aes(x =follow_up_months, y = logPV3_IgG_serum, group = pid)) +
  geom_line(alpha = 0.5, size = 0.2) +
  geom_point(alpha = 0.1) +
  geom_ribbon(data = mean_highpv3baseline, aes(x = summ_visit_month, ymin = lq_logpv3_serum_IgG, ymax = uq_logpv3_serum_IgG), fill = "#BC3908", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = mean_highpv3baseline, aes(x = summ_visit_month, y = median_logpv3_serum_IgG), color = "#BC3908", inherit.aes = FALSE, size=1) +
  geom_point(data = mean_highpv3baseline, aes(x = summ_visit_month, y = median_logpv3_serum_IgG), shape = 17, size = 2, color = "#BC3908", inherit.aes = FALSE) +
  scale_y_continuous(breaks = c(0:4), limits = c(-1,4)) +
  scale_x_continuous(breaks = c(0, 1, 3, 6), limits = c(0,7.4), labels=visit_labels) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 18),
        axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
ggsave("Figures/Serum_PV3_IgG_highbaseline_all_medians.png", width = 6.3, height = 9, dpi = 300)
s3igghb

#Serum PV3 neuts low baseline
s3nlb<- ggplot(lowpv3baseline, aes(x =follow_up_months, y = logPV3_serum, group = pid)) +
  geom_line(alpha = 0.5, size = 0.2) +
  geom_point(alpha = 0.1) +
  geom_hline(yintercept = 3, linetype = "dashed")+
  geom_ribbon(data = mean_lowpv3baseline, aes(x = summ_visit_month, ymin = lq_logpv3_serum, ymax = uq_logpv3_serum), fill = "#BC3908", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = mean_lowpv3baseline, aes(x = summ_visit_month, y = median_logpv3_serum), color = "#BC3908", inherit.aes = FALSE, size=1) +
  geom_point(data = mean_lowpv3baseline, aes(x = summ_visit_month, y = median_logpv3_serum), shape = 17, size = 2, color = "#BC3908", inherit.aes = FALSE) +
  scale_y_continuous(breaks = c(1, 3, 5, 7, 9, 11, 13), limits = c(1,14)) +
  scale_x_continuous(breaks = c(0, 1, 3, 6), limits = c(0,7.4), labels=visit_labels) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 18),
        #axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
        axis.text.x = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
ggsave("Figures/Serum_PV3_Neuts_lowbaseline_all_medians.png", width = 6.3, height = 4, dpi = 300)
s3nlb

#Serum PV3 IgA low baseline
s3igalb<- ggplot(lowpv3baseline, aes(x =follow_up_months, y = logPV3_IgA_serum, group = pid)) +
  geom_line(alpha = 0.5, size = 0.2) +
  geom_point(alpha = 0.1) +
  geom_ribbon(data = mean_lowpv3baseline, aes(x = summ_visit_month, ymin = lq_logpv3_serum_IgA, ymax = uq_logpv3_serum_IgA), fill = "#BC3908", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = mean_lowpv3baseline, aes(x = summ_visit_month, y = median_logpv3_serum_IgA), color = "#BC3908", inherit.aes = FALSE, size=1) +
  geom_point(data = mean_lowpv3baseline, aes(x = summ_visit_month, y = median_logpv3_serum_IgA), shape = 17, size = 2, color = "#BC3908", inherit.aes = FALSE) +
  scale_y_continuous(breaks = c(0:4), limits = c(-1,4)) +
  scale_x_continuous(breaks = c(0, 1, 3, 6), limits = c(0,7.4), labels=visit_labels) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
ggsave("Figures/Serum_PV3_IgA_lowbaseline_all_medians.png", width = 6.3, height = 8, dpi = 300)
s3igalb 

#Serum PV3 IgG low baseline
s3igglb<- ggplot(lowpv3baseline, aes(x =follow_up_months, y = logPV3_IgG_serum, group = pid)) +
  geom_line(alpha = 0.5, size = 0.2) +
  geom_point(alpha = 0.1) +
  geom_ribbon(data = mean_lowpv3baseline, aes(x = summ_visit_month, ymin = lq_logpv3_serum_IgG, ymax = uq_logpv3_serum_IgG), fill = "#BC3908", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = mean_lowpv3baseline, aes(x = summ_visit_month, y = median_logpv3_serum_IgG), color = "#BC3908", inherit.aes = FALSE, size=1) +
  geom_point(data = mean_lowpv3baseline, aes(x = summ_visit_month, y = median_logpv3_serum_IgG), shape = 17, size = 2, color = "#BC3908", inherit.aes = FALSE) +
  scale_y_continuous(breaks = c(0:4), limits = c(-1,4)) +
  scale_x_continuous(breaks = c(0, 1, 3, 6), limits = c(0,7.4), labels=visit_labels) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 18),
        axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
ggsave("Figures/Serum_PV3_IgG_lowbaseline_all_medians.png", width = 6.3, height = 9, dpi = 300)
s3igglb


#####################################################################################
#####################################################################################
#####################################################################################
#Complete data:
#Stool neuts plots

#High PV1 baseline, complete data
mean_completedata_highpv1baseline <- completedata_highpv1baseline %>%
  group_by(summ_visit_month) %>%
  summarise(mean_logpv1_stool = mean(logPV1_stool, na.rm = TRUE),
            mean_logpv1_serum = mean(logPV1_serum, na.rm = TRUE),
            median_logpv1_stool = median(logPV1_stool, na.rm = TRUE),
            median_logpv1_serum = median(logPV1_serum, na.rm = TRUE),
            lq_logpv1_stool = quantile(logPV1_stool, 0.25, na.rm = TRUE),
            lq_logpv1_serum = quantile(logPV1_serum, 0.25, na.rm = TRUE),
            uq_logpv1_stool = quantile(logPV1_stool, 0.75, na.rm = TRUE),
            uq_logpv1_serum = quantile(logPV1_serum, 0.75, na.rm = TRUE))


m1nchb <-ggplot(completedata_highpv1baseline, aes(x =follow_up_months, y = logPV1_stool, group = pid)) +
  geom_line(alpha = 0.5, size=0.2) +
  geom_point(alpha = 0.1) +
  geom_hline(yintercept = 3, linetype = "dashed")+
  geom_ribbon(data = mean_completedata_highpv1baseline, aes(x = summ_visit_month, ymin = lq_logpv1_stool, ymax = uq_logpv1_stool), fill = "#725AC1", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = mean_completedata_highpv1baseline, aes(x = summ_visit_month, y = median_logpv1_stool), color = "#725AC1", inherit.aes = FALSE, size=1) +
  geom_point(data = mean_completedata_highpv1baseline, aes(x = summ_visit_month, y = median_logpv1_stool), shape = 17, size = 2, color = "#725AC1", inherit.aes = FALSE) +
  scale_y_continuous(breaks = c(1, 3, 5, 7, 9, 11, 13), limits = c(1,14)) +
  scale_x_continuous(breaks = c(0, 1, 3, 6), limits = c(0,7.4), labels=visit_labels) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 18),
        #axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
        axis.text.x = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
m1nchb
ggsave("Figures/Complete data/Stool_PV1_Neuts_highbaseline_completedata_medians.png", width = 6.3, height = 5.6, dpi = 300)

#Low PV1 baseline, complete data
mean_completedata_lowpv1baseline <- completedata_lowpv1baseline %>%
  group_by(summ_visit_month) %>%
  summarise(mean_logpv1_stool = mean(logPV1_stool, na.rm = TRUE),
            mean_logpv1_serum = mean(logPV1_serum, na.rm = TRUE),
            median_logpv1_stool = median(logPV1_stool, na.rm = TRUE),
            median_logpv1_serum = median(logPV1_serum, na.rm = TRUE),
            lq_logpv1_stool = quantile(logPV1_stool, 0.25, na.rm = TRUE),
            lq_logpv1_serum = quantile(logPV1_serum, 0.25, na.rm = TRUE),
            uq_logpv1_stool = quantile(logPV1_stool, 0.75, na.rm = TRUE),
            uq_logpv1_serum = quantile(logPV1_serum, 0.75, na.rm = TRUE))


m1nclb <-ggplot(completedata_lowpv1baseline, aes(x =follow_up_months, y = logPV1_stool, group = pid)) +
  geom_line(alpha = 0.5, size=0.2) +
  geom_point(alpha = 0.1) +
  geom_hline(yintercept = 3, linetype = "dashed")+
  geom_ribbon(data = mean_completedata_lowpv1baseline, aes(x = summ_visit_month, ymin = lq_logpv1_stool, ymax = uq_logpv1_stool), fill = "#725AC1", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = mean_completedata_lowpv1baseline, aes(x = summ_visit_month, y = median_logpv1_stool), color = "#725AC1", inherit.aes = FALSE, size=1) +
  geom_point(data = mean_completedata_lowpv1baseline, aes(x = summ_visit_month, y = median_logpv1_stool), shape = 17, size = 2, color = "#725AC1", inherit.aes = FALSE) +
  scale_y_continuous(breaks = c(1, 3, 5, 7, 9, 11, 13), limits = c(1,14)) +
  scale_x_continuous(breaks = c(0, 1, 3, 6), limits = c(0,7.4), labels=visit_labels) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 18),
        #axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
        axis.text.x = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
m1nclb
ggsave("Figures/Complete data/Stool_PV1_Neuts_lowbaseline_completedata_medians.png", width = 6.3, height = 5.6, dpi = 300)


#High PV3 baseline, complete data
mean_completedata_highpv3baseline <- completedata_highpv3baseline %>%
  group_by(summ_visit_month) %>%
  summarise(mean_logpv3_stool = mean(logPV3_stool, na.rm = TRUE),
            mean_logpv3_serum = mean(logPV3_serum, na.rm = TRUE),    
            median_logpv3_stool = median(logPV3_stool, na.rm = TRUE),
            median_logpv3_serum = median(logPV3_serum, na.rm = TRUE),
            lq_logpv3_stool = quantile(logPV3_stool, 0.25, na.rm = TRUE),
            lq_logpv3_serum = quantile(logPV3_serum, 0.25, na.rm = TRUE),
            uq_logpv3_stool = quantile(logPV3_stool, 0.75, na.rm = TRUE),
            uq_logpv3_serum = quantile(logPV3_serum, 0.75, na.rm = TRUE))

m3nchb <-ggplot(completedata_highpv3baseline , aes(x =follow_up_months, y = logPV3_stool, group = pid)) +
  geom_line(alpha = 0.5, size=0.2) +
  geom_point(alpha = 0.1) +
  geom_hline(yintercept = 3, linetype = "dashed")+
  geom_ribbon(data = mean_completedata_highpv3baseline, aes(x = summ_visit_month, ymin = lq_logpv3_stool, ymax = uq_logpv3_stool), fill = "#725AC1", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = mean_completedata_highpv3baseline, aes(x = summ_visit_month, y = median_logpv3_stool), color = "#725AC1", inherit.aes = FALSE, size=1) +
  geom_point(data = mean_completedata_highpv3baseline, aes(x = summ_visit_month, y = median_logpv3_stool), shape = 17, size = 2, color = "#725AC1", inherit.aes = FALSE) +
  scale_y_continuous(breaks = c(1, 3, 5, 7, 9, 11, 13), limits = c(1,14)) +
  scale_x_continuous(breaks = c(0, 1, 3, 6), limits = c(0,7.4), labels=visit_labels) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 18),
        #axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
        axis.text.x = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
ggsave("Figures/Complete data/Stool_PV3_Neuts_highbaseline_completedata_medians.png", width = 6.3, height = 5.6, dpi = 300)
m3nchb

#Low PV3 baseline, complete data
mean_completedata_lowpv3baseline <- completedata_lowpv3baseline %>%
  group_by(summ_visit_month) %>%
  summarise(mean_logpv3_stool = mean(logPV3_stool, na.rm = TRUE),
            mean_logpv3_serum = mean(logPV3_serum, na.rm = TRUE),    
            median_logpv3_stool = median(logPV3_stool, na.rm = TRUE),
            median_logpv3_serum = median(logPV3_serum, na.rm = TRUE),
            lq_logpv3_stool = quantile(logPV3_stool, 0.25, na.rm = TRUE),
            lq_logpv3_serum = quantile(logPV3_serum, 0.25, na.rm = TRUE),
            uq_logpv3_stool = quantile(logPV3_stool, 0.75, na.rm = TRUE),
            uq_logpv3_serum = quantile(logPV3_serum, 0.75, na.rm = TRUE))


m3nclb <-ggplot(completedata_lowpv3baseline , aes(x =follow_up_months, y = logPV3_stool, group = pid)) +
  geom_line(alpha = 0.5, size=0.2) +
  geom_point(alpha = 0.1) +
  geom_hline(yintercept = 3, linetype = "dashed")+
  geom_ribbon(data = mean_completedata_lowpv3baseline, aes(x = summ_visit_month, ymin = lq_logpv3_stool, ymax = uq_logpv3_stool), fill = "#725AC1", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = mean_completedata_lowpv3baseline, aes(x = summ_visit_month, y = median_logpv3_stool), color = "#725AC1", inherit.aes = FALSE, size=1) +
  geom_point(data = mean_completedata_lowpv3baseline, aes(x = summ_visit_month, y = median_logpv3_stool), shape = 17, size = 2, color = "#725AC1", inherit.aes = FALSE) +
  scale_y_continuous(breaks = c(1, 3, 5, 7, 9, 11, 13), limits = c(1,14)) +
  scale_x_continuous(breaks = c(0, 1, 3, 6), limits = c(0,7.4), labels=visit_labels) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 18),
        #axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
        axis.text.x = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
ggsave("Figures/Complete data/Stool_PV3_Neuts_lowbaseline_completedata_medians.png", width = 6.3, height = 5.6, dpi = 300)
m3nclb

#serum neuts plots for complete data stratified by mucosal neut baseline

s1nchb<- ggplot(completedata_highpv1baseline, aes(x =follow_up_months, y = logPV1_serum, group = pid)) +
  geom_line(alpha = 0.5, size = 0.2) +
  geom_point(alpha = 0.1) +
  geom_hline(yintercept = 3, linetype = "dashed")+
  geom_ribbon(data = mean_completedata_highpv1baseline, aes(x = summ_visit_month, ymin = lq_logpv1_serum, ymax = uq_logpv1_serum), fill = "#BC3908", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = mean_completedata_highpv1baseline, aes(x = summ_visit_month, y = median_logpv1_serum), color = "#BC3908", inherit.aes = FALSE, size=1) +
  geom_point(data = mean_completedata_highpv1baseline, aes(x = summ_visit_month, y = median_logpv1_serum), shape = 17, size = 2, color = "#BC3908", inherit.aes = FALSE) +
  scale_y_continuous(breaks = c(1, 3, 5, 7, 9, 11, 13), limits = c(1,14)) +
  scale_x_continuous(breaks = c(0, 1, 3, 6), limits = c(0,7.4), labels=visit_labels) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 18),
        axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
ggsave("Figures/Complete data/Serum_PV1_Neuts_completedata_highbaseline_all_medians.png", width = 6.3, height = 5.6, dpi = 300)
s1nchb 

s1nclb<- ggplot(completedata_lowpv1baseline, aes(x =follow_up_months, y = logPV1_serum, group = pid)) +
  geom_line(alpha = 0.5, size = 0.2) +
  geom_point(alpha = 0.1) +
  geom_hline(yintercept = 3, linetype = "dashed")+
  geom_ribbon(data = mean_completedata_lowpv1baseline, aes(x = summ_visit_month, ymin = lq_logpv1_serum, ymax = uq_logpv1_serum), fill = "#BC3908", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = mean_completedata_lowpv1baseline, aes(x = summ_visit_month, y = median_logpv1_serum), color = "#BC3908", inherit.aes = FALSE, size=1) +
  geom_point(data = mean_completedata_lowpv1baseline, aes(x = summ_visit_month, y = median_logpv1_serum), shape = 17, size = 2, color = "#BC3908", inherit.aes = FALSE) +
  scale_y_continuous(breaks = c(1, 3, 5, 7, 9, 11, 13), limits = c(1,14)) +
  scale_x_continuous(breaks = c(0, 1, 3, 6), limits = c(0,7.4), labels=visit_labels) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 18),
        axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
ggsave("Figures/Complete data/Serum_PV1_Neuts_completedata_lowbaseline_all_medians.png", width = 6.3, height = 5.6, dpi = 300)
s1nclb 

s3nchb<- ggplot(completedata_highpv3baseline, aes(x =follow_up_months, y = logPV3_serum, group = pid)) +
  geom_line(alpha = 0.5, size = 0.2) +
  geom_point(alpha = 0.1) +
  geom_hline(yintercept = 3, linetype = "dashed")+
  geom_ribbon(data = mean_completedata_highpv3baseline, aes(x = summ_visit_month, ymin = lq_logpv3_serum, ymax = uq_logpv3_serum), fill = "#BC3908", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = mean_completedata_highpv3baseline, aes(x = summ_visit_month, y = median_logpv3_serum), color = "#BC3908", inherit.aes = FALSE, size=1) +
  geom_point(data = mean_completedata_highpv3baseline, aes(x = summ_visit_month, y = median_logpv3_serum), shape = 17, size = 2, color = "#BC3908", inherit.aes = FALSE) +
  scale_y_continuous(breaks = c(1, 3, 5, 7, 9, 11, 13), limits = c(1,14)) +
  scale_x_continuous(breaks = c(0, 1, 3, 6), limits = c(0,7.4), labels=visit_labels) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 18),
        axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
ggsave("Figures/Complete data/Serum_PV3_Neuts_completedata_highbaseline_all_medians.png", width = 6.3, height = 5.6, dpi = 300)
s3nchb 


s3nclb<- ggplot(completedata_lowpv3baseline, aes(x =follow_up_months, y = logPV3_serum, group = pid)) +
  geom_line(alpha = 0.5, size = 0.2) +
  geom_point(alpha = 0.1) +
  geom_hline(yintercept = 3, linetype = "dashed")+
  geom_ribbon(data = mean_completedata_lowpv3baseline, aes(x = summ_visit_month, ymin = lq_logpv3_serum, ymax = uq_logpv3_serum), fill = "#BC3908", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = mean_completedata_lowpv3baseline, aes(x = summ_visit_month, y = median_logpv3_serum), color = "#BC3908", inherit.aes = FALSE, size=1) +
  geom_point(data = mean_completedata_lowpv3baseline, aes(x = summ_visit_month, y = median_logpv3_serum), shape = 17, size = 2, color = "#BC3908", inherit.aes = FALSE) +
  scale_y_continuous(breaks = c(1, 3, 5, 7, 9, 11, 13), limits = c(1,14)) +
  scale_x_continuous(breaks = c(0, 1, 3, 6), limits = c(0,7.4), labels=visit_labels) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 18),
        axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
ggsave("Figures/Complete data/Serum_PV3_Neuts_completedata_lowbaseline_all_medians.png", width = 6.3, height = 5.6, dpi = 300)
s3nclb 

