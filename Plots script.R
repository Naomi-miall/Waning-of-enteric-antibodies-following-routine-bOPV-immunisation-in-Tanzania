###R script creating scatter plots for TZA duration of immunity paper
#Naomi and Alex
#October 2025

#Set up: Libaries; working directory; theme:
### scatter plot ####
setwd ("C:/Users/NaomiMiall/Documents/Polio work/Tanzania")

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
#vis <- subset (longitudinal, visit=="bl"|visit=="w2"|visit=="w4"|visit=="w6"|visit=="w8")
#vis_ser <- subset (longitudinal, vis=="bl"|vis=="w4"|vis=="w8")

longitudinal$visit <- factor(longitudinal$visit, levels = c("1", "2", "3", "4"),
                             labels = c("Baseline", "One month", "3 months", "6 months"))

#Remove the two age outlier participants
longitudinal <- filter(longitudinal, flag_ageoutlier==0)

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

# m1n<- ggplot(longitudinal, aes(x =follow_up_days, y = logPV1_stool, group = pid)) +
#   geom_line(alpha = 0.5, size = 0.2) +
#   geom_point(alpha = 0.1) +
#   geom_hline(yintercept = 3, linetype = "dashed")+
#   geom_point(data = mean_longitudinal, aes(x = summ_visit_day, y = mean_logpv1_stool), shape = 17, size = 2, color = "red", inherit.aes = FALSE) +
#   geom_line(data = mean_longitudinal, aes(x = summ_visit_day, y = mean_logpv1_stool), color = "red", inherit.aes = FALSE) +
#   scale_y_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12))
# m1n
# 
# m3n<- ggplot(longitudinal, aes(x =follow_up_days, y = logPV3_stool, group = pid)) +
#   geom_line(alpha = 0.5, size = 0.2) +
#   geom_point(alpha = 0.1) +
#   geom_hline(yintercept = 3, linetype = "dashed")+
#   geom_point(data = mean_longitudinal, aes(x = summ_visit_day, y = mean_logpv3_stool), shape = 17, size = 2, color = "red", inherit.aes = FALSE) +
#   geom_line(data = mean_longitudinal, aes(x = summ_visit_day, y = mean_logpv3_stool), color = "red", inherit.aes = FALSE) +
#   scale_y_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12))
# m3n
# 
# m2n <- ggplot(longitudinal, aes(x =follow_up_days, y = logPV2_stool, group = pid)) +
#   geom_line(alpha = 0.5, size = 0.2) +
#   geom_point(alpha = 0.1) +
#   geom_hline(yintercept = 3, linetype = "dashed")+
#   geom_point(data = mean_longitudinal, aes(x = summ_visit_day, y = mean_logpv2_stool), shape = 17, size = 2, color = "red", inherit.aes = FALSE) +
#   geom_line(data = mean_longitudinal, aes(x = summ_visit_day, y = mean_logpv2_stool), color = "red", inherit.aes = FALSE) +
#   scale_y_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12))
# ggsave("Figures/Stool_PV2_all.png", width = 6, height = 4, dpi = 300)
# m2n

m2n <- ggplot(longitudinal, aes(x =follow_up_months, y = logPV2_stool, group = pid)) +
    geom_line(alpha = 0.5, size = 0.2) +
    geom_point(alpha = 0.1) +
    geom_hline(yintercept = 3, linetype = "dashed")+
  geom_ribbon(data = mean_longitudinal, aes(x = summ_visit_month, ymin = lq_logpv2_stool, ymax = uq_logpv2_stool), fill = "red", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = mean_longitudinal, aes(x = summ_visit_month, y = median_logpv2_stool), color = "red", inherit.aes = FALSE) +
  geom_point(data = mean_longitudinal, aes(x = summ_visit_month, y = median_logpv2_stool), shape = 17, size = 2, color = "red", inherit.aes = FALSE) +
  scale_y_continuous(breaks = 1:10, limits= c(1,10)) +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text = element_text(size = 18))
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
            uq_logpv1_serum = quantile(logPV1_serum, 0.75, na.rm = TRUE))

# m1nhb<- ggplot(highpv1baseline, aes(x =follow_up_days, y = logPV1_stool, group = pid)) +
#   geom_line(alpha = 0.5, size = 0.2) +
#   geom_point(alpha = 0.1) +
#   geom_hline(yintercept = 3, linetype = "dashed")+
#   geom_point(data = mean_highpv1baseline, aes(x = summ_visit_day, y = mean_logpv1_stool), shape = 17, size = 2, color = "red", inherit.aes = FALSE) +
#   geom_line(data = mean_highpv1baseline, aes(x = summ_visit_day, y = mean_logpv1_stool), color = "red", inherit.aes = FALSE) +
#   scale_y_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12))
# ggsave("Figures/Stool_PV1_Neuts_highbaseline_all.png", width = 6, height = 4, dpi = 300)
# m1nhb

m1nhb <- ggplot(highpv1baseline, aes(x = follow_up_months, y = logPV1_stool, group = pid)) +
  geom_line(alpha = 0.5, size = 0.2) + geom_point(alpha = 0.1) +
  geom_hline(yintercept = 3, linetype = "dashed") +
  geom_ribbon(data = mean_highpv1baseline, aes(x = summ_visit_month, ymin = lq_logpv1_stool, ymax = uq_logpv1_stool), fill = "red", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = mean_highpv1baseline, aes(x = summ_visit_month, y = median_logpv1_stool), color = "red", inherit.aes = FALSE) +
  geom_point(data = mean_highpv1baseline, aes(x = summ_visit_month, y = median_logpv1_stool), shape = 17, size = 2, color = "red", inherit.aes = FALSE) +
  scale_y_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)) +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text = element_text(size = 18))
ggsave("Figures/Stool_PV1_Neuts_highbaseline_all_medians.png", width = 6, height = 4, dpi = 300)
m1nhb

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
            uq_logpv1_serum = quantile(logPV1_serum, 0.75, na.rm = TRUE))

# m1nlb <- ggplot(lowpv1baseline, aes(x =follow_up_days, y = logPV1_stool, group = pid)) +
#   geom_line(alpha = 0.5, size = 0.2) +
#   geom_point(alpha = 0.1) +
#   geom_hline(yintercept = 3, linetype = "dashed")+
#   geom_point(data = mean_lowpv1baseline, aes(x = summ_visit_day, y = mean_logpv1_stool), shape = 17, size = 2, color = "red", inherit.aes = FALSE) +
#   geom_line(data = mean_lowpv1baseline, aes(x = summ_visit_day, y = mean_logpv1_stool), color = "red", inherit.aes = FALSE) +
#   scale_y_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12))
# ggsave("Figures/Stool_PV1_Neuts_lowbaseline_all.png", width = 6, height = 4, dpi = 300)
# m1nlb

m1nlb <- ggplot(lowpv1baseline, aes(x = follow_up_months, y = logPV1_stool, group = pid)) +
  geom_line(alpha = 0.5, size = 0.2) + geom_point(alpha = 0.1) +
  geom_hline(yintercept = 3, linetype = "dashed") +
  geom_ribbon(data = mean_lowpv1baseline, aes(x = summ_visit_month, ymin = lq_logpv1_stool, ymax = uq_logpv1_stool), fill = "red", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = mean_lowpv1baseline, aes(x = summ_visit_month, y = median_logpv1_stool), color = "red", inherit.aes = FALSE) +
  geom_point(data = mean_lowpv1baseline, aes(x = summ_visit_month, y = median_logpv1_stool), shape = 17, size = 2, color = "red", inherit.aes = FALSE) +
  scale_y_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text = element_text(size = 18))
ggsave("Figures/Stool_PV1_Neuts_lowbaseline_all_medians.png", width = 6, height = 4, dpi = 300)
m1nlb

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
            uq_logpv3_serum = quantile(logPV3_serum, 0.75, na.rm = TRUE))

# m3nhb <- ggplot(highpv3baseline, aes(x =follow_up_days, y = logPV3_stool, group = pid)) +
#   geom_line(alpha = 0.5, size = 0.2) +
#   geom_point(alpha = 0.1) +
#   geom_hline(yintercept = 3, linetype = "dashed")+
#   geom_point(data = mean_highpv3baseline, aes(x = summ_visit_day, y = mean_logpv3_stool), shape = 17, size = 2, color = "red", inherit.aes = FALSE) +
#   geom_line(data = mean_highpv3baseline, aes(x = summ_visit_day, y = mean_logpv3_stool), color = "red", inherit.aes = FALSE) +
#   scale_y_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12))
# ggsave("Figures/Stool_PV3_Neuts_highbaseline_all.png", width = 6, height = 4, dpi = 300)
# m3nhb 

m3nhb <- ggplot(highpv3baseline, aes(x = follow_up_months, y = logPV3_stool, group = pid)) +
  geom_line(alpha = 0.5, size = 0.2) + geom_point(alpha = 0.1) +
  geom_hline(yintercept = 3, linetype = "dashed") +
  geom_ribbon(data = mean_highpv3baseline, aes(x = summ_visit_month, ymin = lq_logpv3_stool, ymax = uq_logpv3_stool), fill = "red", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = mean_highpv3baseline, aes(x = summ_visit_month, y = median_logpv3_stool), color = "red", inherit.aes = FALSE) +
  geom_point(data = mean_highpv3baseline, aes(x = summ_visit_month, y = median_logpv3_stool), shape = 17, size = 2, color = "red", inherit.aes = FALSE) +
  scale_y_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text = element_text(size = 18))
ggsave("Figures/Stool_PV3_Neuts_highbaseline_all_medians.png", width = 6, height = 4, dpi = 300)
m3nhb 

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
            uq_logpv3_serum = quantile(logPV3_serum, 0.75, na.rm = TRUE))

# m3nlb <- ggplot(lowpv3baseline, aes(x =follow_up_days, y = logPV3_stool, group = pid)) +
#   geom_line(alpha = 0.5, size = 0.2) +
#   geom_point(alpha = 0.1) +
#   geom_hline(yintercept = 3, linetype = "dashed")+
#   geom_point(data = mean_lowpv3baseline, aes(x = summ_visit_day, y = mean_logpv3_stool), shape = 17, size = 2, color = "red", inherit.aes = FALSE) +
#   geom_line(data = mean_lowpv3baseline, aes(x = summ_visit_day, y = mean_logpv3_stool), color = "red", inherit.aes = FALSE) +
#   scale_y_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12))
# ggsave("Figures/Stool_PV3_Neuts_lowbaseline_all.png", width = 6, height = 4, dpi = 300)
# m3nlb

m3nlb <- ggplot(lowpv3baseline, aes(x =follow_up_months, y = logPV3_stool, group = pid)) +
  geom_line(alpha = 0.5, size = 0.2) +
  geom_point(alpha = 0.1) +
  geom_hline(yintercept = 3, linetype = "dashed")+
  geom_ribbon(data = mean_lowpv3baseline, aes(x = summ_visit_month, ymin = lq_logpv3_stool, ymax = uq_logpv3_stool), fill = "red", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = mean_lowpv3baseline, aes(x = summ_visit_month, y = median_logpv3_stool), color = "red", inherit.aes = FALSE) +
  geom_point(data = mean_lowpv3baseline, aes(x = summ_visit_month, y = median_logpv3_stool), shape = 17, size = 2, color = "red", inherit.aes = FALSE) +
  scale_y_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text = element_text(size = 18))
ggsave("Figures/Stool_PV3_Neuts_lowbaseline_all_medians.png", width = 6, height = 4, dpi = 300)
m3nlb 

#####################################################################################
#Stool neut plots for complete data:

#   #Complete data
# mean_completedata <- completedata %>%
#   group_by(summ_visit_month) %>%
#   summarise(mean_logpv1_stool = mean(logPV1_stool, na.rm = TRUE),
#             mean_logpv3_stool = mean(logPV3_stool, na.rm = TRUE),
#             mean_logpv1_serum = mean(logPV1_serum, na.rm = TRUE),
#             mean_logpv3_serum = mean(logPV3_serum, na.rm = TRUE))

# m1nc <-ggplot(completedata, aes(x =follow_up_days, y = logPV1_stool, group = pid)) +
#   geom_line(alpha = 0.5, size=0.2) +
#   geom_point(alpha = 0.1) +
#   geom_hline(yintercept = 3, linetype = "dashed")+
#   geom_point(data = mean_completedata, aes(x = summ_visit_day, y = mean_logpv1_stool), shape = 17, size = 2, color = "red", inherit.aes = FALSE) +
#   geom_line(data = mean_completedata, aes(x = summ_visit_day, y = mean_logpv1_stool), color = "red", inherit.aes = FALSE) +
#   scale_y_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12))
# m1nc
# 
# m3nc <-ggplot(completedata, aes(x =follow_up_days, y = logPV3_stool, group = pid)) +
#   geom_line(alpha = 0.5, size=0.2) +
#   geom_point(alpha = 0.1) +
#   geom_hline(yintercept = 3, linetype = "dashed")+
#   geom_point(data = mean_completedata, aes(x = summ_visit_day, y = mean_logpv3_stool), shape = 17, size = 2, color = "red", inherit.aes = FALSE) +
#   geom_line(data = mean_completedata, aes(x = summ_visit_day, y = mean_logpv3_stool), color = "red", inherit.aes = FALSE) +
#   scale_y_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12))
# m3nc

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

# m1nchb <-ggplot(completedata_highpv1baseline, aes(x =follow_up_days, y = logPV1_stool, group = pid)) +
#   geom_line(alpha = 0.5, size=0.2) +
#   geom_point(alpha = 0.1) +
#   geom_hline(yintercept = 3, linetype = "dashed")+
#   geom_point(data = mean_completedata_highpv1baseline, aes(x = summ_visit_day, y = mean_logpv1_stool), shape = 17, size = 2, color = "red", inherit.aes = FALSE) +
#   geom_line(data = mean_completedata_highpv1baseline, aes(x = summ_visit_day, y = mean_logpv1_stool), color = "red", inherit.aes = FALSE) +
#   scale_y_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12))
# ggsave("Figures/Stool_PV1_Neuts_highbaseline_completedata.png", width = 6.3, height = 5.6, dpi = 300)
# m1nchb 

m1nchb <-ggplot(completedata_highpv1baseline, aes(x =follow_up_months, y = logPV1_stool, group = pid)) +
  geom_line(alpha = 0.5, size=0.2) +
  geom_point(alpha = 0.1) +
  geom_hline(yintercept = 3, linetype = "dashed")+
  geom_ribbon(data = mean_completedata_highpv1baseline, aes(x = summ_visit_month, ymin = lq_logpv1_stool, ymax = uq_logpv1_stool), fill = "red", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = mean_completedata_highpv1baseline, aes(x = summ_visit_month, y = median_logpv1_stool), color = "red", inherit.aes = FALSE) +
  geom_point(data = mean_completedata_highpv1baseline, aes(x = summ_visit_month, y = median_logpv1_stool), shape = 17, size = 2, color = "red", inherit.aes = FALSE) +
  scale_y_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text = element_text(size = 18))
ggsave("Figures/Stool_PV1_Neuts_highbaseline_completedata_medians.png", width = 6.3, height = 5.6, dpi = 300)
m1nchb 

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

# m1nclb <-ggplot(completedata_lowpv1baseline, aes(x =follow_up_days, y = logPV1_stool, group = pid)) +
#   geom_line(alpha = 0.5, size=0.2) +
#   geom_point(alpha = 0.1) +
#   geom_hline(yintercept = 3, linetype = "dashed")+
#   geom_point(data = mean_completedata_lowpv1baseline, aes(x = summ_visit_day, y = mean_logpv1_stool), shape = 17, size = 2, color = "red", inherit.aes = FALSE) +
#   geom_line(data = mean_completedata_lowpv1baseline, aes(x = summ_visit_day, y = mean_logpv1_stool), color = "red", inherit.aes = FALSE) +
#   scale_y_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12))
# ggsave("Figures/Stool_PV1_Neuts_lowbaseline_completedata.png", width = 6.3, height = 5.6, dpi = 300)
# m1nclb

m1nclb <-ggplot(completedata_lowpv1baseline, aes(x =follow_up_months, y = logPV1_stool, group = pid)) +
  geom_line(alpha = 0.5, size=0.2) +
  geom_point(alpha = 0.1) +
  geom_hline(yintercept = 3, linetype = "dashed")+
  geom_ribbon(data = mean_completedata_lowpv1baseline, aes(x = summ_visit_month, ymin = lq_logpv1_stool, ymax = uq_logpv1_stool), fill = "red", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = mean_completedata_lowpv1baseline, aes(x = summ_visit_month, y = median_logpv1_stool), color = "red", inherit.aes = FALSE) +
  geom_point(data = mean_completedata_lowpv1baseline, aes(x = summ_visit_month, y = median_logpv1_stool), shape = 17, size = 2, color = "red", inherit.aes = FALSE) +
  scale_y_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text = element_text(size = 18))
ggsave("Figures/Stool_PV1_Neuts_lowbaseline_completedata_medians.png", width = 6.3, height = 5.6, dpi = 300)
m1nclb


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

# m3nchb <-ggplot(completedata_highpv3baseline , aes(x =follow_up_days, y = logPV3_stool, group = pid)) +
#   geom_line(alpha = 0.5, size=0.2) +
#   geom_point(alpha = 0.1) +
#   geom_hline(yintercept = 3, linetype = "dashed")+
#   geom_point(data = mean_completedata_highpv3baseline , aes(x = summ_visit_day, y = mean_logpv3_stool), shape = 17, size = 2, color = "red", inherit.aes = FALSE) +
#   geom_line(data = mean_completedata_highpv3baseline , aes(x = summ_visit_day, y = mean_logpv3_stool), color = "red", inherit.aes = FALSE) +
#   scale_y_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12))
# ggsave("Figures/Stool_PV3_Neuts_highbaseline_completedata.png", width = 6.3, height = 5.6, dpi = 300)
# m3nchb

m3nchb <-ggplot(completedata_highpv3baseline , aes(x =follow_up_months, y = logPV3_stool, group = pid)) +
  geom_line(alpha = 0.5, size=0.2) +
  geom_point(alpha = 0.1) +
  geom_hline(yintercept = 3, linetype = "dashed")+
  geom_ribbon(data = mean_completedata_highpv3baseline, aes(x = summ_visit_month, ymin = lq_logpv3_stool, ymax = uq_logpv3_stool), fill = "red", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = mean_completedata_highpv3baseline, aes(x = summ_visit_month, y = median_logpv3_stool), color = "red", inherit.aes = FALSE) +
  geom_point(data = mean_completedata_highpv3baseline, aes(x = summ_visit_month, y = median_logpv3_stool), shape = 17, size = 2, color = "red", inherit.aes = FALSE) +
  scale_y_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text = element_text(size = 18))
ggsave("Figures/Stool_PV3_Neuts_highbaseline_completedata_medians.png", width = 6.3, height = 5.6, dpi = 300)
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

# m3nclb <-ggplot(completedata_lowpv3baseline , aes(x =follow_up_days, y = logPV3_stool, group = pid)) +
#   geom_line(alpha = 0.5, size=0.2) +
#   geom_point(alpha = 0.1) +
#   geom_hline(yintercept = 3, linetype = "dashed")+
#   geom_point(data = mean_completedata_lowpv3baseline , aes(x = summ_visit_day, y = mean_logpv3_stool), shape = 17, size = 2, color = "red", inherit.aes = FALSE) +
#   geom_line(data = mean_completedata_lowpv3baseline , aes(x = summ_visit_day, y = mean_logpv3_stool), color = "red", inherit.aes = FALSE) +
#   scale_y_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12))
# ggsave("Figures/Stool_PV3_Neuts_lowbaseline_completedata.png", width = 6.3, height = 5.6, dpi = 300)
# m3nclb

m3nclb <-ggplot(completedata_lowpv3baseline , aes(x =follow_up_months, y = logPV3_stool, group = pid)) +
  geom_line(alpha = 0.5, size=0.2) +
  geom_point(alpha = 0.1) +
  geom_hline(yintercept = 3, linetype = "dashed")+
  geom_ribbon(data = mean_completedata_lowpv3baseline, aes(x = summ_visit_month, ymin = lq_logpv3_stool, ymax = uq_logpv3_stool), fill = "red", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = mean_completedata_lowpv3baseline, aes(x = summ_visit_month, y = median_logpv3_stool), color = "red", inherit.aes = FALSE) +
  geom_point(data = mean_completedata_lowpv3baseline, aes(x = summ_visit_month, y = median_logpv3_stool), shape = 17, size = 2, color = "red", inherit.aes = FALSE) +
  scale_y_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text = element_text(size = 18))
ggsave("Figures/Stool_PV3_Neuts_lowbaseline_completedata_medians.png", width = 6.3, height = 5.6, dpi = 300)
m3nclb

#################################################################################
#################################################################################
## Serum neutr plot: Spagetti plots

#All participants unstratified for PV2
s2n <- ggplot(longitudinal, aes(x =follow_up_months, y = logPV2_serum, group = pid)) +
  geom_line(alpha = 0.5, size = 0.2) +
  geom_point(alpha = 0.1) +
  geom_ribbon(data = mean_longitudinal, aes(x = summ_visit_month, ymin = lq_logpv2_serum, ymax = uq_logpv2_serum), fill = "red", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = mean_longitudinal, aes(x = summ_visit_month, y = median_logpv2_serum), color = "red", inherit.aes = FALSE) +
  geom_point(data = mean_longitudinal, aes(x = summ_visit_month, y = median_logpv2_serum), shape = 17, size = 2, color = "red", inherit.aes = FALSE) +
  scale_y_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)) +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text = element_text(size = 18))
ggsave("Figures/Serum_PV2_all_median.png", width = 6, height = 4, dpi = 300)
s2n

#all participants neuts stratified by mucosal neut baseline
# s1nhb<- ggplot(highpv1baseline, aes(x =follow_up_days, y = logPV1_serum, group = pid)) +
#   geom_line(alpha = 0.5, size = 0.2) +
#   geom_point(alpha = 0.1) +
#   geom_point(data = mean_highpv1baseline, aes(x = summ_visit_day, y = mean_logpv1_serum), shape = 17, size = 2, color = "red", inherit.aes = FALSE) +
#   geom_line(data = mean_highpv1baseline, aes(x = summ_visit_day, y = mean_logpv1_serum), color = "red", inherit.aes = FALSE) +
#   scale_y_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12))
# s1nhb #saved as Serum_PV1_Neuts_highbaseline_all

#High baseline pv1 - all
s1nhb<- ggplot(highpv1baseline, aes(x =follow_up_months, y = logPV1_serum, group = pid)) +
  geom_line(alpha = 0.5, size = 0.2) +
  geom_point(alpha = 0.1) +
  geom_ribbon(data = mean_highpv1baseline, aes(x = summ_visit_month, ymin = lq_logpv1_serum, ymax = uq_logpv1_serum), fill = "red", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = mean_highpv1baseline, aes(x = summ_visit_month, y = median_logpv1_serum), color = "red", inherit.aes = FALSE) +
  geom_point(data = mean_highpv1baseline, aes(x = summ_visit_month, y = median_logpv1_serum), shape = 17, size = 2, color = "red", inherit.aes = FALSE) +
  scale_y_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text = element_text(size = 18))
ggsave("Figures/Serum_PV1_Neuts_highbaseline_all_medians.png", width = 6.3, height = 5.6, dpi = 300)
s1nhb 

# s1nlb <- ggplot(lowpv1baseline, aes(x =follow_up_days, y = logPV1_serum, group = pid)) +
#   geom_line(alpha = 0.5, size = 0.2) +
#   geom_point(alpha = 0.1) +
#   geom_point(data = mean_lowpv1baseline, aes(x = summ_visit_day, y = mean_logpv1_serum), shape = 17, size = 2, color = "red", inherit.aes = FALSE) +
#   geom_line(data = mean_lowpv1baseline, aes(x = summ_visit_day, y = mean_logpv1_serum), color = "red", inherit.aes = FALSE) +
#   scale_y_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12))
# s1nlb #saved as Serum_PV1_Neuts_lowbaseline_all

s1nhb<- ggplot(lowpv1baseline, aes(x =follow_up_months, y = logPV1_serum, group = pid)) +
  geom_line(alpha = 0.5, size = 0.2) +
  geom_point(alpha = 0.1) +
  geom_ribbon(data = mean_lowpv1baseline, aes(x = summ_visit_month, ymin = lq_logpv1_serum, ymax = uq_logpv1_serum), fill = "red", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = mean_lowpv1baseline, aes(x = summ_visit_month, y = median_logpv1_serum), color = "red",  inherit.aes = FALSE) +
  geom_point(data = mean_lowpv1baseline, aes(x = summ_visit_month, y = median_logpv1_serum), shape = 17, size = 2, color = "red", inherit.aes = FALSE) +
  scale_y_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text = element_text(size = 18))
ggsave("Figures/Serum_PV1_Neuts_lowbaseline_all_medians.png", width = 6.3, height = 5.6, dpi = 300)
s1nhb


# s3nhb<- ggplot(highpv3baseline, aes(x =follow_up_days, y = logPV3_serum, group = pid)) +
#   geom_line(alpha = 0.5, size = 0.2) +
#   geom_point(alpha = 0.1) +
#   geom_point(data = mean_highpv3baseline, aes(x = summ_visit_day, y = mean_logpv3_serum), shape = 17, size = 2, color = "red", inherit.aes = FALSE) +
#   geom_line(data = mean_highpv3baseline, aes(x = summ_visit_day, y = mean_logpv3_serum), color = "red", inherit.aes = FALSE) +
#   scale_y_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12))
# s3nhb #saved as Serum_PV3_Neuts_highbaseline_all

s3nhb<- ggplot(highpv3baseline, aes(x =follow_up_months, y = logPV3_serum, group = pid)) +
  geom_line(alpha = 0.5, size = 0.2) +
  geom_point(alpha = 0.1) +
  geom_ribbon(data = mean_highpv3baseline, aes(x = summ_visit_month, ymin = lq_logpv3_serum, ymax = uq_logpv3_serum), fill = "red", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = mean_highpv3baseline, aes(x = summ_visit_month, y = median_logpv3_serum), color = "red", inherit.aes = FALSE) +
  geom_point(data = mean_highpv3baseline, aes(x = summ_visit_month, y = median_logpv3_serum), shape = 17, size = 2, color = "red", inherit.aes = FALSE) +
  scale_y_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text = element_text(size = 18))
ggsave("Figures/Serum_PV3_Neuts_highbaseline_all_medians.png", width = 6.3, height = 5.6, dpi = 300)
s3nhb
# 
# s3nlb <- ggplot(lowpv3baseline, aes(x =follow_up_days, y = logPV3_serum, group = pid)) +
#   geom_line(alpha = 0.5, size = 0.2) +
#   geom_point(alpha = 0.1) +
#   geom_point(data = mean_lowpv3baseline, aes(x = summ_visit_day, y = mean_logpv3_serum), shape = 17, size = 2, color = "red", inherit.aes = FALSE) +
#   geom_line(data = mean_lowpv3baseline, aes(x = summ_visit_day, y = mean_logpv3_serum), color = "red", inherit.aes = FALSE) +
#   scale_y_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12))
# s3nlb #saved as Serum_PV3_Neuts_lowbaseline_all

s3nhb<- ggplot(lowpv3baseline, aes(x =follow_up_months, y = logPV3_serum, group = pid)) +
  geom_line(alpha = 0.5, size = 0.2) +
  geom_point(alpha = 0.1) +
  geom_ribbon(data = mean_lowpv3baseline, aes(x = summ_visit_month, ymin = lq_logpv3_serum, ymax = uq_logpv3_serum), fill = "red", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = mean_lowpv3baseline, aes(x = summ_visit_month, y = median_logpv3_serum), color = "red",  inherit.aes = FALSE) +
  geom_point(data = mean_lowpv3baseline, aes(x = summ_visit_month, y = median_logpv3_serum), shape = 17, size = 2, color = "red", inherit.aes = FALSE) +
  scale_y_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text = element_text(size = 18))
ggsave("Figures/Serum_PV3_Neuts_lowbaseline_all_medians.png", width = 6.3, height = 5.6, dpi = 300)
s1nhb


#complete data participants neuts stratified by mucosal neut baseline
# s1nchb<- ggplot(completedata_highpv1baseline, aes(x =follow_up_days, y = logPV1_serum, group = pid)) +
#   geom_line(alpha = 0.5, size = 0.2) +
#   geom_point(alpha = 0.1) +
#   geom_hline(yintercept = 3, linetype = "dashed")+
#   geom_point(data = mean_completedata_highpv1baseline, aes(x = summ_visit_day, y = mean_logpv1_serum), shape = 17, size = 2, color = "red", inherit.aes = FALSE) +
#   geom_line(data = mean_completedata_highpv1baseline, aes(x = summ_visit_day, y = mean_logpv1_serum), color = "red", inherit.aes = FALSE) +
#   scale_y_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12))
# ggsave("Figures/Serum_PV1_Neuts_highbaseline_completedata.png", width = 6.3, height = 5.6, dpi = 300)
# s1nchb

s1nchb<- ggplot(completedata_highpv1baseline, aes(x =follow_up_months, y = logPV1_serum, group = pid)) +
  geom_line(alpha = 0.5, size = 0.2) +
  geom_point(alpha = 0.1) +
  geom_ribbon(data = mean_completedata_highpv1baseline, aes(x = summ_visit_month, ymin = lq_logpv1_serum, ymax = uq_logpv1_serum), fill = "red", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = mean_completedata_highpv1baseline, aes(x = summ_visit_month, y = median_logpv1_serum), color = "red", inherit.aes = FALSE) +
  geom_point(data = mean_completedata_highpv1baseline, aes(x = summ_visit_month, y = median_logpv1_serum), shape = 17, size = 2, color = "red", inherit.aes = FALSE) +
  scale_y_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text = element_text(size = 18))
ggsave("Figures/Serum_PV1_Neuts_completedata_highbaseline_all_medians.png", width = 6.3, height = 5.6, dpi = 300)
s1nchb 

# s1nclb <- ggplot(completedata_lowpv1baseline, aes(x =follow_up_days, y = logPV1_serum, group = pid)) +
#   geom_line(alpha = 0.5, size = 0.2) +
#   geom_point(alpha = 0.1) +
#   geom_hline(yintercept = 3, linetype = "dashed")+
#   geom_point(data = mean_completedata_lowpv1baseline, aes(x = summ_visit_day, y = mean_logpv1_serum), shape = 17, size = 2, color = "red", inherit.aes = FALSE) +
#   geom_line(data = mean_completedata_lowpv1baseline, aes(x = summ_visit_day, y = mean_logpv1_serum), color = "red", inherit.aes = FALSE) +
#   scale_y_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12))
# ggsave("Figures/Serum_PV1_Neuts_lowbaseline_completedata.png", width = 6.3, height = 5.6, dpi = 300)
# s1nclb


s1nclb<- ggplot(completedata_lowpv1baseline, aes(x =follow_up_months, y = logPV1_serum, group = pid)) +
  geom_line(alpha = 0.5, size = 0.2) +
  geom_point(alpha = 0.1) +
  geom_ribbon(data = mean_completedata_lowpv1baseline, aes(x = summ_visit_month, ymin = lq_logpv1_serum, ymax = uq_logpv1_serum), fill = "red", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = mean_completedata_lowpv1baseline, aes(x = summ_visit_month, y = median_logpv1_serum), color = "red", inherit.aes = FALSE) +
  geom_point(data = mean_completedata_lowpv1baseline, aes(x = summ_visit_month, y = median_logpv1_serum), shape = 17, size = 2, color = "red", inherit.aes = FALSE) +
  scale_y_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text = element_text(size = 18))
ggsave("Figures/Serum_PV1_Neuts_completedata_lowbaseline_all_medians.png", width = 6.3, height = 5.6, dpi = 300)
s1nclb 

# s3nchb<- ggplot(completedata_highpv3baseline, aes(x =follow_up_days, y = logPV3_serum, group = pid)) +
#   geom_line(alpha = 0.5, size = 0.2) +
#   geom_point(alpha = 0.1) +
#   geom_hline(yintercept = 3, linetype = "dashed")+
#   geom_point(data = mean_completedata_highpv3baseline, aes(x = summ_visit_day, y = mean_logpv3_serum), shape = 17, size = 2, color = "red", inherit.aes = FALSE) +
#   geom_line(data = mean_completedata_highpv3baseline, aes(x = summ_visit_day, y = mean_logpv3_serum), color = "red", inherit.aes = FALSE) +
#   scale_y_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12))
# ggsave("Figures/Serum_PV3_Neuts_highbaseline_completedata.png", width = 6.3, height = 5.6, dpi = 300)
# s3nchb 

s3nchb<- ggplot(completedata_highpv3baseline, aes(x =follow_up_months, y = logPV3_serum, group = pid)) +
  geom_line(alpha = 0.5, size = 0.2) +
  geom_point(alpha = 0.1) +
  geom_ribbon(data = mean_completedata_highpv3baseline, aes(x = summ_visit_month, ymin = lq_logpv3_serum, ymax = uq_logpv3_serum), fill = "red", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = mean_completedata_highpv3baseline, aes(x = summ_visit_month, y = median_logpv3_serum), color = "red", inherit.aes = FALSE) +
  geom_point(data = mean_completedata_highpv3baseline, aes(x = summ_visit_month, y = median_logpv3_serum), shape = 17, size = 2, color = "red", inherit.aes = FALSE) +
  scale_y_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text = element_text(size = 18))
ggsave("Figures/Serum_PV3_Neuts_completedata_highbaseline_all_medians.png", width = 6.3, height = 5.6, dpi = 300)
s3nchb 

# s3nclb <- ggplot(completedata_lowpv3baseline, aes(x =follow_up_days, y = logPV3_serum, group = pid)) +
#   geom_line(alpha = 0.5, size = 0.2) +
#   geom_point(alpha = 0.1) +
#   geom_hline(yintercept = 3, linetype = "dashed")+
#   geom_point(data = mean_completedata_lowpv3baseline, aes(x = summ_visit_day, y = mean_logpv3_serum), shape = 17, size = 2, color = "red", inherit.aes = FALSE) +
#   geom_line(data = mean_completedata_lowpv3baseline, aes(x = summ_visit_day, y = mean_logpv3_serum), color = "red", inherit.aes = FALSE) +
#   scale_y_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12))
# ggsave("Figures/Serum_PV3_Neuts_lowbaseline_completedata.png", width = 6.3, height = 5.6, dpi = 300)
# s3nclb

s3nclb<- ggplot(completedata_lowpv3baseline, aes(x =follow_up_months, y = logPV3_serum, group = pid)) +
  geom_line(alpha = 0.5, size = 0.2) +
  geom_point(alpha = 0.1) +
  geom_ribbon(data = mean_completedata_lowpv3baseline, aes(x = summ_visit_month, ymin = lq_logpv3_serum, ymax = uq_logpv3_serum), fill = "red", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = mean_completedata_lowpv3baseline, aes(x = summ_visit_month, y = median_logpv3_serum), color = "red", inherit.aes = FALSE) +
  geom_point(data = mean_completedata_lowpv3baseline, aes(x = summ_visit_month, y = median_logpv3_serum), shape = 17, size = 2, color = "red", inherit.aes = FALSE) +
  scale_y_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text = element_text(size = 18))
ggsave("Figures/Serum_PV3_Neuts_completedata_lowbaseline_all_medians.png", width = 6.3, height = 5.6, dpi = 300)
s3nclb 




