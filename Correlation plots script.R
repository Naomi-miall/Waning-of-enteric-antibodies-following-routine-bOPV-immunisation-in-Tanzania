###R script creating box plots for TZA duration of immunity paper
#Naomi and Alex
#October 2025

#Set up: Libaries; working directory; theme:
### scatter plot ####
setwd ("C:/Users/NaomiMiall/Documents/Polio work/Tanzania")

## library
library(haven)
library(ggplot2)
library(ggtext)
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

correlations <- read_dta("correlations_dataset.dta")
#vis <- subset (longitudinal, visit=="bl"|visit=="w2"|visit=="w4"|visit=="w6"|visit=="w8")
#vis_ser <- subset (longitudinal, vis=="bl"|vis=="w4"|vis=="w8")

correlations$visit <- factor(correlations$visit, levels = c("1", "2", "3", "4", "5", "6", "7", "8"),
                             labels = c("Baseline", "One month", "3 months", "6 months", "1 year", "2 years", "3 years", "8-9 years"))

##################################################################################################
#Stool and serum neuts
  df_cor_pv1 <-
 correlations|>
  summarise(N=n(),
            p=cor.test(logPV1_stool,logPV1_serum,method='kendall',exact=F)$p.value,
            tau_b=cor.test(logPV1_stool,logPV1_serum,method='kendall',exact=F)$est,
            p_num=p,
            p=ifelse(p<0.001,"p<0.001",paste0("p=",round(p,3))))

correlations |>
  ggplot()+geom_jitter(aes(x=logPV1_stool,y=logPV1_serum,colour=factor(sampletype)),shape=1, size=3)+
  theme(axis.title.x = element_markdown(size = 14), axis.title.y = element_markdown(size = 14), axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14),plot.title = element_markdown(hjust = 0.5, size = 16, face = "bold"),plot.title.position = "plot")+
  scale_color_manual(name = "Sample", values = c("1" = "#1f78b4", "2" = "#E66218"),
                     labels = c("1" = "Cross sectional", "2" = "Longitudinal")) +
  geom_text(data=df_cor_pv1 |>
              mutate(est=paste0("Tau_b=",round(tau_b,2)),
                     lab=paste0(est, ", ", p)),
            aes(x=8,y=5.5,label=lab), size=5, fontface="bold")+
  labs(y="Serum log<sub>2</sub> neutralisation titre",x="Stool log<sub>2</sub> neutralisation titre",title="Correlation between stool and serum <br> PV1-specific log<sub>2</sub> neutralisation titre")
ggsave("Figures/Correlations_pv1_sample.png", width = 6, height = 4, dpi = 300)

  df_cor_pv3 <-
    correlations|>
    summarise(N=n(),
              p=cor.test(logPV3_stool,logPV3_serum,method='kendall',exact=F)$p.value,
              tau_b=cor.test(logPV3_stool,logPV3_serum,method='kendall',exact=F)$est,
              p_num=p,
              p=ifelse(p<0.001,"p<0.001",paste0("p=",round(p,3))))
  

  correlations |>
    ggplot()+geom_jitter(aes(x=logPV3_stool,y=logPV3_serum,colour=factor(sampletype)),shape=1, size=3)+
    theme(axis.title.x = element_markdown(size = 14), axis.title.y = element_markdown(size = 14), axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14),plot.title = element_markdown(hjust = 0.5, size = 16, face = "bold"),plot.title.position = "plot")+
    scale_color_manual(name = "Sample", values = c("1" = "#1f78b4", "2" = "#E66218"),
                       labels = c("1" = "Cross sectional", "2" = "Longitudinal")) +
    geom_text(data=df_cor_pv3 |>
                mutate(est=paste0("Tau_b=",round(tau_b,2)),
                       lab=paste0(est, ", ", p)),
              aes(x=8,y=5.5,label=lab), size=5, fontface="bold")+
    labs(y="Serum log<sub>2</sub> neutralisation titre",x="Stool log<sub>2</sub> neutralisation titre",title="Correlation between stool and serum <br> PV3-specific log<sub>2</sub> neutralisation titre")
  ggsave("Figures/Correlations_pv3_sample.png", width = 6, height = 4, dpi = 300)
  
  ###########################################################################
  
 #PV1 and 3 neuts
  df_cor_stool <- 
    correlations |>
    summarise(N=n(),
              p=cor.test(logPV1_stool,logPV3_stool,method='kendall',exact=F)$p.value,
              tau_b=cor.test(logPV1_stool,logPV3_stool,method='kendall',exact=F)$est,
              p_num=p,
              p=ifelse(p<0.001,"p<0.001",paste0("p=",round(p,3))))
  
  correlations |>
    ggplot()+geom_jitter(aes(x=logPV1_stool,y=logPV3_stool,colour=factor(sampletype)),shape=1, size=3)+
    theme(axis.title.x = element_markdown(size = 14), axis.title.y = element_markdown(size = 14), axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14),plot.title = element_markdown(hjust = 0.5, size = 16, face = "bold"),plot.title.position = "plot")+
    scale_color_manual(name = "Sample", values = c("1" = "#1f78b4", "2" = "#E66218"),
                       labels = c("1" = "Cross sectional", "2" = "Longitudinal")) +
    geom_text(data=df_cor_stool |>
                mutate(est=paste0("Tau_b=",round(tau_b,2)),
                       lab=paste0(est, ", ", p)),
              aes(x=8,y=1.5,label=lab), size=5, fontface="bold")+
    labs(y="PV3 log<sub>2</sub> neutralisation titre",x="PV1 log<sub>2</sub> neutralisation titre",title="Correlation between PV1 and PV3 <br> log<sub>2</sub> neutralisation titre in stool samples")
  ggsave("Figures/Correlations_stool_pv.png", width = 6, height = 4, dpi = 300)
  
  
  
  df_cor_serum <-
    correlations |>
    summarise(N=n(),
              p=cor.test(logPV1_serum,logPV3_serum,method='kendall',exact=F)$p.value,
              tau_b=cor.test(logPV1_serum,logPV3_serum,method='kendall',exact=F)$est,
              p_num=p,
              p=ifelse(p<0.001,"p<0.001",paste0("p=",round(p,3))))
  
  
  correlations |>
    ggplot()+geom_jitter(aes(x=logPV1_serum,y=logPV3_serum,colour=factor(sampletype)),shape=1, size=3)+
    theme(axis.title.x = element_markdown(size = 14), axis.title.y = element_markdown(size = 14), axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14),plot.title = element_markdown(hjust = 0.5, size = 16, face = "bold"),plot.title.position = "plot")+
    scale_color_manual(name = "Sample", values = c("1" = "#1f78b4", "2" = "#E66218"),
                       labels = c("1" = "Cross sectional", "2" = "Longitudinal")) +
    geom_text(data=df_cor_serum|>
                mutate(est=paste0("Tau_b=",round(tau_b,2)),
                       lab=paste0(est, ", ", p)),
              aes(x=11,y=5.3,label=lab), size=5, fontface="bold")+
    labs(y="PV3 log<sub>2</sub> neutralisation titre",x="PV1 log<sub>2</sub> neutralisation titre",title="Correlation between PV1 and PV3 <br> log<sub>2</sub> neutralisation titre in serum samples")
  ggsave("Figures/Correlations_serum_pv.png", width = 6, height = 4, dpi = 300)
  
  
#########################################################################
  #Neuts and IgA in stool
  df_cor_assay_stool_pv1 <- 
    correlations |>
    summarise(N=n(),
              p=cor.test(logPV1_stool,logPV1_IgA_stool,method='kendall',exact=F)$p.value,
              tau_b=cor.test(logPV1_stool,logPV1_IgA_stool,method='kendall',exact=F)$est,
              p_num=p,
              p=ifelse(p<0.001,"p<0.001",paste0("p=",round(p,3))))
  
  
  correlations |>
    ggplot()+geom_jitter(aes(x=logPV1_stool,y=logPV1_IgA_stool,colour=factor(sampletype)),shape=1, size=3)+
    theme(axis.title.x = element_markdown(size = 14), axis.title.y = element_markdown(size = 14), axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14),plot.title = element_markdown(hjust = 0.5, size = 16, face = "bold"),plot.title.position = "plot")+
    scale_color_manual(name = "Sample", values = c("1" = "#1f78b4", "2" = "#E66218"),
                       labels = c("1" = "Cross sectional", "2" = "Longitudinal")) +
    geom_text(data=df_cor_assay_stool_pv1 |>
                mutate(est=paste0("Tau_b=",round(tau_b,2)),
                       lab=paste0(est, ", ", p)),
              aes(x=8,y=1.25,label=lab), size=5, fontface="bold")+  
    scale_y_continuous(breaks = c(1.0, 1.5, 2.0, 2.5, 3.0, 3.5), limits = c(1, 3.5)) +
    labs(y="PV1 log<sub>10</sub> IgA MFI",x="PV1 log<sub>2</sub> neutralisation titre",title="Correlation between PV1 <br> log<sub>2</sub> neutralisation titre and <br> log<sub>10</sub> IgA MFI in stool samples")
  ggsave("Figures/Correlations_stool_assay_pv1.png", width = 6, height = 4, dpi = 300)
  
  df_cor_assay_stool_pv3 <- 
    correlations |>
    summarise(N=n(),
              p=cor.test(logPV3_stool,logPV3_IgA_stool,method='kendall',exact=F)$p.value,
              tau_b=cor.test(logPV3_stool,logPV3_IgA_stool,method='kendall',exact=F)$est,
              p_num=p,
              p=ifelse(p<0.001,"p<0.001",paste0("p=",round(p,3))))
  
  
  correlations |>
    ggplot()+geom_jitter(aes(x=logPV3_stool,y=logPV3_IgA_stool,colour=factor(sampletype)),shape=1, size=3)+
    theme(axis.title.x = element_markdown(size = 14), axis.title.y = element_markdown(size = 14), axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14),plot.title = element_markdown(hjust = 0.5, size = 16, face = "bold"),plot.title.position = "plot")+
    scale_color_manual(name = "Sample", values = c("1" = "#1f78b4", "2" = "#E66218"),
                       labels = c("1" = "Cross sectional", "2" = "Longitudinal")) +
    geom_text(data=df_cor_assay_stool_pv3 |>
                mutate(est=paste0("Tau_b=",round(tau_b,2)),
                       lab=paste0(est, ", ", p)),
              aes(x=8,y=1.25,label=lab), size=5, fontface="bold")+
    scale_y_continuous(breaks = c(1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0), limits = c(1, 3.5)) +
    labs(y="PV3 log<sub>10</sub> IgA MFI",x="PV3 log<sub>2</sub> neutralisation titre",title="Correlation between PV3 <br> log<sub>2</sub> neutralisation titre and <br> log<sub>10</sub> IgA MFI in stool samples")
  ggsave("Figures/Correlations_stool_assay_pv3.png", width = 6, height = 4, dpi = 300)
  
#Neuts and IgA in serum
  df_cor_assay_serum_pv1 <- 
    correlations |>
    summarise(N=n(),
              p=cor.test(logPV1_serum,logPV1_IgA_serum,method='kendall',exact=F)$p.value,
              tau_b=cor.test(logPV1_serum,logPV1_IgA_serum,method='kendall',exact=F)$estimate,
              p_num=p,
              p=ifelse(p<0.001,"p<0.001",paste0("p=",round(p,3))))
  
  correlations |>
    ggplot()+geom_jitter(aes(x=logPV1_serum,y=logPV1_IgA_serum,colour=factor(sampletype)),shape=1, size=3)+
    theme(axis.title.x = element_markdown(size = 14), axis.title.y = element_markdown(size = 14), axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14),plot.title = element_markdown(hjust = 0.5, size = 16, face = "bold"),plot.title.position = "plot")+
    scale_color_manual(name = "Sample", values = c("1" = "#1f78b4", "2" = "#E66218"),
                       labels = c("1" = "Cross sectional", "2" = "Longitudinal")) +
    geom_text(data=df_cor_assay_serum_pv1 |>
                mutate(est=paste0("Tau_b=",round(tau_b,2)),
                       lab=paste0(est, ", ", p)),
              aes(x=8,y=1.75,label=lab), size=5, fontface="bold")+
    scale_y_continuous(breaks = c(1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0), limits = c(1.5, 4)) +
    labs(y="PV1 log<sub>10</sub> IgA MFI",x="PV1 log<sub>2</sub> neutralisation titre",title="Correlation between PV1 <br> log<sub>2</sub> neutralisation titre and <br> log<sub>10</sub> IgA MFI in serum samples")
  ggsave("Figures/Correlations_serum_assay_pv1.png", width = 6, height = 4, dpi = 300)
  
  df_cor_assay_serum_pv3 <- 
    correlations |>
    summarise(N=n(),
              p=cor.test(logPV3_serum,logPV3_IgA_serum,method='kendall',exact=F)$p.value,
              tau_b=cor.test(logPV3_serum,logPV3_IgA_serum,method='kendall',exact=F)$est,
              p_num=p,
              p=ifelse(p<0.001,"p<0.001",paste0("p=",round(p,3))))
  
  
  correlations |>
    ggplot()+geom_jitter(aes(x=logPV3_serum,y=logPV3_IgA_serum,colour=factor(sampletype)),shape=1, size=3)+
    theme(axis.title.x = element_markdown(size = 14), axis.title.y = element_markdown(size = 14), axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14),plot.title = element_markdown(hjust = 0.5, size = 16, face = "bold"),plot.title.position = "plot")+
    scale_color_manual(name = "Sample", values = c("1" = "#1f78b4", "2" = "#E66218"),
                       labels = c("1" = "Cross sectional", "2" = "Longitudinal")) +
    geom_text(data=df_cor_assay_serum_pv3 |>
                mutate(est=paste0("Tau_b=",round(tau_b,2)),
                       lab=paste0(est, ", ", p)),
              aes(x=8,y=1.75,label=lab), size=5, fontface="bold")+   
    scale_y_continuous(breaks = c(1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0), limits = c(1.5, 4)) +
    labs(y="PV3 log<sub>10</sub> IgA MFI",x="PV3 log<sub>2</sub> neutralisation titre",title="Correlation between PV3 <br> log<sub>2</sub> neutralisation titre and <br> log<sub>10</sub> IgA MFI in serum samples")
  ggsave("Figures/Correlations_serum_assay_pv3.png", width = 6, height = 4, dpi = 300)
  
  
  #Neuts and IgG in serum
  df_cor_IgG_pv1 <- 
    correlations |>
    summarise(N=n(),
              p=cor.test(logPV1_serum,logPV1_IgG_serum,method='kendall',exact=F)$p.value,
              tau_b=cor.test(logPV1_serum,logPV1_IgG_serum,method='kendall',exact=F)$estimate,
              p_num=p,
              p=ifelse(p<0.001,"p<0.001",paste0("p=",round(p,3))))
  
  correlations |>
    ggplot()+geom_jitter(aes(x=logPV1_serum,y=logPV1_IgG_serum,colour=factor(sampletype)),shape=1, size=3)+
    theme(axis.title.x = element_markdown(size = 14), axis.title.y = element_markdown(size = 14), axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14),plot.title = element_markdown(hjust = 0.5, size = 16, face = "bold"),plot.title.position = "plot")+
    scale_color_manual(name = "Sample", values = c("1" = "#1f78b4", "2" = "#E66218"),
                       labels = c("1" = "Cross sectional", "2" = "Longitudinal")) +
    geom_text(data=df_cor_IgG_pv1 |>
                mutate(est=paste0("Tau_b=",round(tau_b,2)),
                       lab=paste0(est, ", ", p)),
              aes(x=8,y=1.75,label=lab), size=5, fontface="bold")+
    scale_y_continuous(breaks = c(1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0), limits = c(1.5, 4)) +
    labs(y="PV1 log<sub>10</sub> IgG MFI",x="PV1 log<sub>2</sub> neutralisation titre",title="Correlation between PV1 <br> log<sub>2</sub> neutralisation titre and <br> log<sub>10</sub> IgG MFI in serum samples")
  ggsave("Figures/Correlations_IgG_pv1.png", width = 6, height = 4, dpi = 300)
  
  df_cor_IgG_pv3 <- 
    correlations |>
    summarise(N=n(),
              p=cor.test(logPV3_serum,logPV3_IgG_serum,method='kendall',exact=F)$p.value,
              tau_b=cor.test(logPV3_serum,logPV3_IgG_serum,method='kendall',exact=F)$estimate,
              p_num=p,
              p=ifelse(p<0.001,"p<0.001",paste0("p=",round(p,3))))
  
  correlations |>
    ggplot()+geom_jitter(aes(x=logPV3_serum,y=logPV3_IgG_serum,colour=factor(sampletype)),shape=1, size=3)+
    theme(axis.title.x = element_markdown(size = 14), axis.title.y = element_markdown(size = 14), axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14),plot.title = element_markdown(hjust = 0.5, size = 16, face = "bold"),plot.title.position = "plot")+
    scale_color_manual(name = "Sample", values = c("1" = "#1f78b4", "2" = "#E66218"),
                       labels = c("1" = "Cross sectional", "2" = "Longitudinal")) +
    geom_text(data=df_cor_IgG_pv3 |>
                mutate(est=paste0("Tau_b=",round(tau_b,2)),
                       lab=paste0(est, ", ", p)),
              aes(x=8,y=1.75,label=lab), size=5, fontface="bold")+
    scale_y_continuous(breaks = c(1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0), limits = c(1.5, 4)) +
    labs(y="PV3 log<sub>10</sub> IgG MFI",x="PV3 log<sub>2</sub> neutralisation titre",title="Correlation between PV3 <br> log<sub>2</sub> neutralisation titre and <br> log<sub>10</sub> IgG MFI in serum samples")
  ggsave("Figures/Correlations_IgG_pv3.png", width = 6, height = 4, dpi = 300)
  
#IgA and IgG in serum
  df_cor_IgGA_pv1 <- 
    correlations |>
    summarise(N=n(),
              p=cor.test(logPV1_IgA_serum,logPV1_IgG_serum,method='kendall',exact=F)$p.value,
              tau_b=cor.test(logPV1_IgA_serum,logPV1_IgG_serum,method='kendall',exact=F)$estimate,
              p_num=p,
              p=ifelse(p<0.001,"p<0.001",paste0("p=",round(p,3))))
  
  correlations |>
    ggplot()+geom_jitter(aes(x=logPV1_IgA_serum,y=logPV1_IgG_serum,colour=factor(sampletype)),shape=1, size=3)+
    theme(axis.title.x = element_markdown(size = 14), axis.title.y = element_markdown(size = 14), axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14),plot.title = element_markdown(hjust = 0.5, size = 16, face = "bold"),plot.title.position = "plot")+
    scale_color_manual(name = "Sample", values = c("1" = "#1f78b4", "2" = "#E66218"),
                       labels = c("1" = "Cross sectional", "2" = "Longitudinal")) +
    geom_text(data=df_cor_IgGA_pv1 |>
                mutate(est=paste0("Tau_b=",round(tau_b,2)),
                       lab=paste0(est, ", ", p)),
              aes(x=3,y=1.75,label=lab), size=5, fontface="bold")+
    scale_y_continuous(breaks = c(1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0), limits = c(1.5, 4)) +
    scale_x_continuous(breaks = c(1.5, 2.0, 2.5, 3.0, 3.5), limits = c(1.5,3.75)) +
    labs(y="PV1 log<sub>10</sub> IgG MFI",x="PV1 log<sub>10</sub> IgA MFI",title="Correlation between PV1 <br> log<sub>10</sub> IgA MFI and <br> log<sub>10</sub> IgG MFI in serum samples")
  ggsave("Figures/Correlations_IgG_IgA_pv1.png", width = 6, height = 4, dpi = 300)
  
  df_cor_IgGA_pv3 <- 
    correlations |>
    summarise(N=n(),
              p=cor.test(logPV3_IgA_serum,logPV3_IgG_serum,method='kendall',exact=F)$p.value,
              tau_b=cor.test(logPV3_IgA_serum,logPV3_IgG_serum,method='kendall',exact=F)$estimate,
              p_num=p,
              p=ifelse(p<0.001,"p<0.001",paste0("p=",round(p,3))))
  
  correlations |>
    ggplot()+geom_jitter(aes(x=logPV3_IgA_serum,y=logPV3_IgG_serum,colour=factor(sampletype)),shape=1, size=3)+
    theme(axis.title.x = element_markdown(size = 14), axis.title.y = element_markdown(size = 14), axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14),plot.title = element_markdown(hjust = 0.5, size = 16, face = "bold"),plot.title.position = "plot")+
    scale_color_manual(name = "Sample", values = c("1" = "#1f78b4", "2" = "#E66218"),
                       labels = c("1" = "Cross sectional", "2" = "Longitudinal")) +
    geom_text(data=df_cor_IgGA_pv3 |>
                mutate(est=paste0("Tau_b=",round(tau_b,2)),
                       lab=paste0(est, ", ", p)),
              aes(x=3,y=1.75,label=lab), size=5, fontface="bold")+
    scale_y_continuous(breaks = c(1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0), limits = c(1.5, 4)) +
    scale_x_continuous(breaks = c(1.5, 2.0, 2.5, 3.0, 3.5), limits = c(1.5,3.75)) +
    labs(y="PV3 log<sub>10</sub> IgG MFI",x="PV3 log<sub>10</sub> IgA MFI",title="Correlation between PV3 <br> log<sub>10</sub> IgA MFI and <br> log<sub>10</sub> IgG MFI in serum samples")
  ggsave("Figures/Correlations_IgG_IgA_pv3.png", width = 6, height = 4, dpi = 300)
  
  
