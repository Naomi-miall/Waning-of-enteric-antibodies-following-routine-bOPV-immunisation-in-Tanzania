##Forest plot code

setwd ("C:/Users/NaomiMiall/Documents/Polio work/Tanzania")

library(forestplot)
library(dplyr)
library(readxl)
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


df <- read_dta("Forest_plot_2.dta")

df$var <- recode(df$var,
                 "bf_4months" = "Non-exclusive breastfeeding",
                 "stunting" = "Stunting",
                 "watersource_bin" = "At least 15 minutes to water-source",
                 "hiv_exposure" = "HIV exposed",
                 "coinfect_enteroinvas" = "Enteroinvasive bacteria",
                 "coinfect_otherbacter" = "Other bacteria",
                 "coinfect_virus" = "Virus",
                 "coinfect_caa" = "Schistosomiasis CAA")

df$var <- factor(df$var, levels = c("Non-exclusive breastfeeding", "At least 15 minutes to water-source","Stunting","HIV exposed", "Schistosomiasis CAA","Enteroinvasive bacteria","Other bacteria","Virus"),
                 labels = c("Non-exclusive breastfeeding", "At least 15 minutes to water-source","Stunting","HIV exposed", "Schistosomiasis CAA","Enteroinvasive bacteria","Other bacteria","Virus"))
df <- df[order(df$var), ]
df$var <- as.character(df$var)


label_matrix <- rbind(
  c("Variable", "PV1 OR [95% Ci], p", "PV3 OR [95% Ci], p"),  # updated header
  cbind(
    df$var,
    ifelse(is.na(df$OR1),
           "NA",
           sprintf("%.2f [%.2f–%.2f], p=%.2f", df$OR1, df$Lower_confidence_interval1, df$Upper_confidence_interval1, df$p1)),
    ifelse(is.na(df$OR3),
           "NA",
           sprintf("%.2f [%.2f–%.2f], p=%.2f", df$OR3, df$Lower_confidence_interval3, df$Upper_confidence_interval3, df$p3))
  )
)


or_matrix <- rbind(
  c(NA, NA),  # header row
  cbind(df$OR1, df$OR3)
)

lower_matrix <- rbind(
  c(NA, NA),
  cbind(df$Lower_confidence_interval1, df$Lower_confidence_interval3)
)

upper_matrix <- rbind(
  c(NA, NA),
  cbind(df$Upper_confidence_interval1, df$Upper_confidence_interval3)
)

forestplot <- forestplot(labeltext = label_matrix,
           mean = or_matrix,
           lower = lower_matrix,
           upper = upper_matrix,
           xlog = TRUE,
           vertices = TRUE,
           is.summary = c(TRUE, rep(FALSE, nrow(df))),
           zero = 1,
           col = fpColors(box = c("#9F2042", "#38618C"),
                          lines = c("#9F2042", "#38618C")),
           boxsize = 0.2,
           clip = c(0.25, 4),
           xticks = c(0.25, 0.5, 1, 2, 3, 4),
           legend = c("Positive neutralisation for poliovirus type-1", "Positive neutralisation for poliovirus type-3"),
           txt_gp = fpTxtGp(label = gpar(fontfamily = "", cex = 0.9)))
forestplot

