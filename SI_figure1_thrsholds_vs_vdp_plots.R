#If needed, install library/package
#install.packages("ggpubr")

#Load library
library(ggpubr)
library(tidyverse)
library(rstatix)
library(ggsignif)
library(ggplot2)
library(dplyr)
library(rlang)

################################################################################
# #Set the working directory to the path where your CSV files are located
setwd("C:/Users/HUSDQ4/OneDrive - cchmc/cincy_work/all_projects_data_work/vdp_analysis/analysis_comparisons")

################################################################################
# ##Load the data from both CSV files and arrange it to plot
FA_vdp <- read.csv("./IRC740H_visit1_spiral_cf/vdp_analysis_results_December2024/FA_meanVDP_thresholds.csv")
N4_vdp <- read.csv("./IRC740H_visit1_spiral_cf/vdp_analysis_results_December2024/N4_meanVDP_thresholds.csv")

FA_plot <- ggplot(FA_vdp, aes(x = Threshold, y = VDP, color = Category, 
                              shape = Category, linetype = Category)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = VDP - SD, ymax = VDP + SD), width = 0.02) +
  geom_line() +
  geom_vline(xintercept = 40, linetype = "dashed", color = "black", linewidth = 0.5) +
  labs(title = " ", x = "Threshold", y = expression(bold(Mean *" VDP"[{"FA "}] * "(%)"))) +
  ylim(-10, 59) + xlim(3, 100) +
  theme_bw() +
  theme(panel.border = element_rect(color = "#000000", fill = NA, linewidth = 1),
        axis.text = element_text(size = 22, color = "#000000", face = "bold"),
        axis.title = element_text(size = 22, color = "#000000", face = "bold"),
        axis.line.x = element_line(linewidth = 1),
        axis.line.y = element_line(linewidth = 1),
        legend.title = element_blank(),
        legend.position = c(0.2, 0.8),
        legend.background = element_blank(),
        legend.text = element_text(size = 14, face = "bold")) +
  scale_color_manual(values = c("Healthy" = "darkgreen", "CF" = "#0072b2", "Difference" = "black"),
                    breaks = c("Healthy", "CF", "Difference")) +
  scale_shape_manual(values = c("Healthy" = 16, "CF" = 5, "Difference" = 17),
                    breaks = c("Healthy", "CF", "Difference")) +
  scale_linetype_manual(values = c("Healthy" = "solid", "CF" = "dotted", "Difference" = "dashed"),
                    breaks = c("Healthy", "CF", "Difference"))
print(FA_plot)

N4_plot <- ggplot(N4_vdp, aes(x = Threshold, y = VDP, color = Category, 
                              shape = Category, linetype = Category)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = VDP - SD, ymax = VDP + SD), width = 0.02) +
  geom_line() +
  geom_vline(xintercept = 60, linetype = "dashed", color = "black", linewidth = 0.5) +
  labs(title = " ", x = "Threshold", y = expression(bold(Mean *" VDP"[{"N4 "}] * "(%)"))) +
  ylim(-10, 59) + xlim(3, 100) +
  theme_bw() +
  theme(panel.border = element_rect(color = "#000000", fill = NA, linewidth = 1),
        axis.text = element_text(size = 22, color = "#000000", face = "bold"),
        axis.title = element_text(size = 22, color = "#000000", face = "bold"),
        axis.line.x = element_line(linewidth = 1),
        axis.line.y = element_line(linewidth = 1),
        legend.title = element_blank(),
        legend.position = c(0.2, 0.8),
        legend.background = element_blank(),
        legend.text = element_text(size = 14, face = "bold")) +
  scale_color_manual(values = c("Healthy" = "darkgreen", "CF" = "#cc79a7", "Difference" = "black"),
                     breaks = c("Healthy", "CF", "Difference")) +
  scale_shape_manual(values = c("Healthy" = 16, "CF" = 5, "Difference" = 17),
                     breaks = c("Healthy", "CF", "Difference")) +
  scale_linetype_manual(values = c("Healthy" = "solid", "CF" = "dotted", "Difference" = "dashed"),
                        breaks = c("Healthy", "CF", "Difference"))
print(N4_plot)


ggsave("./zR_plots_4ppr/si_FA_meanVDP_Thresholds_comparison.png", plot = FA_plot, width = 4.9, height = 3.2, dpi = 300)
ggsave("./zR_plots_4ppr/si_N4_meanVDP_Thresholds_comparison.png", plot = N4_plot, width = 4.9, height = 3.2, dpi = 300)

