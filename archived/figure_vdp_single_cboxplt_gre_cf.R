#load ggplot2
library(ggplot2)
library(ggpubr)
library(blandr)
library(dplyr)
library(tidyr)
library(rstatix)

# Set the working directory to the path where your CSV files are located
setwd("C:/Users/HUSDQ4/OneDrive - cchmc/cincy_work/all_projects_data_work/vdp_analysis/analysis_comparisons")
# linetypes: 0 = blank, 1 = solid, 2 = dashed, 3 = dotted, 4 = dotdash, 5 = longdash, 6 = twodash

#Colorblind friendly palette:
cbPalette <- c("#999999","#E69F00","#CC79A7","#F0E442","#56B4E9","#009E73","#0072B2","#D55E00")
# Define order of the data
plt_order <- c("Thresholding","Hierarchical","Adaptive","Percentile","Median","Reader")
################################################################################
# ##Load the data from both CSV files and arrange it to plot
N4_gre_cf_vdps <- read.csv("./HyPOINT_phase1_cart_cf/vdp_analysis_results_September2024/N4_combined_vdps.csv")

# Reshape Data to Long Format for ggplot2
N4_data_long <- N4_gre_cf_vdps %>%
  pivot_longer(cols = -Subject_id, names_to = "Analysis_Method", values_to = "VDP")

# Calculate Mean and SD for each Analysis_Method
summary_N4 <- N4_data_long %>%
  group_by(Analysis_Method) %>%
  summarise(Mean_VDP = mean(VDP), SD_VDP = sd(VDP) )
summary_N4
##Make the desired order a factor
summary_N4$Analysis_Method <- factor(summary_N4$Analysis_Method, levels = plt_order)

##############################################
y_label <- expression(bold(VDP[GRE-N4] * "(%)"))
# #Box plots with p-values VDP - N4
vdp_bxp_nop <- ggpaired(N4_data_long, x = "Analysis_Method", y = "VDP", fill = "Analysis_Method",
                        palette = c("#999999","#E69F00","#CC79A7","#F0E442","#56B4E9","#009E73"),
                        id = "Subject_id", width = 0.5,
                        ylim = c(0, 70), line.color = "black", line.size = 0.5,
                        legend = "none", xlab = "") +
  ylab(y_label) +
  geom_point(aes(color = "black"), shape = 19, size = 3) +
  # geom_line(aes(group = paste(Correction, Category)), color = "#000000", size = 0.5) +
  scale_color_manual(values = c("#000000", "#009e73")) +
  theme(panel.border = element_rect(color = "#000000", fill = NA, linewidth = 1),
        axis.text = element_text(size = 14, color = "#000000", face = "bold"),
        axis.title = element_text(size = 14, color = "#000000", face = "bold"),
        axis.line.x = element_line(linewidth = 1),
        axis.line.y = element_line(linewidth = 1),
        axis.text.x = element_text(angle = 45, hjust = 1))
vdp_bxp_nop

# Save the plot as a png file in the specified directory
ggsave("./zR_plots_4ppr/gre_N4_vdp_cbxp_nop.png", plot = vdp_bxp_nop, width = 4.5, height = 3.7, dpi = 300)

# #Box plots with p-values VDP - FA
vdp_bxp_nop <- ggpaired(df_3visits_FA, x = "VISIT", y = "VDP", fill = "VISIT",
                        palette = c("#004166", "#0072b2", "#00a3ff"), id = "Subject_id", width = 0.5,
                        ylim = c(0, 40), line.color = "black", line.size = 0.5,
                        legend = "none", xlab = "") +
  ylab(y_label) +
  geom_point(aes(color = Category), shape = 19, size = 3) +
  # geom_line(aes(group = paste(Correction, Category)), color = "#000000", size = 0.5) +
  scale_color_manual(values = c("#000000", "#009e73")) +
  theme(panel.border = element_rect(color = "#000000", fill = NA, linewidth = 1),
        axis.text = element_text(size = 22, color = "#000000", face = "bold"),
        axis.title = element_text(size = 22, color = "#000000", face = "bold"),
        axis.line.x = element_line(linewidth = 1), axis.line.y = element_line(linewidth = 1))
vdp_bxp_nop

# Save the plot as a png file in the specified directory
ggsave("./zR_plots_4ppr/vdp_FA_3visits_cbxp_nop.png", plot = vdp_bxp_nop, width = 4.5, height = 3.7, dpi = 300)



