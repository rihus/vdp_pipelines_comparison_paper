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

##Colorblind friendly palette:
cbPalette <- c("#999999","#E69F00","#CC79A7","#F0E442","#56B4E9","#009E73","#0072B2","#D55E00")
##Define order of the data
plt_order <- c("Thresholding","Hierarchical","Adaptive","Percentile","Median","Reader") #
################################################################################
# ##Load the data from both CSV files and arrange it to plot
N4_gre_h_vdps <- read.csv("./age_matched_spiral_healthy/vdp_analysis_results_October2024/N4_combined_vdps.csv")
FA_gre_h_vdps <- read.csv("./age_matched_spiral_healthy/vdp_analysis_results_October2024/FA_combined_vdps.csv")
# ##Reshape Data to Long Format for ggplot2
N4_long <- N4_gre_h_vdps %>%
  pivot_longer(cols = -Subject_id, names_to = "AnalysisMethod", values_to = "VDP")
FA_long <- FA_gre_h_vdps %>%
  pivot_longer(cols = -Subject_id, names_to = "AnalysisMethod", values_to = "VDP")
##Set the order the order of boxes in plots
N4_long$AnalysisMethod <- factor(N4_long$AnalysisMethod, levels = plt_order)
FA_long$AnalysisMethod <- factor(FA_long$AnalysisMethod, levels = plt_order)
##############################################
ylabel_n4 <- expression(bold(VDP[N4]* " " * "(%)")) # Spiral
# #Box plots with p-values VDP - N4
vdpn4_bxp_nop <- ggpaired(N4_long, x = "AnalysisMethod", y = "VDP", fill = "AnalysisMethod",
                        palette = c("#999999","#E69F00","#CC79A7","#F0E442","#56B4E9","#009E73"), 
                        id = "Subject_id", width = 0.5,
                        ylim = c(0, 10), line.color = "black", line.size = 0.5,
                        legend = "none", xlab = "") + #
  ylab(ylabel_n4) +
  geom_point(aes(color = "black"), shape = 19, size = 3) +
  # geom_line(aes(group = paste(Correction, Category)), color = "#000000", size = 0.5) +
  scale_color_manual(values = c("#000000", "#009e73")) +
  theme(panel.border = element_rect(color = "#000000", fill = NA, linewidth = 1),
        axis.text = element_text(size = 14, color = "#000000", face = "bold"),
        axis.title = element_text(size = 14, color = "#000000", face = "bold"),
        axis.line.x = element_line(linewidth = 1),
        axis.line.y = element_line(linewidth = 1),
        axis.text.x = element_text(angle = 45, hjust = 1))
print(vdpn4_bxp_nop)

##############################################
ylabel_fa <- expression(bold(VDP[FA]* " " * "(%)"))
# #Box plots with p-values VDP - FA
vdpfa_bxp_nop <- ggpaired(FA_long, x = "AnalysisMethod", y = "VDP", fill = "AnalysisMethod",
                        palette = c("#999999","#E69F00","#CC79A7","#F0E442","#56B4E9","#009E73"),
                        id = "Subject_id", width = 0.5,
                        ylim = c(0, 50), line.color = "black", line.size = 0.5,
                        legend = "none", xlab = "") +
  ylab(ylabel_fa) +
  geom_point(aes(color = "black"), shape = 19, size = 3) +
  # geom_line(aes(group = paste(Correction, Category)), color = "#000000", size = 0.5) +
  scale_color_manual(values = c("#000000", "#009e73")) +
  theme(panel.border = element_rect(color = "#000000", fill = NA, linewidth = 1),
        axis.text = element_text(size = 14, color = "#000000", face = "bold"),
        axis.title = element_text(size = 14, color = "#000000", face = "bold"),
        axis.line.x = element_line(linewidth = 1),
        axis.line.y = element_line(linewidth = 1),
        axis.text.x = element_text(angle = 45, hjust = 1))
print(vdpfa_bxp_nop)

#Save the plot as a png file in the specified directory
ggsave("./zR_plots_4abs/spir_N4_healthy_vdps_cbxp.png", plot = vdpn4_bxp_nop, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4abs/spir_FA_healthy_vdps_cbxp.png", plot = vdpfa_bxp_nop, width = 4.5, height = 3.7, dpi = 300)

