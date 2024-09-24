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

################################################################################
# ##Load the data from both CSV files and arrange it to plot
FA_spir_cf_vdps <- read.csv("./IRC740H_visit1_spiral/vdp_analysis_results_September2024/FA_combined_vdps.csv")
N4_spir_cf_vdps <- read.csv("./IRC740H_visit1_spiral/vdp_analysis_results_September2024/N4_combined_vdps.csv")

# Reshape Data to Long Format for ggplot2
FA_data_long <- FA_spir_cf_vdps %>%
  pivot_longer(cols = -Subject_id, names_to = "Analysis_Method", values_to = "VDP")
N4_data_long <- N4_spir_cf_vdps %>%
  pivot_longer(cols = -Subject_id, names_to = "Analysis_Method", values_to = "VDP")

# Calculate Mean and SD for each Analysis_Method
summary_FA <- FA_data_long %>%
  group_by(Analysis_Method) %>%
  summarise(Mean_VDP = mean(VDP), SD_VDP = sd(VDP) )
summary_FA
summary_N4 <- N4_data_long %>%
  group_by(Analysis_Method) %>%
  summarise(Mean_VDP = mean(VDP), SD_VDP = sd(VDP) )
summary_N4
# Plot Bar Chart with Error Bars (Mean ± SD) -- FA corrected
y_label_fa <- expression(bold("Mean VDP"[FA]* " (%)"))
bar_fa <- ggplot(summary_FA, aes(x = Analysis_Method, y = Mean_VDP, fill = Analysis_Method)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_errorbar(aes(ymin = Mean_VDP - SD_VDP, ymax = Mean_VDP + SD_VDP), 
                width = 0.2, position = position_dodge(0.7)) +
  ylim(-3, 25) +
  theme_bw() +
  ylab(y_label_fa) +
  xlab("Analysis Method") +
  # labs(title = expression(bold("Mean VDP"[FA]*" by Analysis Method (± SD)"))) +
  theme(text = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  scale_fill_brewer(palette = "Set2")
bar_fa

# Plot Bar Chart with Error Bars (Mean ± SD) -- N4 corrected
y_label_n4 <- expression(bold("Mean VDP"[N4]* " (%)"))
bar_n4 <- ggplot(summary_N4, aes(x = Analysis_Method, y = Mean_VDP, fill = Analysis_Method)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_errorbar(aes(ymin = Mean_VDP - SD_VDP, ymax = Mean_VDP + SD_VDP), 
                width = 0.2, position = position_dodge(0.7)) +
  ylim(-3, 25) +
  theme_bw() +
  ylab(y_label_n4) +
  xlab("Analysis Method") +
  # labs(title = expression(bold("Mean VDP"[N4]*" by Analysis Method (± SD)"))) +
  theme(text = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  scale_fill_brewer(palette = "Set2")
bar_n4

# Save the plot as a png file in the specified directory
ggsave("./zR_plots_4abs/spir_vdps_BarComp_FA.png", plot = bar_fa, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4abs/spir_vdps_BarComp_N4.png", plot = bar_n4, width = 4.5, height = 3.7, dpi = 300)

