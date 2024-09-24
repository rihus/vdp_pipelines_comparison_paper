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
FA_spir_cf_vdps <- read.csv("./IRC740H_visit1_spiral/vdp_analysis_results_September2024/FA_combined_vdps.xlsx")
N4_spir_cf_vdps <- read.csv("./IRC740H_2Dspiral_CF/snr_vdp_hvp_analysis_results_July2024/N4_corr_median_analysis_results.csv")

# Reshape Data to Long Format for ggplot2
data_long <- data %>%
  pivot_longer(cols = everything(), names_to = "Method", values_to = "VDP")

# Calculate Mean and SD for each Method
summary_data <- data_long %>%
  group_by(Method) %>%
  summarise(
    Mean_VDP = mean(VDP),
    SD_VDP = sd(VDP)
  )

# Plot Bar Chart with Error Bars (Mean ± SD)
ggplot(summary_data, aes(x = Method, y = Mean_VDP, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_errorbar(aes(ymin = Mean_VDP - SD_VDP, ymax = Mean_VDP + SD_VDP), 
                width = 0.2, position = position_dodge(0.7)) +
  theme_minimal() +
  labs(title = "Mean VDP by Method (± SD)", 
       x = "Method", 
       y = "Mean VDP (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set2")


# Save the plot as a png file in the specified directory
ggsave("./zR_plots_4abs/cartVSspir_vdp_linreg_plain.png", plot = linreg_vdp, width = 4.5, height = 3.7, dpi = 300)
