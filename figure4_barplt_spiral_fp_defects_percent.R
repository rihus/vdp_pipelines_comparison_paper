#load ggplot2
library(ggplot2)
library(ggpubr)
library(blandr)
library(dplyr)
library(tidyr)
library(rstatix)
library(PMCMRplus)

## Set the working directory to the path where your CSV files are located
setwd("C:/Users/HUSDQ4/OneDrive - cchmc/cincy_work/all_projects_data_work/vdp_analysis/analysis_comparisons")
## linetypes: 0 = blank, 1 = solid, 2 = dashed, 3 = dotted, 4 = dotdash, 5 = longdash, 6 = twodash

##Colorblind friendly palette:
#cbPalette <- c("#999999","#F0E442","#CC79A7","#E69F00","#D55E00","#56B4E9","#0072B2") # "#009E73",
cbPalette <- c("#CC6677", "#882255", "#88CCEE", "#332288", "#999933", "#888888",
               "#AA4499", "#DDCC77", "#661100", "#6699CC", "#44AA99", "#117733")
##Define order of the data
plt_order <- c("Reader", "Adaptive","Hierarchical","Thresholding","Percentile","Median")
create_barplot <- function(input_df, xdata, ydata, y_label, x_label, ylimits,
                           barclr, subjs_colname, plt_order){
  # Convert xdata to a factor in the input dataframe if not already
  input_df[[xdata]] <- factor(input_df[[xdata]], levels = plt_order)
  ## Calculate basic stat properties
  summary_stat <- input_df %>%
    group_by(!!sym(xdata)) %>%
    summarise(Mean_y = mean(!!sym(ydata)), SD_y = sd(!!sym(ydata)))
  print(summary_stat)
  # Ensure xdata is a factor with the given order
  summary_stat[[xdata]] <- factor(summary_stat[[xdata]], levels = plt_order)
  ## Plot Bar Chart with Error Bars (Mean ± SD)
  bar_plt <- ggplot(summary_stat, aes(x = !!sym(xdata), y = Mean_y, fill = !!sym(xdata))) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    geom_errorbar(aes(ymin = Mean_y - SD_y, ymax = Mean_y + SD_y), 
                  width = 0.2, position = position_dodge(0.7)) +
    geom_point(data = input_df, aes(x = !!sym(xdata), y = !!sym(ydata)), color = "#000000", size = 1) +
    ylim(ylimits) +
    theme_bw() +
    ylab(y_label) +
    xlab(x_label) +
    theme(text = element_text(size = 14, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    scale_fill_manual(values = barclr)
  print(bar_plt)
  return(bar_plt)
}
################################################################################
# ##Load the fraction of reader defect data from CSV files
FA_spir_cf_pfp <- read.csv("./IRC740H_visit1_spiral_cf/vdp_analysis_results_October2024/FA_combined_fp_percent.csv")
N4_spir_cf_pfp <- read.csv("./IRC740H_visit1_spiral_cf/vdp_analysis_results_October2024/N4_combined_fp_percent.csv")

FA_spir_ctrl_pfp <- read.csv("./age_matched_spiral_healthy/vdp_analysis_results_October2024/FA_combined_fp_percent.csv")
N4_spir_ctrl_pfp <- read.csv("./age_matched_spiral_healthy/vdp_analysis_results_October2024/N4_combined_fp_percent.csv")

# #Reshape Data to Long Format for ggplot2
FA_cf_long <- FA_spir_cf_pfp %>%
  pivot_longer(cols = -Subject_id, names_to = "Analysis_Method", values_to = "pFP")
N4_cf_long <- N4_spir_cf_pfp %>%
  pivot_longer(cols = -Subject_id, names_to = "Analysis_Method", values_to = "pFP")

FA_ctrl_long <- FA_spir_ctrl_pfp %>%
  pivot_longer(cols = -Subject_id, names_to = "Analysis_Method", values_to = "pFP")
N4_ctrl_long <- N4_spir_ctrl_pfp %>%
  pivot_longer(cols = -Subject_id, names_to = "Analysis_Method", values_to = "pFP")
##To make false postive as percentage
FA_cf_long$pFP <- FA_cf_long$pFP * 100
N4_cf_long$pFP <- N4_cf_long$pFP * 100
FA_ctrl_long$pFP <- FA_ctrl_long$pFP * 100
N4_ctrl_long$pFP <- N4_ctrl_long$pFP * 100

##############################BAR Plots and statistical tests###################
#########################FA
cf_fa_bar <- create_barplot(FA_cf_long, "Analysis_Method", "pFP",
               expression(bold("Mean FP VDP"[FA]* " " *"(%)")),
               "Analysis Method", c(0, 30), cbPalette, "Subject_id", plt_order)
friedman.test(pFP ~ Analysis_Method | Subject_id, data = FA_cf_long)
frdAllPairsNemenyiTest(pFP ~ Analysis_Method | Subject_id, data = FA_cf_long)

ctrl_fa_bar <- create_barplot(FA_ctrl_long, "Analysis_Method", "pFP",
                            expression(bold("Mean FP VDP"[FA]* " " *"(%)")),
                            "Analysis Method", c(-7, 50), cbPalette, "Subject_id", plt_order)
friedman.test(pFP ~ Analysis_Method | Subject_id, data = FA_ctrl_long)
frdAllPairsNemenyiTest(pFP ~ Analysis_Method | Subject_id, data = FA_ctrl_long)

#########################N4: expression(bold("Mean FP VDP"[N4]* " " *"(%)"))
cf_n4_bar <- create_barplot(N4_cf_long, "Analysis_Method", "pFP",
                            "Mean FP VDP (%)",
               "Analysis Method", c(0, 30), cbPalette, "Subject_id", plt_order)
friedman.test(pFP ~ Analysis_Method | Subject_id, data = N4_cf_long)
frdAllPairsNemenyiTest(pFP ~ Analysis_Method | Subject_id, data = N4_cf_long)

ctrl_n4_bar <- create_barplot(N4_ctrl_long, "Analysis_Method", "pFP",
                              "Mean FP VDP (%)",
                            "Analysis Method", c(0, 10), cbPalette, "Subject_id", plt_order)
friedman.test(pFP ~ Analysis_Method | Subject_id, data = N4_ctrl_long)
frdAllPairsNemenyiTest(pFP ~ Analysis_Method | Subject_id, data = N4_ctrl_long)

# Save the plot as a png file in the specified directory
ggsave("./zR_plots_4abs/FP_defect_percent_barplt_cf_FA.png", plot = cf_fa_bar, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4abs/FP_defect_percent_barplt_cf_N4.png", plot = cf_n4_bar, width = 4.5, height = 3.7, dpi = 300)

ggsave("./zR_plots_4abs/FP_defect_percent_barplt_ctrl_FA.png", plot = ctrl_fa_bar, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4abs/FP_defect_percent_barplt_ctrl_N4.png", plot = ctrl_n4_bar, width = 4.5, height = 3.7, dpi = 300)
