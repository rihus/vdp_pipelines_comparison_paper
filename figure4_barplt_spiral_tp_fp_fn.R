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
FA_spir_cf_tp <- read.csv("./IRC740H_visit1_spiral_cf/vdp_analysis_results_October2024/FA_combined_TP_rdr_fractions.csv")
N4_spir_cf_tp <- read.csv("./IRC740H_visit1_spiral_cf/vdp_analysis_results_October2024/N4_combined_TP_rdr_fractions.csv")

FA_spir_cf_fp <- read.csv("./IRC740H_visit1_spiral_cf/vdp_analysis_results_October2024/FA_combined_FP_rdr_fractions.csv")
N4_spir_cf_fp <- read.csv("./IRC740H_visit1_spiral_cf/vdp_analysis_results_October2024/N4_combined_FP_rdr_fractions.csv")

FA_spir_cf_fn <- read.csv("./IRC740H_visit1_spiral_cf/vdp_analysis_results_October2024/FA_combined_FN_rdr_fractions.csv")
N4_spir_cf_fn <- read.csv("./IRC740H_visit1_spiral_cf/vdp_analysis_results_October2024/N4_combined_FN_rdr_fractions.csv")


# #Reshape Data to Long Format for ggplot2
FA_tp_long <- FA_spir_cf_tp %>%
  pivot_longer(cols = -Subject_id, names_to = "Analysis_Method", values_to = "TP")
N4_tp_long <- N4_spir_cf_tp %>%
  pivot_longer(cols = -Subject_id, names_to = "Analysis_Method", values_to = "TP")

FA_fp_long <- FA_spir_cf_fp %>%
  pivot_longer(cols = -Subject_id, names_to = "Analysis_Method", values_to = "FP")
N4_fp_long <- N4_spir_cf_fp %>%
  pivot_longer(cols = -Subject_id, names_to = "Analysis_Method", values_to = "FP")

FA_fn_long <- FA_spir_cf_fn %>%
  pivot_longer(cols = -Subject_id, names_to = "Analysis_Method", values_to = "FN")
N4_fn_long <- N4_spir_cf_fn %>%
  pivot_longer(cols = -Subject_id, names_to = "Analysis_Method", values_to = "FN")


##############################BAR Plots and statistical tests###################
#########################FA
fa_tp_bar <- create_barplot(FA_tp_long, "Analysis_Method", "TP",
               expression(bold("Mean TP"[FA]* " " *"Fraction")),
               "Analysis Method", c(-0.1, 2), cbPalette, "Subject_id", plt_order)
friedman.test(TP ~ Analysis_Method | Subject_id, data = FA_tp_long)
frdAllPairsNemenyiTest(TP ~ Analysis_Method | Subject_id, data = FA_tp_long)

fa_fp_bar <- create_barplot(FA_fp_long, "Analysis_Method", "FP",
                            expression(bold("Mean FP"[FA]* " " *"Fraction")),
                            "Analysis Method", c(-20, 120), cbPalette, "Subject_id", plt_order)
friedman.test(FP ~ Analysis_Method | Subject_id, data = FA_fp_long)
frdAllPairsNemenyiTest(FP ~ Analysis_Method | Subject_id, data = FA_fp_long)

fa_fn_bar <- create_barplot(FA_fn_long, "Analysis_Method", "FN",
                            expression(bold("Mean FN"[FA]* " " *"Fraction")),
                            "Analysis Method", c(-0.1, 2), cbPalette, "Subject_id", plt_order)
friedman.test(FN ~ Analysis_Method | Subject_id, data = FA_fn_long)
frdAllPairsNemenyiTest(FN ~ Analysis_Method | Subject_id, data = FA_fn_long)
#########################N4
n4_tp_bar <- create_barplot(N4_tp_long, "Analysis_Method", "TP",
               expression(bold("Mean TP"[N4]* " " *"Fraction")),
               "Analysis Method", c(-0.1, 2), cbPalette, "Subject_id", plt_order)
friedman.test(TP ~ Analysis_Method | Subject_id, data = N4_tp_long)
frdAllPairsNemenyiTest(TP ~ Analysis_Method | Subject_id, data = N4_tp_long)

n4_fp_bar <- create_barplot(N4_fp_long, "Analysis_Method", "FP",
                            expression(bold("Mean FP"[N4]* " " *"Fraction")),
                            "Analysis Method", c(-10, 70), cbPalette, "Subject_id", plt_order)
friedman.test(FP ~ Analysis_Method | Subject_id, data = N4_fp_long)
frdAllPairsNemenyiTest(FP ~ Analysis_Method | Subject_id, data = N4_fp_long)

n4_fn_bar <- create_barplot(N4_fn_long, "Analysis_Method", "FN",
                            expression(bold("Mean FN"[N4]* " " *"Fraction")),
                            "Analysis Method", c(-0.1, 2), cbPalette, "Subject_id", plt_order)
friedman.test(FN ~ Analysis_Method | Subject_id, data = N4_fn_long)
frdAllPairsNemenyiTest(FN ~ Analysis_Method | Subject_id, data = N4_fn_long)




# Save the plot as a png file in the specified directory
ggsave("./zR_plots_4abs/cf_barplt_spir_rdr_vs_FA_defects_TPs.png", plot = fa_tp_bar, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4abs/cf_barplt_spir_rdr_vs_N4_defects_TPs.png", plot = n4_tp_bar, width = 4.5, height = 3.7, dpi = 300)

ggsave("./zR_plots_4abs/cf_barplt_spir_rdr_vs_FA_defects_FPs.png", plot = fa_fp_bar, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4abs/cf_barplt_spir_rdr_vs_N4_defects_FPs.png", plot = n4_fp_bar, width = 4.5, height = 3.7, dpi = 300)

ggsave("./zR_plots_4abs/cf_barplt_spir_rdr_vs_FA_defects_FNs.png", plot = fa_fn_bar, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4abs/cf_barplt_spir_rdr_vs_N4_defects_FNs.png", plot = n4_fn_bar, width = 4.5, height = 3.7, dpi = 300)
