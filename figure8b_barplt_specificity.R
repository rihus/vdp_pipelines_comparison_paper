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
cbPalette <- c("#882255",  "#CC6677", "#332288", "#6699CC", "#DDCC77", "#999933",
               "#AA4499", "#661100", "#44AA99", "#117733", "#888888", "#88CCEE") #
##Define order of the data
plt_order <- c("Reader", "Hierarchical","Adaptive","Percentile",
               "Mean","Thresholding", "Median")

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
# ##Load the specificity data from CSV files
FA_cf_specificity <- read.csv("./IRC740H_visit1_spiral_cf/vdp_analysis_results_May2025/FA_combined_specificity.csv")
N4_cf_specificity <- read.csv("./IRC740H_visit1_spiral_cf/vdp_analysis_results_May2025/N4_combined_specificity.csv")

# #Reshape Data to Long Format for ggplot2
FA_cf_long <- FA_cf_specificity %>%
  pivot_longer(cols = -Subject_id, names_to = "Analysis_Method", values_to = "Specificity") %>%
  filter(!is.na(Specificity)) %>%
  mutate(Analysis_Method = factor(Analysis_Method, levels = c("Hierarchical", "Adaptive", "Percentile", "Mean", "Thresholding", "Median")))

N4_cf_long <- N4_cf_specificity %>%
  pivot_longer(cols = -Subject_id, names_to = "Analysis_Method", values_to = "Specificity") %>%
  filter(!is.na(Specificity)) %>%
  mutate(Analysis_Method = factor(Analysis_Method, levels = c("Hierarchical", "Adaptive", "Percentile", "Mean", "Thresholding", "Median")))

##To make false postive as percentage
FA_cf_long$Specificity <- FA_cf_long$Specificity * 100
N4_cf_long$Specificity <- N4_cf_long$Specificity * 100

##############################BAR Plots and statistical tests###################
#########################FA
cf_fa_bar <- create_barplot(FA_cf_long, "Analysis_Method", "Specificity",
               expression(bold("Specificity"[FA]* " " *"(%)")),
               "Analysis Method", c(0, 160), cbPalette, "Subject_id", plt_order)
friedman.test(Specificity ~ Analysis_Method | Subject_id, data = FA_cf_long)
frdAllPairsNemenyiTest(Specificity ~ Analysis_Method | Subject_id, data = FA_cf_long)

#########################N4: expression(bold("Mean FP VDP"[N4]* " " *"(%)"))
cf_n4_bar <- create_barplot(N4_cf_long, "Analysis_Method", "Specificity",
                            expression(bold("Specificity"[N4]* " " *"(%)")),
               "Analysis Method", c(0, 160), cbPalette, "Subject_id", plt_order)
friedman.test(Specificity ~ Analysis_Method | Subject_id, data = N4_cf_long)
frdAllPairsNemenyiTest(Specificity ~ Analysis_Method | Subject_id, data = N4_cf_long)

# Save the plot as a png file in the specified directory
ggsave("./zR_plots_4ppr/cf_Specificity_barplt_FA.png", plot = cf_fa_bar, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/cf_Specificity_barplt_N4.png", plot = cf_n4_bar, width = 4.5, height = 3.7, dpi = 300)


