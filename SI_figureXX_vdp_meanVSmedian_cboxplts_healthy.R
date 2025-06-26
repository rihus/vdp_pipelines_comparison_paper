#load ggplot2
library(ggplot2)
library(ggpubr)
library(blandr)
library(dplyr)
library(tidyr)
library(rstatix)

# Set the working directory to the path where your CSV files are located
setwd("C:/Users/HUSDQ4/OneDrive - cchmc/cincy_work/all_projects_data_work/vdp_analysis/analysis_comparisons")

##Colorblind friendly palette:
cbPalette <- c("#888888", "#CC6677", "#882255", "#332288", "#6699CC", "#DDCC77",
               "#999933", "#AA4499", "#661100", "#44AA99", "#117733", "#88CCEE")

##Define order of the data
plt_order <- c("Reader", "Adaptive","Hierarchical","Thresholding","Percentile",
               "Mean","Median", "Mean_GLB")
################################################################################
###Function for plotting connected boxplot
connected_bxp <- function(data_in, x_var, y_var, id = "Subject_id", ylim,
                          palette = NULL, xlab = "", ylab = NULL) {
  require(ggpubr)
  # Check if palette is provided
  if (is.null(palette) || is.null(ylab)) {
    palette <- c("#999999", "#555555")
    ylab <- y_var}
  
  # Create the ggpaired plot
  bxp <- ggpaired(data_in, x = x_var, y = y_var, id = id,
                  fill = x_var, palette = palette,
                  width = 0.5, ylim = ylim, line.color = "#000000",
                  line.size = 0.5, point.size = 2, legend = "none", xlab = xlab) +
    ylab(ylab) +
    theme(panel.border = element_rect(color = "#000000", fill = NA, linewidth = 1),
          axis.text = element_text(size = 22, color = "#000000", face = "bold"),
          axis.title = element_text(size = 22, color = "#000000", face = "bold"),
          axis.line.x = element_line(linewidth = 1), 
          axis.line.y = element_line(linewidth = 1))
  print(bxp)
  return(bxp)
}

calc_pval <- function(data_in, x_var, y_var, tst = "wilcox", paired = TRUE) {
  formula <- as.formula(paste(y_var, "~", x_var))
  ##Mean +- SD
  mean_val <- aggregate(formula, data_in, mean)
  print("Mean")
  print(mean_val)
  sd_val <- aggregate(formula, data_in, sd)
  print("Standard-deviation")
  print(sd_val)
  ## Statistical test
  pval <- NULL
  if (tst == "wilcox") {
    pval <- data_in %>%
      wilcox_test(formula, paired = paired) %>%
      adjust_pvalue(method = 'bonferroni') %>%
      add_significance()
  } else if (tst == "ttest") {
    pval <- data_in %>%
      t_test(formula, paired = paired) %>%
      adjust_pvalue(method = 'bonferroni') %>%
      add_significance()
  } else {
    print("For tst only Wilcoxon (wilcox) or T-test (ttest) are accepted")
  }
  print(pval)
  return(pval)
}

calc_add_p <- function(data_in, x_var, y_var, fig_handle, py_pos,
                       tst = "wilcox", paird = TRUE, addp_eq=FALSE) {
  p_thresh <- calc_pval(data_in, x_var, y_var, tst, paird)
  p_thresh <- p_thresh %>% add_xy_position(x = x_var)
  if (addp_eq == TRUE) {plabel = "P={scales::pvalue(p, accuracy = 0.001)}"}
  else {plabel = "P{scales::pvalue(p, accuracy = 0.0001)}"}
  bxp_p <- fig_handle + stat_pvalue_manual(p_thresh, label = plabel,#"P = {p.adj}",
                                           y.position = py_pos, label.size = 8,
                                           bracket.size = 0.8,
                                           tip.length = 0.03, vjust=-0.35)
  print(bxp_p)
  return(bxp_p)
}

################################################################################

# ##Load the data from both CSV files and arrange it to plot
N4_vdps_mdn <- read.csv("./age_matched_spiral_healthy/vdp_analysis_results_December2024/N4_combined_vdps.csv")
N4_vdps_mean <- read.csv("./age_matched_spiral_healthy/vdp_analysis_results_May2025/N4_combined_vdps.csv")

FA_vdps_mdn <- read.csv("./age_matched_spiral_healthy/vdp_analysis_results_December2024/FA_combined_vdps.csv")
FA_vdps_mean <- read.csv("./age_matched_spiral_healthy/vdp_analysis_results_May2025/FA_combined_vdps.csv")

names(N4_vdps_mean)[names(N4_vdps_mean) == "Median"] <- "Mean_GLB"
names(FA_vdps_mean)[names(FA_vdps_mean) == "Median"] <- "Mean_GLB"

# #Combine the data
combined_N4 <- data.frame(VDP = c(N4_vdps_mdn$Median, N4_vdps_mean$Mean_GLB),
  AnalysisMethod = c(rep("Median", nrow(N4_vdps_mdn)),rep("Mean_GLB", nrow(N4_vdps_mean))))

combined_FA <- data.frame(VDP = c(FA_vdps_mdn$Median, FA_vdps_mean$Mean_GLB),
  AnalysisMethod = c(rep("Median", nrow(FA_vdps_mdn)), rep("Mean_GLB", nrow(FA_vdps_mean))))

# #Combine paired data into one dataframe
paired_N4_data <- na.omit(data.frame(Subject_id = N4_vdps_mdn$Subject_id,
  Median = N4_vdps_mdn$Median, Mean_GLB = N4_vdps_mean$Mean_GLB))

paired_FA_data <- na.omit(data.frame(Subject_id = FA_vdps_mdn$Subject_id,
                                     Median = FA_vdps_mdn$Median, Mean_GLB = FA_vdps_mean$Mean_GLB))

long_N4 <- pivot_longer(paired_N4_data, cols = c("Median", "Mean_GLB"),
                          names_to = "AnalysisMethod", values_to = "VDP")

long_FA <- pivot_longer(paired_FA_data, cols = c("Median", "Mean_GLB"),
                        names_to = "AnalysisMethod", values_to = "VDP")

ylabeln4 <- expression(bold(VDP[N4]*" " *"(%)"))
ylabelfa <- expression(bold(VDP[FA]*" " *"(%)"))

medn_bxp_n4 <- connected_bxp(long_N4, "AnalysisMethod","VDP", "Subject_id",
                             c(0, 7.5), cbPalette[c(7, 8)], "", ylabeln4)

medn_bxp_n4_p <- calc_add_p(long_N4, "AnalysisMethod", "VDP", medn_bxp_n4,
                            6, tst = "wilcox", paird = TRUE, addp_eq=FALSE)

#
medn_bxp_fa <- connected_bxp(long_FA, "AnalysisMethod","VDP", "Subject_id",
                             c(0, 13), cbPalette[c(7, 8)], "", ylabelfa)

medn_bxp_fa_p <- calc_add_p(long_FA, "AnalysisMethod", "VDP", medn_bxp_fa,
                            11.5, tst = "wilcox", paird = TRUE, addp_eq=TRUE)

ggsave("./zR_plots_4ppr/healthy_bxp_medn_mean_N4_vdp.png", plot = medn_bxp_n4, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/healthy_bxp_medn_mean_N4_vdp_p.png", plot = medn_bxp_n4_p, width = 4.5, height = 3.7, dpi = 300)

ggsave("./zR_plots_4ppr/healthy_bxp_medn_mean_FA_vdp.png", plot = medn_bxp_fa, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/healthy_bxp_medn_mean_FA_vdp_p.png", plot = medn_bxp_fa_p, width = 4.5, height = 3.7, dpi = 300)
