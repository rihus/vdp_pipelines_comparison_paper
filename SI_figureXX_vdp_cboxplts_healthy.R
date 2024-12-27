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
cbPalette <- c("#888888", "#CC6677", "#882255", "#88CCEE", "#332288", "#6699CC",
               "#999933", "#AA4499", "#DDCC77", "#661100", "#44AA99", "#117733")

##Define order of the data
plt_order <- c("Reader", "Adaptive","Hierarchical","Thresholding","Percentile",
               "Mean","Median")
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
  else {plabel = "P{scales::pvalue(p, accuracy = 0.001)}"}
  bxp_p <- fig_handle + stat_pvalue_manual(p_thresh, label = plabel,#"P = {p.adj}",
                                           y.position = py_pos, label.size = 8,
                                           bracket.size = 0.8,
                                           tip.length = 0.03, vjust=-0.35)
  print(bxp_p)
  return(bxp_p)
}

################################################################################
# ##Load the data from both CSV files and arrange it to plot
N4_spir_ctrl_vdps <- read.csv("./age_matched_spiral_healthy/vdp_analysis_results_December2024/N4_combined_vdps.csv")
FA_spir_ctrl_vdps <- read.csv("./age_matched_spiral_healthy/vdp_analysis_results_December2024/FA_combined_vdps.csv")

##Reshape Data to Long Format for ggplot2
N4_long <- N4_spir_ctrl_vdps %>%
  pivot_longer(cols = -Subject_id, names_to = "AnalysisMethod", values_to = "VDP")
FA_long <- FA_spir_ctrl_vdps %>%
  pivot_longer(cols = -Subject_id, names_to = "AnalysisMethod", values_to = "VDP")

# ##Set the order the order of boxes in plots
# N4_long$AnalysisMethod <- factor(N4_long$AnalysisMethod, levels = plt_order)
# FA_long$AnalysisMethod <- factor(FA_long$AnalysisMethod, levels = plt_order)

##Select data for only "Reader" and one analysis plot
N4_rdr_adpt <- N4_long %>% filter(AnalysisMethod %in% c("Reader", "Adaptive"))
N4_rdr_hrar <- N4_long %>% filter(AnalysisMethod %in% c("Reader", "Hierarchical"))
N4_rdr_thrsh <- N4_long %>% filter(AnalysisMethod %in% c("Reader", "Thresholding"))
N4_rdr_prcnt <- N4_long %>% filter(AnalysisMethod %in% c("Reader", "Percentile"))
N4_rdr_mean <- N4_long %>% filter(AnalysisMethod %in% c("Reader", "Mean"))
N4_rdr_medn <- N4_long %>% filter(AnalysisMethod %in% c("Reader", "Median"))

FA_rdr_adpt <- FA_long %>% filter(AnalysisMethod %in% c("Reader", "Adaptive"))
FA_rdr_hrar <- FA_long %>% filter(AnalysisMethod %in% c("Reader", "Hierarchical"))
FA_rdr_thrsh <- FA_long %>% filter(AnalysisMethod %in% c("Reader", "Thresholding"))
FA_rdr_prcnt <- FA_long %>% filter(AnalysisMethod %in% c("Reader", "Percentile"))
FA_rdr_mean <- FA_long %>% filter(AnalysisMethod %in% c("Reader", "Mean"))
FA_rdr_medn <- FA_long %>% filter(AnalysisMethod %in% c("Reader", "Median"))

##Define y-label of the plot
# ylabeln4 <- "VDP (%)"
# ylabelfa <- "VDP (%)"
ylabeln4 <- expression(bold(VDP[N4]*" " *"(%)")) # [Spiral]
ylabelfa <- expression(bold(VDP[FA]*" " *"(%)"))
##############################################
adpt_bxp_n4 <- connected_bxp(N4_rdr_adpt, "AnalysisMethod","VDP", "Subject_id",
                             c(0, 5), cbPalette[c(1, 2)], "", ylabeln4)
hrar_bxp_n4 <- connected_bxp(N4_rdr_hrar, "AnalysisMethod","VDP", "Subject_id",
                             c(0, 5), cbPalette[c(1, 3)], "", ylabeln4)
thresh_bxp_n4 <- connected_bxp(N4_rdr_thrsh, "AnalysisMethod","VDP", "Subject_id",
                               c(0, 5), cbPalette[c(1, 4)], "", ylabeln4)
prcnt_bxp_n4 <- connected_bxp(N4_rdr_prcnt, "AnalysisMethod","VDP", "Subject_id",
                                c(0, 5), cbPalette[c(1, 5)], "", ylabeln4)
mean_bxp_n4 <- connected_bxp(N4_rdr_mean, "AnalysisMethod","VDP", "Subject_id",
                              c(0, 1), cbPalette[c(1, 6)], "", ylabeln4)
medn_bxp_n4 <- connected_bxp(N4_rdr_medn, "AnalysisMethod","VDP", "Subject_id",
                                c(0, 5), cbPalette[c(1, 7)], "", ylabeln4)

adpt_bxp_fa <- connected_bxp(FA_rdr_adpt, "AnalysisMethod","VDP", "Subject_id",
                             c(0, 60), cbPalette[c(1, 2)], "", ylabelfa)
hrar_bxp_fa <- connected_bxp(FA_rdr_hrar, "AnalysisMethod","VDP", "Subject_id",
                             c(0, 22), cbPalette[c(1, 3)], "", ylabelfa)
thresh_bxp_fa <- connected_bxp(FA_rdr_thrsh, "AnalysisMethod","VDP", "Subject_id",
                               c(0, 10), cbPalette[c(1, 4)], "", ylabelfa)
prcnt_bxp_fa <- connected_bxp(FA_rdr_prcnt, "AnalysisMethod","VDP", "Subject_id",
                              c(0, 20), cbPalette[c(1, 5)], "", ylabelfa)
mean_bxp_fa <- connected_bxp(FA_rdr_mean, "AnalysisMethod","VDP", "Subject_id",
                              c(0, 10), cbPalette[c(1, 6)], "", ylabelfa)
medn_bxp_fa <- connected_bxp(FA_rdr_medn, "AnalysisMethod","VDP", "Subject_id",
                             c(0, 15), cbPalette[c(1, 7)], "", ylabelfa)

##Calculate and add p-value to the plot
adpt_bxp_n4_p <- calc_add_p(N4_rdr_adpt, "AnalysisMethod", "VDP", adpt_bxp_n4,
                            4.3, tst = "wilcox", paird = TRUE, addp_eq=FALSE)
hrar_bxp_n4_p <- calc_add_p(N4_rdr_hrar, "AnalysisMethod", "VDP", hrar_bxp_n4,
                            4.3, tst = "wilcox", paird = TRUE, addp_eq=FALSE)
thresh_bxp_n4_p <- calc_add_p(N4_rdr_thrsh, "AnalysisMethod", "VDP", thresh_bxp_n4,
                            4.3, tst = "wilcox", paird = TRUE, addp_eq=FALSE)
prcnt_bxp_n4_p <- calc_add_p(N4_rdr_prcnt, "AnalysisMethod", "VDP", prcnt_bxp_n4,
                            4.3, tst = "wilcox", paird = TRUE, addp_eq=FALSE)
mean_bxp_n4_p <- calc_add_p(N4_rdr_mean, "AnalysisMethod", "VDP", mean_bxp_n4,
                            0.85, tst = "wilcox", paird = TRUE, addp_eq=FALSE)
medn_bxp_n4_p <- calc_add_p(N4_rdr_medn, "AnalysisMethod", "VDP", medn_bxp_n4,
                            4.3, tst = "wilcox", paird = TRUE, addp_eq=FALSE)

adpt_bxp_fa_p <- calc_add_p(FA_rdr_adpt, "AnalysisMethod", "VDP", adpt_bxp_fa,
                            53, tst = "wilcox", paird = TRUE, addp_eq=FALSE)
hrar_bxp_fa_p <- calc_add_p(FA_rdr_hrar, "AnalysisMethod", "VDP", hrar_bxp_fa,
                            19, tst = "wilcox", paird = TRUE, addp_eq=FALSE)
thresh_bxp_fa_p <- calc_add_p(FA_rdr_thrsh, "AnalysisMethod", "VDP", thresh_bxp_fa,
                              8.5, tst = "wilcox", paird = TRUE, addp_eq=FALSE)
prcnt_bxp_fa_p <- calc_add_p(FA_rdr_prcnt, "AnalysisMethod", "VDP", prcnt_bxp_fa,
                             17, tst = "wilcox", paird = TRUE, addp_eq=FALSE)
mean_bxp_fa_p <- calc_add_p(FA_rdr_mean, "AnalysisMethod", "VDP", mean_bxp_fa,
                             8, tst = "wilcox", paird = TRUE, addp_eq=FALSE)
medn_bxp_fa_p <- calc_add_p(FA_rdr_medn, "AnalysisMethod", "VDP", medn_bxp_fa,
                            13, tst = "wilcox", paird = TRUE, addp_eq=FALSE)

##Save figures, with/without p value
ggsave("./zR_plots_4ppr/ctrl_bxp_adpt_spir_N4_rdr_vdp.png", plot = adpt_bxp_n4, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/ctrl_bxp_adpt_spir_N4_rdr_vdp_p.png", plot = adpt_bxp_n4_p, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/ctrl_bxp_hrar_spir_N4_rdr_vdp.png", plot = hrar_bxp_n4, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/ctrl_bxp_hrar_spir_N4_rdr_vdp_p.png", plot = hrar_bxp_n4_p, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/ctrl_bxp_thresh_spir_N4_rdr_vdp.png", plot = thresh_bxp_n4, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/ctrl_bxp_thresh_spir_N4_rdr_vdp_p.png", plot = thresh_bxp_n4_p, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/ctrl_bxp_prcnt_spir_N4_rdr_vdp.png", plot = prcnt_bxp_n4, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/ctrl_bxp_prcnt_spir_N4_rdr_vdp_p.png", plot = prcnt_bxp_n4_p, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/ctrl_bxp_mean_spir_N4_rdr_vdp.png", plot = mean_bxp_n4, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/ctrl_bxp_mean_spir_N4_rdr_vdp_p.png", plot = mean_bxp_n4_p, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/ctrl_bxp_medn_spir_N4_rdr_vdp.png", plot = medn_bxp_n4, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/ctrl_bxp_medn_spir_N4_rdr_vdp_p.png", plot = medn_bxp_n4_p, width = 4.5, height = 3.7, dpi = 300)

ggsave("./zR_plots_4ppr/ctrl_bxp_adpt_spir_FA_rdr_vdp.png", plot = adpt_bxp_fa, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/ctrl_bxp_adpt_spir_FA_rdr_vdp_p.png", plot = adpt_bxp_fa_p, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/ctrl_bxp_hrar_spir_FA_rdr_vdp.png", plot = hrar_bxp_fa, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/ctrl_bxp_hrar_spir_FA_rdr_vdp_p.png", plot = hrar_bxp_fa_p, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/ctrl_bxp_thresh_spir_FA_rdr_vdp.png", plot = thresh_bxp_fa, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/ctrl_bxp_thresh_spir_FA_rdr_vdp_p.png", plot = thresh_bxp_fa_p, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/ctrl_bxp_prcnt_spir_FA_rdr_vdp.png", plot = prcnt_bxp_fa, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/ctrl_bxp_prcnt_spir_FA_rdr_vdp_p.png", plot = prcnt_bxp_fa_p, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/ctrl_bxp_mean_spir_FA_rdr_vdp.png", plot = mean_bxp_fa, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/ctrl_bxp_mean_spir_FA_rdr_vdp_p.png", plot = mean_bxp_fa_p, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/ctrl_bxp_medn_spir_FA_rdr_vdp.png", plot = medn_bxp_fa, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/ctrl_bxp_medn_spir_FA_rdr_vdp_p.png", plot = medn_bxp_fa_p, width = 4.5, height = 3.7, dpi = 300)

##
