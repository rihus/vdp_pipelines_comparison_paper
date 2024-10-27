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
plt_order <- c("Thresholding","Hierarchical","Adaptive","Percentile","Median","Reader")
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
                line.size = 0.5, legend = "none", xlab = xlab) +
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
  ## Statistical test
  pval <- NULL
  if (tst == "wilcox") {
    pval <- data_in %>%
      wilcox_test(formula, paired = paired) %>%
      add_significance()
  } else if (tst == "ttest") {
    pval <- data_in %>%
      t_test(formula, paired = paired) %>%
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
  if (addp_eq == TRUE) {plabel = "p={scales::pvalue(p)}"}
  else {plabel = "p{scales::pvalue(p)}"}
  bxp_p <- fig_handle + stat_pvalue_manual(p_thresh, label = plabel,
                                           y.position = py_pos, label.size = 8,
                                           bracket.size = 0.8,
                                           tip.length = 0.03, vjust=-0.35)
  print(bxp_p)
  return(bxp_p)
}

################################################################################
# ##Load the data from both CSV files and arrange it to plot
N4_spir_cf_vdps <- read.csv("./IRC740H_visit1_spiral_cf/vdp_analysis_results_September2024/N4_combined_vdps.csv")
FA_spir_cf_vdps <- read.csv("./IRC740H_visit1_spiral_cf/vdp_analysis_results_September2024/FA_combined_vdps.csv")

##Reshape Data to Long Format for ggplot2
N4_long <- N4_spir_cf_vdps %>%
  pivot_longer(cols = -Subject_id, names_to = "AnalysisMethod", values_to = "VDP")
FA_long <- FA_spir_cf_vdps %>%
  pivot_longer(cols = -Subject_id, names_to = "AnalysisMethod", values_to = "VDP")

# ##Set the order the order of boxes in plots
# N4_long$AnalysisMethod <- factor(N4_long$AnalysisMethod, levels = plt_order)
# FA_long$AnalysisMethod <- factor(FA_long$AnalysisMethod, levels = plt_order)

##Select data for only "Reader" and one analysis plot
N4_rdr_thrsh <- N4_long %>% filter(AnalysisMethod %in% c("Reader", "Thresholding"))
N4_rdr_hrar <- N4_long %>% filter(AnalysisMethod %in% c("Reader", "Hierarchical"))
N4_rdr_adpt <- N4_long %>% filter(AnalysisMethod %in% c("Reader", "Adaptive"))
N4_rdr_prcnt <- N4_long %>% filter(AnalysisMethod %in% c("Reader", "Percentile"))
N4_rdr_medn <- N4_long %>% filter(AnalysisMethod %in% c("Reader", "Median"))

FA_rdr_thrsh <- FA_long %>% filter(AnalysisMethod %in% c("Reader", "Thresholding"))
FA_rdr_hrar <- FA_long %>% filter(AnalysisMethod %in% c("Reader", "Hierarchical"))
FA_rdr_adpt <- FA_long %>% filter(AnalysisMethod %in% c("Reader", "Adaptive"))
FA_rdr_prcnt <- FA_long %>% filter(AnalysisMethod %in% c("Reader", "Percentile"))
FA_rdr_medn <- FA_long %>% filter(AnalysisMethod %in% c("Reader", "Median"))

##Define y-label of the plot
ylabeln4 <- expression(bold(VDP[Spiral-N4]* "(%)"))
ylabelfa <- expression(bold(VDP[Spiral-FA]* "(%)"))
##############################################
thresh_bxp_n4 <- connected_bxp(N4_rdr_thrsh, "AnalysisMethod","VDP", "Subject_id",
                               c(0, 40), cbPalette[c(6, 1)], "", ylabeln4)
hrar_bxp_n4 <- connected_bxp(N4_rdr_hrar, "AnalysisMethod","VDP", "Subject_id",
                                c(0, 40), cbPalette[c(6, 2)], "", ylabeln4)
adpt_bxp_n4 <- connected_bxp(N4_rdr_adpt, "AnalysisMethod","VDP", "Subject_id",
                                c(0, 40), cbPalette[c(6, 3)], "", ylabeln4)
prcnt_bxp_n4 <- connected_bxp(N4_rdr_prcnt, "AnalysisMethod","VDP", "Subject_id",
                                c(0, 40), cbPalette[c(6, 4)], "", ylabeln4)
medn_bxp_n4 <- connected_bxp(N4_rdr_medn, "AnalysisMethod","VDP", "Subject_id",
                                c(0, 40), cbPalette[c(6, 5)], "", ylabeln4)
thresh_bxp_fa <- connected_bxp(FA_rdr_thrsh, "AnalysisMethod","VDP", "Subject_id",
                               c(0, 40), cbPalette[c(6, 1)], "", ylabelfa)
hrar_bxp_fa <- connected_bxp(FA_rdr_hrar, "AnalysisMethod","VDP", "Subject_id",
                             c(0, 40), cbPalette[c(6, 2)], "", ylabelfa)
adpt_bxp_fa <- connected_bxp(FA_rdr_adpt, "AnalysisMethod","VDP", "Subject_id",
                             c(0, 40), cbPalette[c(6, 3)], "", ylabelfa)
prcnt_bxp_fa <- connected_bxp(FA_rdr_prcnt, "AnalysisMethod","VDP", "Subject_id",
                              c(0, 40), cbPalette[c(6, 4)], "", ylabelfa)
medn_bxp_fa <- connected_bxp(FA_rdr_medn, "AnalysisMethod","VDP", "Subject_id",
                             c(0, 40), cbPalette[c(6, 5)], "", ylabelfa)

##Calculate and add p-value to the plot
thresh_bxp_n4_p <- calc_add_p(N4_rdr_thrsh, "AnalysisMethod", "VDP", thresh_bxp_n4,
                              35, tst = "wilcox", paird = TRUE, addp_eq=FALSE)
hrar_bxp_n4_p <- calc_add_p(N4_rdr_hrar, "AnalysisMethod", "VDP", hrar_bxp_n4,
                            35, tst = "wilcox", paird = TRUE, addp_eq=TRUE)
adpt_bxp_n4_p <- calc_add_p(N4_rdr_adpt, "AnalysisMethod", "VDP", adpt_bxp_n4,
                            35, tst = "wilcox", paird = TRUE, addp_eq=TRUE)
prcnt_bxp_n4_p <- calc_add_p(N4_rdr_prcnt, "AnalysisMethod", "VDP", prcnt_bxp_n4,
                             35, tst = "wilcox", paird = TRUE, addp_eq=FALSE)
medn_bxp_n4_p <- calc_add_p(N4_rdr_medn, "AnalysisMethod", "VDP", medn_bxp_n4,
                            35, tst = "wilcox", paird = TRUE, addp_eq=FALSE)

thresh_bxp_fa_p <- calc_add_p(FA_rdr_thrsh, "AnalysisMethod", "VDP", thresh_bxp_fa,
                              35, tst = "wilcox", paird = TRUE, addp_eq=TRUE)
hrar_bxp_fa_p <- calc_add_p(FA_rdr_hrar, "AnalysisMethod", "VDP", hrar_bxp_fa,
                            35, tst = "wilcox", paird = TRUE, addp_eq=TRUE)
adpt_bxp_fa_p <- calc_add_p(FA_rdr_adpt, "AnalysisMethod", "VDP", adpt_bxp_fa,
                            35, tst = "wilcox", paird = TRUE, addp_eq=FALSE)
prcnt_bxp_fa_p <- calc_add_p(FA_rdr_prcnt, "AnalysisMethod", "VDP", prcnt_bxp_fa,
                             35, tst = "wilcox", paird = TRUE, addp_eq=TRUE)
medn_bxp_fa_p <- calc_add_p(FA_rdr_medn, "AnalysisMethod", "VDP", medn_bxp_fa,
                            35, tst = "wilcox", paird = TRUE, addp_eq=FALSE)

##Save figures, with/without p value
ggsave("./zR_plots_4abs/spir_vdp_rdr_thresh_bxp_n4.png", plot = thresh_bxp_n4, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4abs/spir_vdp_rdr_thresh_bxp_n4_p.png", plot = thresh_bxp_n4_p, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4abs/spir_vdp_rdr_hrar_bxp_n4.png", plot = hrar_bxp_n4, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4abs/spir_vdp_rdr_hrar_bxp_n4_p.png", plot = hrar_bxp_n4_p, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4abs/spir_vdp_rdr_adpt_bxp_n4.png", plot = adpt_bxp_n4, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4abs/spir_vdp_rdr_adpt_bxp_n4_p.png", plot = adpt_bxp_n4_p, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4abs/spir_vdp_rdr_prcnt_bxp_n4.png", plot = prcnt_bxp_n4, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4abs/spir_vdp_rdr_prcnt_bxp_n4_p.png", plot = prcnt_bxp_n4_p, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4abs/spir_vdp_rdr_medn_bxp_n4.png", plot = medn_bxp_n4, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4abs/spir_vdp_rdr_medn_bxp_n4_p.png", plot = medn_bxp_n4_p, width = 4.5, height = 3.7, dpi = 300)

ggsave("./zR_plots_4abs/spir_vdp_rdr_thresh_bxp_fa.png", plot = thresh_bxp_fa, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4abs/spir_vdp_rdr_thresh_bxp_fa_p.png", plot = thresh_bxp_fa_p, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4abs/spir_vdp_rdr_hrar_bxp_fa.png", plot = hrar_bxp_fa, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4abs/spir_vdp_rdr_hrar_bxp_fa_p.png", plot = hrar_bxp_fa_p, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4abs/spir_vdp_rdr_adpt_bxp_fa.png", plot = adpt_bxp_fa, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4abs/spir_vdp_rdr_adpt_bxp_fa_p.png", plot = adpt_bxp_fa_p, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4abs/spir_vdp_rdr_prcnt_bxp_fa.png", plot = prcnt_bxp_fa, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4abs/spir_vdp_rdr_prcnt_bxp_fa_p.png", plot = prcnt_bxp_fa_p, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4abs/spir_vdp_rdr_medn_bxp_fa.png", plot = medn_bxp_fa, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4abs/spir_vdp_rdr_medn_bxp_fa_p.png", plot = medn_bxp_fa_p, width = 4.5, height = 3.7, dpi = 300)

