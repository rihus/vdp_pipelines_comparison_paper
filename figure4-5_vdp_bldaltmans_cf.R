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
cbPalette <- c("#888888", "#882255", "#CC6677", "#332288", "#6699CC", "#DDCC77",
               "#999933", "#AA4499", "#661100", "#44AA99", "#117733", "#88CCEE")

# ##Order of the colors for data
# order <- c("Reader", "Hierarchical","Adaptive","Percentile","Mean","Thresholding",
#            "Median")
################################################################################

BlandAltman <- function(main_df, column1, column2, us1, us2, x_lim, y_lim,
                        dot_clr = "#000000", l_pos=1, u_pos=1.0) {
  main_df$mean_val <- rowMeans(main_df[, c(column1, column2)])
  main_df$diff_val <- as.numeric(main_df[[column1]]) - as.numeric(main_df[[column2]])
  mean_diff <- mean(main_df$diff_val)
  # Lower and Upper 95% confidence interval limits
  l_lim <- mean_diff - 1.96 * sd(main_df$diff_val)
  u_lim <- mean_diff + 1.96 * sd(main_df$diff_val)
  # Insert x-axis and y-axis labels
  y_label <- bquote(bold(VDP[.(us1)] - VDP[.(us2)] ~ "(%)"))
  x_label <- bquote(bold((VDP[.(us1)] + VDP[.(us2)]) / 2 ~ "(%)"))
  
  # Bland-Altman plot
  theme_set(theme_bw())
  baltman <- ggplot(main_df, aes(x = mean_val, y = diff_val)) +
    xlim(x_lim[1], x_lim[2]) +
    ylim(y_lim[1], y_lim[2]) +
    geom_point(color = dot_clr, size=3, alpha=1) +
    geom_hline(yintercept = mean_diff, linewidth = 1) +
    geom_hline(yintercept = l_lim, color = "#000000", linetype = 4, linewidth = 1.25) +
    geom_hline(yintercept = u_lim, color = "#000000", linetype = 4, linewidth = 1.25) +
    ylab(y_label) + xlab(x_label) +
    guides(color = "none") +
    theme(axis.title.x = element_text(vjust = 0, size = 22, face = "bold"),
          axis.title.y = element_text(vjust = 0, size = 22, face = "bold"),
          axis.text = element_text(size = 22, face = "bold", color = "#000000"))
  print(baltman)
  # Adding labels
  baltman_labld <- baltman +
    annotate("text", x = 0.8 * x_lim[2], y = mean_diff,
             label = sprintf("Mean=%.3f", mean_diff), size = 5) +
    annotate("text", x = 0.78 * x_lim[2], y = u_pos * u_lim,
             label = sprintf("+1.96xSD=%.3f", u_lim), size = 5) +
    annotate("text", x = 0.78 * x_lim[2], y = l_pos * l_lim,
             label = sprintf("-1.96xSD=%.3f", l_lim), size = 5)
  print(baltman_labld)
  return(list(baltman, baltman_labld))
}


################################################################################
# ##Load the data from both CSV files and arrange it to plot
N4_spir_cf_vdps <- read.csv("./IRC740H_visit1_spiral_cf/vdp_analysis_results_May2025/N4_combined_vdps.csv")
FA_spir_cf_vdps <- read.csv("./IRC740H_visit1_spiral_cf/vdp_analysis_results_May2025/FA_combined_vdps.csv")

##Bland-Altman plots - FA_corrected
fa_hrar <- BlandAltman(main_df = FA_spir_cf_vdps, column1 = "Reader",
                       column2 = "Hierarchical", us1 = "Reader", us2 = "FA",
                       x_lim = c(0, 30), y_lim = c(-20, 12),
                       dot_clr = cbPalette[2], l_pos=1.1, u_pos = 1.5)

fa_adpt <- BlandAltman(main_df = FA_spir_cf_vdps, column1 = "Reader",
                       column2 = "Adaptive", us1 = "Reader", us2 = "FA",
                       x_lim = c(0, 30), y_lim = c(-20, 12),
                       dot_clr = cbPalette[3], l_pos=1.1, u_pos = 1.5)

fa_prcnt <- BlandAltman(main_df = FA_spir_cf_vdps, column1 = "Reader",
                       column2 = "Percentile", us1 = "Reader", us2 = "FA",
                       x_lim = c(0, 30), y_lim = c(-20, 12),
                       dot_clr = cbPalette[4], l_pos=1.1, u_pos = 1.5)

fa_mean <- BlandAltman(main_df = FA_spir_cf_vdps, column1 = "Reader",
                        column2 = "Mean", us1 = "Reader", us2 = "FA",
                        x_lim = c(0, 30), y_lim = c(-20, 12),
                        dot_clr = cbPalette[5], l_pos=1.1, u_pos = 1.5)

fa_thrsh <- BlandAltman(main_df = FA_spir_cf_vdps, column1 = "Reader",
                       column2 = "Thresholding", us1 = "Reader", us2 = "FA",
                       x_lim = c(0, 30), y_lim = c(-20, 12),
                       dot_clr = cbPalette[6], l_pos=1.1, u_pos = 1.5)

fa_mdn <- BlandAltman(main_df = FA_spir_cf_vdps, column1 = "Reader",
                        column2 = "Median", us1 = "Reader", us2 = "FA",
                        x_lim = c(0, 30), y_lim = c(-20, 12),
                        dot_clr = cbPalette[7], l_pos=1.1, u_pos = 1.5)

##Bland-Altman plots - N4_corrected
n4_hrar <- BlandAltman(main_df = N4_spir_cf_vdps, column1 = "Reader",
                       column2 = "Hierarchical", us1 = "Reader", us2 = "N4",
                       x_lim = c(0, 30), y_lim = c(-20, 12),
                       dot_clr = cbPalette[2], l_pos=1.1, u_pos = 1.1)

n4_adpt <- BlandAltman(main_df = N4_spir_cf_vdps, column1 = "Reader",
                       column2 = "Adaptive", us1 = "Reader", us2 = "N4",
                       x_lim = c(0, 30), y_lim = c(-20, 12),
                       dot_clr = cbPalette[3], l_pos=1.1, u_pos = 1.5)

n4_prcnt <- BlandAltman(main_df = N4_spir_cf_vdps, column1 = "Reader",
                        column2 = "Percentile", us1 = "Reader", us2 = "N4",
                        x_lim = c(0, 30), y_lim = c(-20, 12),
                        dot_clr = cbPalette[4], l_pos=1.1, u_pos = 1.5)

n4_mean <- BlandAltman(main_df = N4_spir_cf_vdps, column1 = "Reader",
                       column2 = "Mean", us1 = "Reader", us2 = "N4",
                       x_lim = c(0, 30), y_lim = c(-20, 12),
                       dot_clr = cbPalette[5], l_pos=1.1, u_pos = 1.5)

n4_thrsh <- BlandAltman(main_df = N4_spir_cf_vdps, column1 = "Reader",
                        column2 = "Thresholding", us1 = "Reader", us2 = "N4",
                        x_lim = c(0, 30), y_lim = c(-20, 12),
                        dot_clr = cbPalette[6], l_pos=1.1, u_pos = 1.5)

n4_mdn <- BlandAltman(main_df = N4_spir_cf_vdps, column1 = "Reader",
                      column2 = "Median", us1 = "Reader", us2 = "N4",
                      x_lim = c(0, 30), y_lim = c(-20, 12),
                      dot_clr = cbPalette[7], l_pos=1.1, u_pos = 1.5)

# #Save the plot as a png file in the specified directory
ggsave("./zR_plots_4ppr/FA_baltman_hrar_rdr_vdp.png", plot = fa_hrar[[1]], width = 5.3, height = 3.9, dpi = 300)
ggsave("./zR_plots_4ppr/FA_baltman_hrar_rdr_vdp_lbl.png", plot = fa_hrar[[2]], width = 5.3, height = 3.9, dpi = 300)
ggsave("./zR_plots_4ppr/FA_baltman_adpt_rdr_vdp.png", plot = fa_adpt[[1]], width = 5.3, height = 3.9, dpi = 300)
ggsave("./zR_plots_4ppr/FA_baltman_adpt_rdr_vdp_lbl.png", plot = fa_adpt[[2]], width = 5.3, height = 3.9, dpi = 300)
ggsave("./zR_plots_4ppr/FA_baltman_prcnt_rdr_vdp.png", plot = fa_prcnt[[1]], width = 5.3, height = 3.9, dpi = 300)
ggsave("./zR_plots_4ppr/FA_baltman_prcnt_rdr_vdp_lbl.png", plot = fa_prcnt[[2]], width = 5.3, height = 3.9, dpi = 300)
ggsave("./zR_plots_4ppr/FA_baltman_mean_rdr_vdp.png", plot = fa_mean[[1]], width = 5.3, height = 3.9, dpi = 300)
ggsave("./zR_plots_4ppr/FA_baltman_mean_rdr_vdp_lbl.png", plot = fa_mean[[2]], width = 5.3, height = 3.9, dpi = 300)
ggsave("./zR_plots_4ppr/FA_baltman_thrsh_rdr_vdp.png", plot = fa_thrsh[[1]], width = 5.3, height = 3.9, dpi = 300)
ggsave("./zR_plots_4ppr/FA_baltman_thrsh_rdr_vdp_lbl.png", plot = fa_thrsh[[2]], width = 5.3, height = 3.9, dpi = 300)
ggsave("./zR_plots_4ppr/FA_baltman_mdn_rdr_vdp.png", plot = fa_mdn[[1]], width = 5.3, height = 3.9, dpi = 300)
ggsave("./zR_plots_4ppr/FA_baltman_mdn_rdr_vdp_lbl.png", plot = fa_mdn[[2]], width = 5.3, height = 3.9, dpi = 300)

ggsave("./zR_plots_4ppr/N4_baltman_hrar_rdr_vdp.png", plot = n4_hrar[[1]], width = 5.3, height = 3.9, dpi = 300)
ggsave("./zR_plots_4ppr/N4_baltman_hrar_rdr_vdp_lbl.png", plot = n4_hrar[[2]], width = 5.3, height = 3.9, dpi = 300)
ggsave("./zR_plots_4ppr/N4_baltman_adpt_rdr_vdp.png", plot = n4_adpt[[1]], width = 5.3, height = 3.9, dpi = 300)
ggsave("./zR_plots_4ppr/N4_baltman_adpt_rdr_vdp_lbl.png", plot = n4_adpt[[2]], width = 5.3, height = 3.9, dpi = 300)
ggsave("./zR_plots_4ppr/N4_baltman_prcnt_rdr_vdp.png", plot = n4_prcnt[[1]], width = 5.3, height = 3.9, dpi = 300)
ggsave("./zR_plots_4ppr/N4_baltman_prcnt_rdr_vdp_lbl.png", plot = n4_prcnt[[2]], width = 5.3, height = 3.9, dpi = 300)
ggsave("./zR_plots_4ppr/N4_baltman_mean_rdr_vdp.png", plot = n4_mean[[1]], width = 5.3, height = 3.9, dpi = 300)
ggsave("./zR_plots_4ppr/N4_baltman_mean_rdr_vdp_lbl.png", plot = n4_mean[[2]], width = 5.3, height = 3.9, dpi = 300)
ggsave("./zR_plots_4ppr/N4_baltman_thrsh_rdr_vdp.png", plot = n4_thrsh[[1]], width = 5.3, height = 3.9, dpi = 300)
ggsave("./zR_plots_4ppr/N4_baltman_thrsh_rdr_vdp_lbl.png", plot = n4_thrsh[[2]], width = 5.3, height = 3.9, dpi = 300)
ggsave("./zR_plots_4ppr/N4_baltman_mdn_rdr_vdp.png", plot = n4_mdn[[1]], width = 5.3, height = 3.9, dpi = 300)
ggsave("./zR_plots_4ppr/N4_baltman_mdn_rdr_vdp_lbl.png", plot = n4_mdn[[2]], width = 5.3, height = 3.9, dpi = 300)




