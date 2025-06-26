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

linreg_plt <- function(df_in, x_in, y_in, x_label, y_label, x_lim = c(0, 40),
                       y_lim = c(50, 150), corrltn_method = "pearson") {
  # Check if correlation method is valid
  if (!corrltn_method %in% c("pearson", "spearman")) {
    stop("Only pearson or spearman correlations are accepted.")
  }
  
  # Create formula for linear regression
  yx_formula <- as.formula(paste(y_in, "~", x_in))
  fitt <- lm(yx_formula, data = df_in)
  slope <- coef(fitt)[2]
  intercept <- coef(fitt)[1]
  print(fitt)
  
  # Base scatter plot with regression line
  linreg <- ggpubr::ggscatter(df_in, x = x_in, y = y_in,
                              add = "reg.line",  # Add regression line
                              xlim = x_lim, ylim = y_lim,
                              add.params = list(color = "#000000", size = 2, fill = "red", alpha = 0.15),
                              conf.int = TRUE, fullrange = TRUE,
                              cor.method = corrltn_method,
                              color = "#000000", shape = 19, size = 4, # Points color, shape and size
                              ggtheme = theme_bw(),
                              font.x = c(22, "bold", "#000000"),
                              font.y = c(22, "bold", "#000000"),
                              font.tickslab = c(22, "bold", "#000000")) +
    xlab(x_label) + ylab(y_label) +
    geom_abline(intercept = intercept, slope = slope,
                linetype = "dotted", linewidth = 2, color = "#000000") +
    theme(legend.position = "none")
  
  # Calculate correlation coefficient (R) and p-value
  cor_test <- cor.test(df_in[[x_in]], df_in[[y_in]], method = corrltn_method)
  r_value <- cor_test$estimate
  p_value <- cor_test$p.value
  
  # Format R and p-value labels
  r_label <- ifelse(corrltn_method == "pearson", "R", "rho")
  r_text <- paste0(r_label, " == ", round(r_value, 2))
  
  # Format p-value correctly for parsing
  if (p_value < 0.001) {
    p_text <- "P < 0.001"
  } else {
    p_text <- paste0("P == ", round(p_value, 3))
  }
  
  # Linear regression equation
  fit_eq <- parse(text = sprintf("y == %.2f * x + %.2f", slope, intercept))
  
  # Add correlation coefficient (R), p-value, and equation to the plot
  linreg_txt <- linreg + 
    annotate("text", x = Inf, y = Inf, 
             label = paste(r_text, p_text, sep = "~`,`~"), 
             hjust = 1.1, vjust = 2, color = "#000000", size = 6, parse = TRUE) +
    annotate("text", x = Inf, y = Inf, label = fit_eq, 
             hjust = 1.1, vjust = 3.5, color = "#000000", size = 6, parse = TRUE)
  
  print(linreg_txt)
  
  # Return handles to both plots
  return(list(linreg, linreg_txt))
}
################################################################################
# ##Load the data from both CSV files and arrange it to plot
N4_spir_cf_vdps <- read.csv("./IRC740H_visit1_spiral_cf/vdp_analysis_results_May2025/N4_combined_vdps.csv")
FA_spir_cf_vdps <- read.csv("./IRC740H_visit1_spiral_cf/vdp_analysis_results_May2025/FA_combined_vdps.csv")

FEV1_rdr <- linreg_plt(N4_spir_cf_vdps, "FEV1", "Reader",
                      expression(bold(FEV1*" "*"(%)")), expression(bold(VDP[Reader]*" "*"(%)")),
                      x_lim=c(70, 150), y_lim=c(-10, 35), corrltn_method="pearson")

N4_FEV1_hrar <- linreg_plt(N4_spir_cf_vdps, "FEV1", "Hierarchical",
                        expression(bold(FEV1*" "*"(%)")), expression(bold(Hierarchical*" "*VDP[N4]*" "*"(%)")),
                        x_lim=c(70, 150), y_lim=c(-10, 35), corrltn_method="pearson")
N4_FEV1_adpt <- linreg_plt(N4_spir_cf_vdps, "FEV1", "Adaptive",
                          expression(bold(FEV1*" "*"(%)")), expression(bold(Adaptive*" "*VDP[N4]*" "*"(%)")),
                          x_lim=c(70, 150), y_lim=c(-10, 35), corrltn_method="pearson")
N4_FEV1_prctl <- linreg_plt(N4_spir_cf_vdps, "FEV1", "Percentile",
                          expression(bold(FEV1*" "*"(%)")), expression(bold(Percentile*" "*VDP[N4]*" "*"(%)")),
                          x_lim=c(70, 150), y_lim=c(-10, 35), corrltn_method="pearson")
N4_FEV1_mnlb <- linreg_plt(N4_spir_cf_vdps, "FEV1", "Mean",
                           expression(bold(FEV1*" "*"(%)")), expression(bold(Mean[LB]*" "*VDP[N4]*" "*"(%)")),
                           x_lim=c(70, 150), y_lim=c(-10, 35), corrltn_method="pearson")
N4_FEV1_thrsh <- linreg_plt(N4_spir_cf_vdps, "FEV1", "Thresholding",
                          expression(bold(FEV1*" "*"(%)")), expression(bold(Thresholding*" "*VDP[N4]*" "*"(%)")),
                          x_lim=c(70, 150), y_lim=c(-10, 35), corrltn_method="pearson")
N4_FEV1_mnglb <- linreg_plt(N4_spir_cf_vdps, "FEV1", "Median",
                          expression(bold(FEV1*" "*"(%)")), expression(bold(Mean[GLB]*" "*VDP[N4]*" "*"(%)")),
                          x_lim=c(70, 150), y_lim=c(-10, 35), corrltn_method="pearson")


FA_FEV1_hrar <- linreg_plt(FA_spir_cf_vdps, "FEV1", "Hierarchical",
                          expression(bold(FEV1*" "*"(%)")), expression(bold(Hierarchical*" "*VDP[FA]*" "*"(%)")),
                          x_lim=c(70, 150), y_lim=c(-10, 35), corrltn_method="pearson")
FA_FEV1_adpt <- linreg_plt(FA_spir_cf_vdps, "FEV1", "Adaptive",
                          expression(bold(FEV1*" "*"(%)")), expression(bold(Adaptive*" "*VDP[FA]*" "*"(%)")),
                          x_lim=c(70, 150), y_lim=c(-10, 35), corrltn_method="pearson")
FA_FEV1_prctl <- linreg_plt(FA_spir_cf_vdps, "FEV1", "Percentile",
                           expression(bold(FEV1*" "*"(%)")), expression(bold(Percentile*" "*VDP[FA]*" "*"(%)")),
                           x_lim=c(70, 150), y_lim=c(-10, 35), corrltn_method="pearson")
FA_FEV1_mnlb <- linreg_plt(FA_spir_cf_vdps, "FEV1", "Mean",
                          expression(bold(FEV1*" "*"(%)")), expression(bold(Mean[LB]*" "*VDP[FA]*" "*"(%)")),
                          x_lim=c(70, 150), y_lim=c(-10, 35), corrltn_method="pearson")
FA_FEV1_thrsh <- linreg_plt(FA_spir_cf_vdps, "FEV1", "Thresholding",
                           expression(bold(FEV1*" "*"(%)")), expression(bold(Thresholding*" "*VDP[FA]*" "*"(%)")),
                           x_lim=c(70, 150), y_lim=c(-10, 35), corrltn_method="pearson")
FA_FEV1_mnglb <- linreg_plt(FA_spir_cf_vdps, "FEV1", "Median",
                           expression(bold(FEV1*" "*"(%)")), expression(bold(Mean[GLB]*" "*VDP[FA]*" "*"(%)")),
                           x_lim=c(70, 150), y_lim=c(-10, 35), corrltn_method="pearson")


# #Save the plot as a png file in the specified directory
ggsave("./zR_plots_4ppr/cf_FEV1_rdr_vdp_linreg_plain.png", plot = FEV1_rdr[[1]], width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/cf_FEV1_rdr_vdp_linreg_p.png", plot = FEV1_rdr[[2]], width = 4.5, height = 3.7, dpi = 300)


ggsave("./zR_plots_4ppr/cf_FEV1_N4_hrar_vdp_linreg_plain.png", plot = N4_FEV1_hrar[[1]], width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/cf_FEV1_N4_hrar_vdp_linreg_p.png", plot = N4_FEV1_hrar[[2]], width = 4.5, height = 3.7, dpi = 300)

ggsave("./zR_plots_4ppr/cf_FEV1_N4_adpt_vdp_linreg_plain.png", plot = N4_FEV1_adpt[[1]], width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/cf_FEV1_N4_adpt_vdp_linreg_p.png", plot = N4_FEV1_adpt[[2]], width = 4.5, height = 3.7, dpi = 300)

ggsave("./zR_plots_4ppr/cf_FEV1_N4_prctl_vdp_linreg_plain.png", plot = N4_FEV1_prctl[[1]], width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/cf_FEV1_N4_prctl_vdp_linreg_p.png", plot = N4_FEV1_prctl[[2]], width = 4.5, height = 3.7, dpi = 300)

ggsave("./zR_plots_4ppr/cf_FEV1_N4_mnlb_vdp_linreg_plain.png", plot = N4_FEV1_mnlb[[1]], width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/cf_FEV1_N4_mnlb_vdp_linreg_p.png", plot = N4_FEV1_mnlb[[2]], width = 4.5, height = 3.7, dpi = 300)

ggsave("./zR_plots_4ppr/cf_FEV1_N4_thrsh_vdp_linreg_plain.png", plot = N4_FEV1_thrsh[[1]], width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/cf_FEV1_N4_thrsh_vdp_linreg_p.png", plot = N4_FEV1_thrsh[[2]], width = 4.5, height = 3.7, dpi = 300)

ggsave("./zR_plots_4ppr/cf_FEV1_N4_mnglb_vdp_linreg_plain.png", plot = N4_FEV1_mnglb[[1]], width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/cf_FEV1_N4_mnglb_vdp_linreg_p.png", plot = N4_FEV1_mnglb[[2]], width = 4.5, height = 3.7, dpi = 300)


ggsave("./zR_plots_4ppr/cf_FEV1_FA_hrar_vdp_linreg_plain.png", plot = FA_FEV1_hrar[[1]], width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/cf_FEV1_FA_hrar_vdp_linreg_p.png", plot = FA_FEV1_hrar[[2]], width = 4.5, height = 3.7, dpi = 300)

ggsave("./zR_plots_4ppr/cf_FEV1_FA_adpt_vdp_linreg_plain.png", plot = FA_FEV1_adpt[[1]], width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/cf_FEV1_FA_adpt_vdp_linreg_p.png", plot = FA_FEV1_adpt[[2]], width = 4.5, height = 3.7, dpi = 300)

ggsave("./zR_plots_4ppr/cf_FEV1_FA_prctl_vdp_linreg_plain.png", plot = FA_FEV1_prctl[[1]], width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/cf_FEV1_FA_prctl_vdp_linreg_p.png", plot = FA_FEV1_prctl[[2]], width = 4.5, height = 3.7, dpi = 300)

ggsave("./zR_plots_4ppr/cf_FEV1_FA_mnlb_vdp_linreg_plain.png", plot = FA_FEV1_mnlb[[1]], width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/cf_FEV1_FA_mnlb_vdp_linreg_p.png", plot = FA_FEV1_mnlb[[2]], width = 4.5, height = 3.7, dpi = 300)

ggsave("./zR_plots_4ppr/cf_FEV1_FA_thrsh_vdp_linreg_plain.png", plot = FA_FEV1_thrsh[[1]], width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/cf_FEV1_FA_thrsh_vdp_linreg_p.png", plot = FA_FEV1_thrsh[[2]], width = 4.5, height = 3.7, dpi = 300)

ggsave("./zR_plots_4ppr/cf_FEV1_FA_mnglb_vdp_linreg_plain.png", plot = FA_FEV1_mnglb[[1]], width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/cf_FEV1_FA_mnglb_vdp_linreg_p.png", plot = FA_FEV1_mnglb[[2]], width = 4.5, height = 3.7, dpi = 300)

