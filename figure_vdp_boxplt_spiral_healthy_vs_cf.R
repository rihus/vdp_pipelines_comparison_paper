#load libraries
library(ggplot2)
library(ggpubr)
library(blandr)
library(dplyr)
library(tidyr)
library(rstatix)
library(PMCMRplus)

# Set the working directory to the path where your CSV files are located
setwd("C:/Users/HUSDQ4/OneDrive - cchmc/cincy_work/all_projects_data_work/vdp_analysis/analysis_comparisons")

boxplt <- function(df_in, x_in, y_in, ylabel=NULL, xlabel="", ylimit=NULL,
                   bxfill = "white", bxpallet=NULL, pt_clrs=c("#000000", "#009e73"),
                   stat_tst="wilcox", paired=FALSE, py_pos=NULL){
  require(ggpubr)
  # Check if ylabel, ylimit and py_pos are provided
  if (is.null(ylabel) || is.null(ylimit) || is.null(py_pos)){
    ylabel <- y_in
    ylimit <- c(0, 1.5 * max(df_in[[y_in]], na.rm = TRUE))
    py_pos <- 1.35 * max(df_in[[y_in]], na.rm = TRUE) }
  ##Paste x, y inputs for statistical tests
  yx_formula <- as.formula(paste(y_in, "~", x_in))
  mean_val <- aggregate(yx_formula, df_in, mean)
  print(mean_val)
  sd_val <- aggregate(yx_formula, df_in, sd)
  print(sd_val)
  ## Statistical test
  pval <- NULL
  if (stat_tst == "wilcox") {
    pval <- df_in %>%
      wilcox_test(yx_formula, paired = paired) %>%
      adjust_pvalue(method = 'bonferroni') %>%
      add_significance()
    print(pval)
  } else if (stat_tst == "ttest") {
    pval <- df_in %>%
      t_test(yx_formula, paired = paired) %>%
      adjust_pvalue(method = 'bonferroni') %>%
      add_significance()
    print(pval)
  } else {
    stop("For stat_tst, only Wilcoxon (wilcox) or T-test (ttest) are accepted")
  }
  ##Print adjusted p-value if test is paired, otherwise p
  # if (paired) {plabel= "P = {p.adj}"} else {plabel= "P = {p}"}
  ##Unconnected box-plot 
  bxp_ <- ggboxplot(df_in, x = x_in, y = y_in,
                    ylim = ylimit, fill = bxfill, pallete = bxpallet,
                    outlier.shape = NA,
                    font.x = c(22, "bold", "#000000"),
                    font.y = c(24, "bold", "#000000"), 
                    font.tickslab = c(22, "bold", "#000000")) +
    xlab(xlabel) +
    ylab(ylabel) +
    geom_jitter(aes(color = as.factor(df_in[[x_in]])), width=0.2, size=2, alpha=1) +
    scale_color_manual(values = pt_clrs) +
    theme(legend.position = "none")
  bxp_ <- bxp_ + geom_vline(xintercept = Inf, linetype = "solid")
  bxp_ <- bxp_ + geom_hline(yintercept = Inf, linetype = "solid")
  print(bxp_)
  pval <- pval %>% add_xy_position(x = x_in)
  bxp_p <-  bxp_ + stat_pvalue_manual(pval, label = "P = {p.adj}", ##"P={scales::pvalue(p)}"
                                      y.position=py_pos, label.size = 8,
                                      bracket.size = 0.8, tip.length = 0.025,
                                      vjust=-0.25)
  print(bxp_p)
  return(list(bxp_, bxp_p))
}


################################################################################
# ##Load the data from CSV files and arrange for plotting
##CF
FA_spir_cf <- read.csv("./IRC740H_visit1_spiral_cf/vdp_analysis_results_October2024/FA_combined_vdps.csv")
N4_spir_cf <- read.csv("./IRC740H_visit1_spiral_cf/vdp_analysis_results_October2024/N4_combined_vdps.csv")
##Healthy
FA_spir_ctrl <- read.csv("./age_matched_spiral_healthy/vdp_analysis_results_October2024/FA_combined_vdps.csv")
N4_spir_ctrl <- read.csv("./age_matched_spiral_healthy/vdp_analysis_results_October2024/N4_combined_vdps.csv")
# #Create a combined data frame with an indicator for Correction and Category
combined_fa <- rbind(transform(FA_spir_ctrl, Category="Healthy"),
                            transform(FA_spir_cf, Category="CF"))
combined_n4 <- rbind(transform(N4_spir_ctrl, Category="Healthy"),
                            transform(N4_spir_cf, Category="CF")) # Correction = "N4",

##############################Reader VDP
rdr_vdp <- boxplt(combined_fa, "Category", "Reader", ylabel=expression(bold(VDP[Reader]* " "* "(%)")),
                  xlabel="", ylimit=c(0, 35), py_pos=31)
##############################Box plot: FA VDP (CF vs Healthy)
fa_th <- boxplt(combined_fa, "Category", "Thresholding", ylabel=expression(bold(VDP[Thresholding]* " "* "(%)")),
                  xlabel="", ylimit=c(0, 35), py_pos=31)
fa_pct <- boxplt(combined_fa, "Category", "Percentile", ylabel=expression(bold(VDP[Percentile]* " "* "(%)")),
                xlabel="", ylimit=c(0, 35), py_pos=31)
fa_md <- boxplt(combined_fa, "Category", "Median", ylabel=expression(bold(VDP[Median]* " "* "(%)")),
                 xlabel="", ylimit=c(0, 35), py_pos=31)
fa_hr <- boxplt(combined_fa, "Category", "Hierarchical", ylabel=expression(bold(VDP[Hierarchical]* " "* "(%)")),
                xlabel="", ylimit=c(0, 35), py_pos=31)
fa_adp <- boxplt(combined_fa, "Category", "Adaptive", ylabel=expression(bold(VDP[Adaptive]* " "* "(%)")),
                xlabel="", ylimit=c(0, 55), py_pos=31)

# #Save the plot as a png file in the specified directory
ggsave("./zR_plots_4abs/spir_vdp_fa_rdr_Healthy_vs_CF_nop.png", plot = rdr_vdp[[1]], width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4abs/spir_vdp_fa_rdr_Healthy_vs_CF_p.png", plot = rdr_vdp[[2]], width = 4.5, height = 3.7, dpi = 300)

ggsave("./zR_plots_4abs/spir_vdp_fa_th_Healthy_vs_CF_nop.png", plot = fa_th[[1]], width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4abs/spir_vdp_fa_th_Healthy_vs_CF_p.png", plot = fa_th[[2]], width = 4.5, height = 3.7, dpi = 300)

ggsave("./zR_plots_4abs/spir_vdp_fa_pct_Healthy_vs_CF_nop.png", plot = fa_pct[[1]], width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4abs/spir_vdp_fa_pct_Healthy_vs_CF_p.png", plot = fa_pct[[2]], width = 4.5, height = 3.7, dpi = 300)

ggsave("./zR_plots_4abs/spir_vdp_fa_md_Healthy_vs_CF_nop.png", plot = fa_md[[1]], width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4abs/spir_vdp_fa_md_Healthy_vs_CF_p.png", plot = fa_md[[2]], width = 4.5, height = 3.7, dpi = 300)

ggsave("./zR_plots_4abs/spir_vdp_fa_hr_Healthy_vs_CF_nop.png", plot = fa_hr[[1]], width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4abs/spir_vdp_fa_hr_Healthy_vs_CF_p.png", plot = fa_hr[[2]], width = 4.5, height = 3.7, dpi = 300)

ggsave("./zR_plots_4abs/spir_vdp_fa_adp_Healthy_vs_CF_nop.png", plot = fa_adp[[1]], width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4abs/spir_vdp_fa_adp_Healthy_vs_CF_p.png", plot = fa_adp[[2]], width = 4.5, height = 3.7, dpi = 300)

##############################Box plot: N4 VDP (CF vs Healthy)
n4_th <- boxplt(combined_n4, "Category", "Thresholding", ylabel=expression(bold(VDP[Thresholding]* " "* "(%)")),
                xlabel="", ylimit=c(0, 40), py_pos=36)
n4_pct <- boxplt(combined_n4, "Category", "Percentile", ylabel=expression(bold(VDP[Percentile]* " "* "(%)")),
                 xlabel="", ylimit=c(0, 40), py_pos=36)
n4_md <- boxplt(combined_n4, "Category", "Median", ylabel=expression(bold(VDP[Median]* " "* "(%)")),
                xlabel="", ylimit=c(0, 40), py_pos=36)
n4_hr <- boxplt(combined_n4, "Category", "Hierarchical", ylabel=expression(bold(VDP[Hierarchical]* " "* "(%)")),
                xlabel="", ylimit=c(0, 40), py_pos=36)
n4_adp <- boxplt(combined_n4, "Category", "Adaptive", ylabel=expression(bold(VDP[Adaptive]* " "* "(%)")),
                 xlabel="", ylimit=c(0, 40), py_pos=36)

# #Save the plot as a png file in the specified directory
ggsave("./zR_plots_4abs/spir_vdp_n4_th_Healthy_vs_CF_nop.png", plot = n4_th[[1]], width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4abs/spir_vdp_n4_th_Healthy_vs_CF_p.png", plot = n4_th[[2]], width = 4.5, height = 3.7, dpi = 300)

ggsave("./zR_plots_4abs/spir_vdp_n4_pct_Healthy_vs_CF_nop.png", plot = n4_pct[[1]], width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4abs/spir_vdp_n4_pct_Healthy_vs_CF_p.png", plot = n4_pct[[2]], width = 4.5, height = 3.7, dpi = 300)

ggsave("./zR_plots_4abs/spir_vdp_n4_md_Healthy_vs_CF_nop.png", plot = n4_md[[1]], width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4abs/spir_vdp_n4_md_Healthy_vs_CF_p.png", plot = n4_md[[2]], width = 4.5, height = 3.7, dpi = 300)

ggsave("./zR_plots_4abs/spir_vdp_n4_hr_Healthy_vs_CF_nop.png", plot = n4_hr[[1]], width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4abs/spir_vdp_n4_hr_Healthy_vs_CF_p.png", plot = n4_hr[[2]], width = 4.5, height = 3.7, dpi = 300)

ggsave("./zR_plots_4abs/spir_vdp_n4_adp_Healthy_vs_CF_nop.png", plot = n4_adp[[1]], width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4abs/spir_vdp_n4_adp_Healthy_vs_CF_p.png", plot = n4_adp[[2]], width = 4.5, height = 3.7, dpi = 300)


# ##If want to add mean line on box plot
# stat_summary(fun=mean, geom="crossbar", aes(fill = x_in), linewidth=0.25, linetype="dotted",
#              width=0.6, color = "black", position=position_dodge(0.75))

