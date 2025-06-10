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

##Colorblind friendly palette:
cbPalette <- c("#888888", "#CC6677", "#882255", "#332288", "#6699CC", "#DDCC77",
               "#999933", "#AA4499", "#661100", "#44AA99", "#117733", "#88CCEE")
##Define order of the data -- if needed
# plt_order <- c("Reader", "Adaptive","Hierarchical", "Mean", "Percentile",
#                "Thresholding", "Median")

boxplt <- function(df_in, x_in, y_in, ylabel=NULL, xlabel="", ylimit=NULL,
                   bxfill = "white", bxpallet=NULL, pt_clrs=c("#000000", "#009e73"),
                   stat_tst="wilcox", paired=FALSE, py_pos=NULL, addp_eq=FALSE){
  require(ggpubr)
  # Check if ylabel, ylimit and py_pos are provided
  if (is.null(ylabel) || is.null(ylimit) || is.null(py_pos)){
    ylabel <- y_in
    ylimit <- c(0, 1.5 * max(df_in[[y_in]], na.rm = TRUE))
    py_pos <- 1.35 * max(df_in[[y_in]], na.rm = TRUE) }
  ##Paste x, y inputs for statistical tests
  yx_formula <- as.formula(paste(y_in, "~", x_in))
  mean_val <- aggregate(yx_formula, df_in, mean)
  print("Mean")
  print(mean_val)
  sd_val <- aggregate(yx_formula, df_in, sd)
  print("Standard-deviation")
  print(sd_val)
  range_val <- aggregate(yx_formula, df_in, function(x) {
    paste0("[", min(x, na.rm = TRUE), ", ", max(x, na.rm = TRUE), "]")
  })
  print("Range")
  print(range_val)
  
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
  ##Add p-value on the plot
  pval <- pval %>% add_xy_position(x = x_in)
  if (addp_eq == TRUE) {plabel = "P={scales::pvalue(p, accuracy = 0.001)}"}
  else {plabel = "P{scales::pvalue(p, accuracy = 0.001)}"}
  bxp_p <-  bxp_ + stat_pvalue_manual(pval, label = plabel, #"P = {p.adj}"
                                      y.position=py_pos, label.size = 8,
                                      bracket.size = 0.8, tip.length = 0.025,
                                      vjust=-0.25)
  print(bxp_p)
  return(list(bxp_, bxp_p))
}


################################################################################
# ##Load the data from CSV files and arrange for plotting
##CF
FA_spir_cf <- read.csv("./IRC740H_visit1_spiral_cf/vdp_analysis_results_May2025/FA_combined_vdps.csv")
N4_spir_cf <- read.csv("./IRC740H_visit1_spiral_cf/vdp_analysis_results_May2025/N4_combined_vdps.csv")
##Healthy
FA_spir_ctrl <- read.csv("./age_matched_spiral_healthy/vdp_analysis_results_May2025/FA_combined_vdps.csv")
N4_spir_ctrl <- read.csv("./age_matched_spiral_healthy/vdp_analysis_results_May2025/N4_combined_vdps.csv")
# #Create a combined data frame with an indicator for Correction and Category
combined_fa <- rbind(transform(FA_spir_ctrl, Category="Healthy"),
                            transform(FA_spir_cf, Category="CF"))
combined_n4 <- rbind(transform(N4_spir_ctrl, Category="Healthy"),
                            transform(N4_spir_cf, Category="CF")) # Correction = "N4",

# (df_in, x_in, y_in, ylabel=NULL, xlabel="", ylimit=NULL,
#   bxfill = "white", bxpallet=NULL, pt_clrs=c("#000000", "#009e73"),
#   stat_tst="wilcox", paired=FALSE, py_pos=NULL, addp_eq=FALSE)

##############################Reader VDP
rdr_vdp <- boxplt(combined_fa, "Category", "Reader", ylabel=expression(bold(Reader* " " *VDP* " "* "(%)")),
                  xlabel="", ylimit=c(0, 40), bxfill= cbPalette[1], py_pos=36, addp_eq=FALSE)
##############################Box plot: FA VDP (CF vs Healthy)
fa_hr <- boxplt(combined_fa, "Category", "Hierarchical", ylabel=expression(bold(VDP[FA]* " "* "(%)")),
                xlabel="", ylimit=c(0, 35), bxfill= cbPalette[3], py_pos=31, addp_eq=TRUE)
fa_adp <- boxplt(combined_fa, "Category", "Adaptive", ylabel=expression(bold(VDP[FA]*" "* "(%)")),
                 xlabel="", ylimit=c(0, 55), bxfill= cbPalette[2], py_pos=49, addp_eq=TRUE)
fa_pct <- boxplt(combined_fa, "Category", "Percentile", ylabel=expression(bold(VDP[FA]* " "* "(%)")),
                xlabel="", ylimit=c(0, 35), bxfill= cbPalette[4], py_pos=31, addp_eq=TRUE)
fa_mn <- boxplt(combined_fa, "Category", "Mean", ylabel=expression(bold(VDP[FA]* " "* "(%)")),
                 xlabel="", ylimit=c(0, 35), bxfill= cbPalette[5], py_pos=31, addp_eq=TRUE)
fa_th <- boxplt(combined_fa, "Category", "Thresholding", ylabel=expression(bold(VDP[FA]* " "* "(%)")),
                xlabel="", ylimit=c(0, 35), bxfill= cbPalette[6], py_pos=31, addp_eq=TRUE)
fa_md <- boxplt(combined_fa, "Category", "Median", ylabel=expression(bold(VDP[FA]* " "* "(%)")),
                 xlabel="", ylimit=c(0, 35), bxfill= cbPalette[7], py_pos=31, addp_eq=TRUE)



# #Save the plot as a png file in the specified directory
ggsave("./zR_plots_4ppr/Healthy_vs_CF_spir_rdr_vdp.png", plot = rdr_vdp[[1]], width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/Healthy_vs_CF_spir_rdr_vdp_p.png", plot = rdr_vdp[[2]], width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/Healthy_vs_CF_spir_FA_adp_vdp.png", plot = fa_adp[[1]], width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/Healthy_vs_CF_spir_FA_adp_vdp_p.png", plot = fa_adp[[2]], width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/Healthy_vs_CF_spir_FA_hr_vdp.png", plot = fa_hr[[1]], width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/Healthy_vs_CF_spir_FA_hr_vdp_p.png", plot = fa_hr[[2]], width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/Healthy_vs_CF_spir_FA_thr_vdp.png", plot = fa_th[[1]], width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/Healthy_vs_CF_spir_FA_thr_vdp_p.png", plot = fa_th[[2]], width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/Healthy_vs_CF_spir_FA_pct_vdp.png", plot = fa_pct[[1]], width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/Healthy_vs_CF_spir_FA_pct_vdp_p.png", plot = fa_pct[[2]], width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/Healthy_vs_CF_spir_FA_mn_vdp.png", plot = fa_mn[[1]], width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/Healthy_vs_CF_spir_FA_mn_vdp_p.png", plot = fa_mn[[2]], width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/Healthy_vs_CF_spir_FA_mdn_vdp.png", plot = fa_md[[1]], width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/Healthy_vs_CF_spir_FA_mdn_vdp_p.png", plot = fa_md[[2]], width = 4.5, height = 3.7, dpi = 300)



##############################Box plot: N4 VDP (CF vs Healthy)
n4_hr <- boxplt(combined_n4, "Category", "Hierarchical", ylabel=expression(bold(VDP[N4]* " "* "(%)")),
                xlabel="", ylimit=c(0, 40), bxfill= cbPalette[3], py_pos=34, addp_eq=TRUE)
n4_adp <- boxplt(combined_n4, "Category", "Adaptive", ylabel=expression(bold(VDP[N4]* " "* "(%)")),
                 xlabel="", ylimit=c(0, 40), bxfill= cbPalette[2], py_pos=34, addp_eq=TRUE)
n4_pct <- boxplt(combined_n4, "Category", "Percentile", ylabel=expression(bold(VDP[N4]* " "* "(%)")),
                 xlabel="", ylimit=c(0, 40), bxfill= cbPalette[4], py_pos=34, addp_eq=FALSE)
n4_mn <- boxplt(combined_n4, "Category", "Mean", ylabel=expression(bold(VDP[N4]* " "* "(%)")),
                 xlabel="", ylimit=c(0, 40), bxfill= cbPalette[5], py_pos=34, addp_eq=FALSE)
n4_th <- boxplt(combined_n4, "Category", "Thresholding", ylabel=expression(bold(VDP[N4]* " "* "(%)")),
                xlabel="", ylimit=c(0, 40), bxfill= cbPalette[6], py_pos=34, addp_eq=FALSE)
n4_md <- boxplt(combined_n4, "Category", "Median", ylabel=expression(bold(VDP[N4]* " "* "(%)")),
                xlabel="", ylimit=c(0, 40), bxfill= cbPalette[7], py_pos=34, addp_eq=FALSE)
# #Save the plot as a png file in the specified directory

ggsave("./zR_plots_4ppr/Healthy_vs_CF_spir_N4_adp_vdp.png", plot = n4_adp[[1]], width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/Healthy_vs_CF_spir_N4_adp_vdp_p.png", plot = n4_adp[[2]], width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/Healthy_vs_CF_spir_N4_hr_vdp.png", plot = n4_hr[[1]], width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/Healthy_vs_CF_spir_N4_hr_vdp_p.png", plot = n4_hr[[2]], width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/Healthy_vs_CF_spir_N4_thr_vdp.png", plot = n4_th[[1]], width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/Healthy_vs_CF_spir_N4_thr_vdp_p.png", plot = n4_th[[2]], width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/Healthy_vs_CF_spir_N4_pct_vdp.png", plot = n4_pct[[1]], width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/Healthy_vs_CF_spir_N4_pct_vdp_p.png", plot = n4_pct[[2]], width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/Healthy_vs_CF_spir_N4_mn_vdp.png", plot = n4_mn[[1]], width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/Healthy_vs_CF_spir_N4_mn_vdp_p.png", plot = n4_mn[[2]], width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/Healthy_vs_CF_spir_N4_mdn_vdp.png", plot = n4_md[[1]], width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/Healthy_vs_CF_spir_N4_mdn_vdp_p.png", plot = n4_md[[2]], width = 4.5, height = 3.7, dpi = 300)


# ##If want to add mean line on box plot
# stat_summary(fun=mean, geom="crossbar", aes(fill = x_in), linewidth=0.25, linetype="dotted",
#              width=0.6, color = "black", position=position_dodge(0.75))

