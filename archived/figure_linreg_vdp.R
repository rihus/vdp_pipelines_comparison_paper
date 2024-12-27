#load ggplot2
library(ggplot2)
library(ggpubr)
library(blandr)
library(dplyr)
library(tidyr)
library(rstatix)

# #Set the working directory to the path where your CSV files are located
setwd("C:/Users/HUSDQ4/OneDrive - cchmc/cincy_work/all_projects_data_work/vdp_analysis/CFNonCF_Bronch")

################################################################################
# ##Load the data from both CSV files and arrange it to plot
N4_corr_spir_cf <- read.csv("./IRC740H_2Dspiral_CF/snr_vdp_hvp_analysis_results_July2024/N4_corr_median_analysis_results.csv")
FA_corr_spir_cf <- read.csv("./IRC740H_2Dspiral_CF/snr_vdp_hvp_analysis_results_July2024/FA_corr_median_analysis_results.csv")

N4_corr_spir_ctrl <- read.csv("./IRC740H_2Dspiral_noBininng_healthy/snr_vdp_hvp_analysis_results_Apr2024/N4_corr_median_analysis_results.csv")
FA_corr_spir_ctrl <- read.csv("./IRC740H_2Dspiral_noBininng_healthy/snr_vdp_hvp_analysis_results_Apr2024/FA_corr_median_analysis_results.csv")

# #Create a combined data frame with an indicator for Correction and Category
combined_spir_data <- rbind(transform(N4_corr_spir_cf, Correction = "N4", Category = "CF"),
                               transform(FA_corr_spir_cf, Correction = "FA", Category = "CF"),
                               transform(N4_corr_spir_ctrl, Correction = "N4", Category = "Healthy"),
                               transform(FA_corr_spir_ctrl, Correction = "FA", Category = "Healthy"))

# combined_spir_data_CF <- rbind(transform(N4_corr_spir_cf, Correction = "N4", Category = "CF"),
#                            transform(FA_corr_spir_cf, Correction = "FA", Category = "CF"))
# combined_spir_data_H <- rbind(transform(N4_corr_spir_ctrl, Correction = "N4", Category = "Healthy"),
#                                transform(FA_corr_spir_ctrl, Correction = "FA", Category = "Healthy"))

##Create a new data frame selecting shared subjects for connected box plots
common_subjs <- combined_spir_data$Subject_id[duplicated(combined_spir_data$Subject_id) | 
                                                      duplicated(combined_spir_data$Subject_id, fromLast = TRUE)]
common_subjs_spir <- combined_spir_data[combined_spir_data$Subject_id %in% common_subjs, , drop = FALSE]

# common_subjs_CF <- combined_spir_data_CF$Subject_id[duplicated(combined_spir_data_CF$Subject_id) | 
#                                                  duplicated(combined_spir_data_CF$Subject_id, fromLast = TRUE)]
# common_subjs_spir_CF <- combined_spir_data_CF[combined_spir_data_CF$Subject_id %in% common_subjs_CF, , drop = FALSE]
# common_subjs_H <- combined_spir_data_H$Subject_id[duplicated(combined_spir_data_H$Subject_id) | 
#                                                       duplicated(combined_spir_data_H$Subject_id, fromLast = TRUE)]
# common_subjs_spir_H <- combined_spir_data_H[combined_spir_data_H$Subject_id %in% common_subjs_H, , drop = FALSE]

# new_df_cf <- data.frame(N4_vdp = common_subjs_spir_CF$VDP[common_subjs_spir_CF$Correction == "N4"],
#   FA_vdp = common_subjs_spir_CF$VDP[common_subjs_spir_CF$Correction == "FA"],
#   Category = common_subjs_spir_CF$Category[common_subjs_spir_CF$Correction == "N4"])
# new_df_h <- data.frame(N4_vdp = common_subjs_spir_H$VDP[common_subjs_spir_H$Correction == "N4"],
#   FA_vdp = common_subjs_spir_H$VDP[common_subjs_spir_H$Correction == "FA"],
#   Category = common_subjs_spir_H$Category[common_subjs_spir_H$Correction == "N4"])
new_df <- data.frame(N4_vdp = common_subjs_spir$VDP[common_subjs_spir$Correction == "N4"],
  FA_vdp = common_subjs_spir$VDP[common_subjs_spir$Correction == "FA"],
  Category = common_subjs_spir$Category[common_subjs_spir$Correction == "N4"])

###########Bland-Altman plot: VDP
new_df$mean_vdp <- rowMeans(new_df[, c("N4_vdp", "FA_vdp")])
new_df$diff_vdp <- as.numeric(new_df[["N4_vdp"]]) - as.numeric(new_df[["FA_vdp"]])
mean_diff <- mean(new_df$diff_vdp)
#Lower 95% confidence interval limits
l_lim <- mean_diff - 1.96*sd(new_df$diff_vdp)
#Upper 95% confidence interval limits
u_lim <- mean_diff + 1.96*sd(new_df$diff_vdp)
xlim_vdp <- 36
y_label <- expression(bold(VDP[N4] * "-" * VDP[FA] * " "  * "(%)"))
x_label <- expression(bold("(" * VDP[N4] * "+" * VDP[FA] * ")/2 " * "(%)"))
#Bland-Altman plot
theme_set(theme_bw())
vdp_bldaltman <- ggplot(new_df, aes(x = mean_vdp, y = diff_vdp)) +
  xlim(0, xlim_vdp) + ylim(-10, 12) +
  annotate("rect", xmin= -Inf, xmax= Inf, ymin=l_lim, ymax=u_lim,
           fill="red", alpha = 0.2) +
  geom_point(aes(color = Category), size=3, alpha=1) +
  scale_color_manual(values = c("#000000", "#009e73")) +
  geom_hline(yintercept = mean_diff, linewidth = 1) +
  geom_hline(yintercept = l_lim, color = "#000000", linetype=4, linewidth = 1.25) +
  geom_hline(yintercept = u_lim, color = "#000000", linetype=4, linewidth = 1.25) +
  ylab(y_label) +
  xlab(x_label) +
  guides(color = "none") +
  theme(axis.title.x = element_text(vjust = 0, size = 22, face = "bold"),
        axis.title.y = element_text(vjust = 0, size = 22, face = "bold"),
        axis.text = element_text(size = 22, face = "bold", color = "#000000"))
vdp_bldaltman
vdp_bldaltman_txt <- vdp_bldaltman + annotate("text", x=0.8*xlim_vdp, y=2*mean_diff,label=sprintf("Mean=%.3f", mean_diff),
                                                    size=5) + # ,fontface="bold"
  annotate("text", x=0.78*xlim_vdp, y=1.5*u_lim,label=sprintf("+1.96xSD=%.3f", u_lim),size=5) + # ,fontface="bold"
  annotate("text", x=0.78*xlim_vdp, y=1.5*l_lim,label=sprintf("-1.96xSD=%.3f", l_lim),size=5)
print(vdp_bldaltman_txt)
# Save the plot as a png file in the specified directory
ggsave("./zR_plots_4ppr/vdp_N4_keyhole_baltman_plain.png", plot = vdp_bldaltman, width = 5.2, height = 3.8, dpi = 300)
ggsave("./zR_plots_4ppr/vdp_N4_keyhole_baltman_p.png", plot = vdp_bldaltman_txt, width = 5.2, height = 3.8, dpi = 300)

###########Regression line plot: VDP threshold 
vdp_fit_parameters <- lm(FA_vdp ~ N4_vdp, data = new_df)
vdp_slope <- coef(vdp_fit_parameters)[2]
vdp_intercept <- coef(vdp_fit_parameters)[1]
x_label <- expression(bold(VDP[N4] * " " * "(%)"))
y_label <- expression(bold(VDP[FA] * " " * "(%)"))
linreg_vdp <- ggscatter(new_df, x = "N4_vdp", y = "FA_vdp",
                         add = "reg.line",  # Add regression line
                         xlim = c(0, 30),
                         ylim = c(0, 30), 
                         add.params = list(color = "#000000", size = 2, fill = "red", alpha = 0.15),
                         conf.int = TRUE, fullrange = TRUE,
                         cor.method = "pearson",
                         color = "Category", shape = 19, size = 4, # Points color, shape and size
                         ggtheme = theme_bw(),
                         font.x = c(22, "bold", "#000000"),
                         font.y = c(22, "bold", "#000000"),
                         font.tickslab = c(22, "bold", "#000000")) +
  xlab(x_label) +
  ylab(y_label) +
  font("caption", size = 12, color = "gray", face = "bold.italic") +
  scale_color_manual(values = c("#000000", "#009e73")) +
  geom_abline(intercept = vdp_intercept, slope = vdp_slope, linetype = "dotted", linewidth=2, color = "#000000") +
  # geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth=1, color = "blue") +
  theme(legend.position = "none")
linreg_vdp

fit_eq_vdp <- parse(text = sprintf("VDP[FA] == %.2f + %.2f * VDP[N4]", vdp_intercept, vdp_slope))
linreg_vdp_txt <- linreg_vdp + stat_cor(p.accuracy = .001, method = "pearson", cor.coef.name = c("r", "P"),
         aes(label = paste(after_stat(rr.label),after_stat(p.label), sep = "~`,`~")),
         label.x = 0, label.y = 23, size = 6) +
          annotate("text", x = 13, y = 28, label = fit_eq_vdp, color = "#000000", size = 6)
print(linreg_vdp_txt)

# Save the plot as a png file in the specified directory
ggsave("./zR_plots_4ppr/vdp_N4_keyhole_linreg_plain.png", plot = linreg_vdp, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/vdp_N4_keyhole_linreg_p.png", plot = linreg_vdp_txt, width = 4.5, height = 3.7, dpi = 300)

