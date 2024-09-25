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

#Colorblind friendly palette:
cbPalette <- c("#999999","#E69F00","#CC79A7","#F0E442","#56B4E9","#009E73","#0072B2","#D55E00")
# Color order to use
# plt_order <- c("Thresholding","Hierarchical","Adaptive","Percentile","Median","Reader")

###Function for plotting connected boxplot
connected_bxp <- function(data_in, id = "Subject_id", x_var, y_var, ylim,
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

################################################################################
# ##Load the data from both CSV files and arrange it to plot
N4_gre_cf_vdps <- read.csv("./HyPOINT_phase1_cart_cf/vdp_analysis_results_September2024/N4_combined_vdps.csv")

##Reshape Data to Long Format for ggplot2
N4_long <- N4_gre_cf_vdps %>%
  pivot_longer(cols = -Subject_id, names_to = "AnalysisMethod", values_to = "VDP")

##Select data for only "Reader" and one analysis plot
N4_rdr_thrsh <- N4_long %>% filter(AnalysisMethod %in% c("Reader", "Thresholding"))
N4_rdr_hrar <- N4_long %>% filter(AnalysisMethod %in% c("Reader", "Hierarchical"))
N4_rdr_adpt <- N4_long %>% filter(AnalysisMethod %in% c("Reader", "Adaptive"))
N4_rdr_prcnt <- N4_long %>% filter(AnalysisMethod %in% c("Reader", "Percentile"))
N4_rdr_medn <- N4_long %>% filter(AnalysisMethod %in% c("Reader", "Median"))

##Define y-label of the plot
y_label <- expression(bold(VDP[GRE-N4]* "(%)"))

##############################################
# ##Statistical test
# vdp_stat.test <- N4_rdr_thresh %>%
#   wilcox_test(VDP ~ AnalysisMethod , paired = TRUE) %>%
#   add_significance()
# vdp_stat.test
##############################################
thresh_vdp_bxp <- connected_bxp(N4_rdr_thrsh, "Subject_id", "AnalysisMethod",
                                "VDP", c(0, 35), cbPalette[c(6, 1)], "", y_label)
hrar_vdp_bxp <- connected_bxp(N4_rdr_hrar, "Subject_id", "AnalysisMethod", "VDP",
                                c(0, 35), cbPalette[c(6, 2)], "", y_label)
adpt_vdp_bxp <- connected_bxp(N4_rdr_adpt, "Subject_id", "AnalysisMethod", "VDP",
                                c(0, 35), cbPalette[c(6, 3)], "", y_label)
prcnt_vdp_bxp <- connected_bxp(N4_rdr_prcnt, "Subject_id", "AnalysisMethod", "VDP",
                                c(0, 70), cbPalette[c(6, 4)], "", y_label)
medn_vdp_bxp <- connected_bxp(N4_rdr_medn, "Subject_id", "AnalysisMethod", "VDP",
                                c(0, 35), cbPalette[c(6, 5)], "", y_label)





# Box plots with p-values SNR
vdpmed_bxp <-  ggpaired(N4_rdr_thresh, x = "Analysis_Method", y = "VDP",
                        fill = "Analysis_Method", palette = c("#e69f00", "#56b4e9"),
                        width = 0.5, ylim = c(0, 35), line.color = "Analysis_Method",
                        line.size = 0.5,legend = "none", xlab = "") + # Sequence
  ylab(y_label) +
  # geom_point(aes(color = Category), shape = 19, size = 3) +
  scale_color_manual(values = c("#000000", "#009e73")) +
  theme(panel.border = element_rect(color = "#000000", fill = NA, linewidth = 1),
        axis.text = element_text(size = 22, color = "#000000", face = "bold"),
        axis.title = element_text(size = 22, color = "#000000", face = "bold"),
        axis.line.x = element_line(linewidth = 1), axis.line.y = element_line(linewidth = 1)
  )
vdpmed_bxp
vdp_stat.test <- vdp_stat.test %>% add_xy_position(x = "Analysis_Method")
vdpmed_bxp_p <-  vdpmed_bxp + stat_pvalue_manual(vdp_stat.test, label = "p{scales::pvalue(p)}",
                                         y.position = 29, label.size = 8, bracket.size = 0.8,
                                         tip.length = 0.03, vjust=-0.35)
vdpmed_bxp_p
# Save the plot as a png file in the specified directory
ggsave("./zR_plots_4ppr/cartVSspir_vdp_bxp_nop.png", plot = vdpmed_bxp, width = 4.5, height = 3.7, dpi = 300)
ggsave("./zR_plots_4ppr/cartVSspir_vdp_bxp_p.png", plot = vdpmed_bxp_p, width = 4.5, height = 3.7, dpi = 300)



