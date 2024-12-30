##Provide data
## Set the working directory to the path where your CSV files are located
setwd("C:/Users/HUSDQ4/OneDrive - cchmc/cincy_work/all_projects_data_work/vdp_analysis/analysis_comparisons")

# ##Load the from CSV files
FA_cf_sensitivity <- read.csv("./IRC740H_visit1_spiral_cf/vdp_analysis_results_December2024/FA_combined_sensitivity.csv")
N4_cf_sensitivity <- read.csv("./IRC740H_visit1_spiral_cf/vdp_analysis_results_December2024/N4_combined_sensitivity.csv")

FA_cf_specificity <- read.csv("./IRC740H_visit1_spiral_cf/vdp_analysis_results_December2024/FA_combined_specificity.csv")
N4_cf_specificity <- read.csv("./IRC740H_visit1_spiral_cf/vdp_analysis_results_December2024/N4_combined_specificity.csv")

FA_cf_accuracy <- read.csv("./IRC740H_visit1_spiral_cf/vdp_analysis_results_December2024/FA_combined_accuracy.csv")
N4_cf_accuracy <- read.csv("./IRC740H_visit1_spiral_cf/vdp_analysis_results_December2024/N4_combined_accuracy.csv")


normality_test <- function(data, columns, alpha = 0.05) {
  # Check if all columns exist in the data
  missing_cols <- setdiff(columns, colnames(data))
  if (length(missing_cols) > 0) {
    stop(paste("These columns are missing in the data:", paste(missing_cols, collapse = ", ")))
  }
  
  # Initialize a list to store results
  results <- list()
  
  # Perform Shapiro-Wilk test for each column
  for (col in columns) {
    test_result <- shapiro.test(data[[col]])
    is_normal <- ifelse(test_result$p.value > alpha, "Yes", "No")
    results[[col]] <- list(
      W_statistic = test_result$statistic,
      p_value = test_result$p.value,
      is_normal = is_normal
    )
  }
  
  # Print results in a readable format
  cat("Shapiro-Wilk Normality Test Results:\n")
  cat("------------------------------------\n")
  for (col in columns) {
    cat(paste0("Column: ", col, "\n"))
    cat(paste0("  W Statistic: ", round(results[[col]]$W_statistic, 4), "\n"))
    cat(paste0("  p-value: ", round(results[[col]]$p_value, 4), "\n"))
    cat(paste0("  Normal Data: ", results[[col]]$is_normal, "\n\n"))
  }
  
  # Return results as a list
  return(results)
}

normality_test(df, c("col1", "col2", "col3"))

