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

##Define the function to perform normality tests on all columns of dataframe
normality_test <- function(data, columns, alpha = 0.05) {
  # Check if all columns exist in the data
  missing_cols <- setdiff(columns, colnames(data))
  if (length(missing_cols) > 0) {
    stop(paste("These columns are missing in the data:", paste(missing_cols, collapse = ", ")))
  }
  
  # Initialize a dataframe to store results
  results <- data.frame(
    Column = character(),
    W_Statistic = numeric(),
    P_Value = numeric(),
    Normal_Data = character(),
    stringsAsFactors = FALSE
  )
  
  # Perform Shapiro-Wilk test for each column
  for (col in columns) {
    test_result <- shapiro.test(data[[col]])
    is_normal <- ifelse(test_result$p.value > alpha, "Yes", "No")
    results <- rbind(results, data.frame(
      Column = col,
      W_Statistic = test_result$statistic,
      P_Value = test_result$p.value,
      Normal_Data = is_normal
    ))
  }
  
  # Print results in a readable format
  cat("Shapiro-Wilk Normality Test Results:\n")
  cat("------------------------------------\n")
  print(results, row.names = FALSE)
  
  # Return results as a dataframe
  return(results)
}


FA_sensitivity <- normality_test(FA_cf_sensitivity, c("Median", "Percentile", "Adaptive", "Hierarchical",
                                    "Thresholding", "Mean"))
N4_sensitivity <- normality_test(N4_cf_sensitivity, c("Median", "Percentile", "Adaptive", "Hierarchical",
                                                      "Thresholding", "Mean"))

FA_specificity <- normality_test(FA_cf_specificity, c("Median", "Percentile", "Adaptive", "Hierarchical",
                                                      "Thresholding", "Mean"))
N4_specificity <- normality_test(N4_cf_specificity, c("Median", "Percentile", "Adaptive", "Hierarchical",
                                                      "Thresholding", "Mean"))

FA_accuracy <- normality_test(FA_cf_accuracy, c("Median", "Percentile", "Adaptive", "Hierarchical",
                                                      "Thresholding", "Mean"))
N4_accuracy <- normality_test(N4_cf_accuracy, c("Median", "Percentile", "Adaptive", "Hierarchical",
                                                      "Thresholding", "Mean"))


