# Load necessary libraries
library("pacman")
p_load("randomForest", "caret", "MLmetrics", "dtw")
set.seed(123)  # For reproducibility

source("utils.r")

k <- 10
VERBOSE <- FALSE


cat(k, "-fold cross-validation\t Random Forests regression analysis\n", sep="")

# Load the data
# dataFile <- "../data/CKD_CVD_journal.pone.0199920.s002_EDITED_IMPUTED_v2.csv"
# response_var <- "TimeToEventMonths"

dataFile <- "../data/sepsis_severity_dataset_col_norm_edited_target-SOFA-score.csv"
response_var <- "SOFA.score"  # Replace with your response variable name

# dataFile <- "../data/EHRs_heart_failure_S1Data_EDITED_v3.csv"
# response_var <- "TIME_DAYS"  # Replace with your response variable name

# dataFile <- "../data/ObesityDataSet_raw_and_data_sinthetic_EDITED.csv"
# response_var <- "obesity_class"  # Replace with your response variable name

# dataFile <- "../data/journal.pone.0216416_Takashi2019_diabetes_type1_dataset_preprocessed.csv"
# response_var <- "duration.of.diabetes"

data <- read.csv(dataFile, header=T)
cat(dataFile,"\n")

data %>% dim() %>% print()

data <- data[sample(nrow(data)), ]

# Ensure the target variable is a factor
# data$SOFA.score <- as.factor(data$"SOFA.score")
predictor_vars <- setdiff(names(data), response_var)

cat("global target mean ", mean(data[[response_var]]),  " ± ", sd(data[[response_var]]), "\n", sep="")

crossval_mean_differences_between_global_targets_and_fold_targets <- c()

crossval_dtw_global_interations <- 20
cat("crossval_dtw_global_interations = ", crossval_dtw_global_interations, "\n", sep="")

cross_dtw_mean_similarities <- 0
crossval_interp_mean_similarities_Pearson <- 0
crossval_interp_mean_similarities_Kendall <- 0
crossval_interp_mean_similarities_Spearman <- 0


for(a in seq(1:crossval_dtw_global_interations)) {

    crossval_dtw_similarities <- c()
    crossval_interp_similarities_Pearson <- c()
    crossval_interp_similarities_Kendall <- c()
    crossval_interp_similarities_Spearman <- c()

    # Set parameters for repeated hold-out validation
    train_fraction <- 1-1/k
    if(VERBOSE) cat("k = ", k,  " somehow like training set = ", dec_five(100*train_fraction), "%\n", sep="")
    r_squared_values <- numeric(k)
    r_squared_values_V2 <- numeric(k)
    r_squared_values_V3 <- numeric(k)

    target_differences_perc_means <- numeric(k)

    folds <- sample(1:k, nrow(data), replace = TRUE)

    # Perform k-fold cross-validation
    for (i in 1:k) {
      # Split data into training and testing sets
      train_data <- data[folds != i, ]
      test_data <- data[folds == i, ]

      # cat("train_data #rows(): ", nrow(train_data), "/", nrow(data), "\t", round(nrow(train_data)*100/nrow(data)), "%\t", sep="")

      # Fit Random Forest model
      rf_model <- randomForest(as.formula(paste(response_var, paste(predictor_vars, collapse = "+"), sep = " ~ ")),
                              data = train_data)

      # Make predictions on the test set
      predictions <- predict(rf_model, newdata = test_data)

      # Calculate R-squared
      actuals <- test_data[[response_var]]
      r_squared <- 1 - (sum((actuals - predictions) ^ 2) / sum((actuals - mean(actuals)) ^ 2))
      r_squared_V2 <- cor(actuals, predictions)^2
      r_squared_V3 <- R2_Score(predictions, actuals)

      if(VERBOSE) cat("a = ", a, ",i=", i, " of ", k, " ", sep="")
      if(VERBOSE) cat("target mean ", mean(actuals),  " ± ", sd(actuals), " ", sep="")

      if(VERBOSE) cat(" mean diff with original = ")
      thisDiff <- abs(computeDiffPerc(mean(data[[response_var]]),mean(actuals)))
      if(VERBOSE) cat(dec_five(thisDiff), "%\n", sep="")

      crossval_dtw_alignment <- dtw(data[[response_var]], actuals)
      crossval_dtw_similarity <- 1 / (1 + crossval_dtw_alignment$"distance")  # Convert distance to similarity
      # DTW similarity (common transforms): values between 0 (very dissimilar) and 1 (identical)

      # crossval_interpolate y to length of x
      y_crossval_interp <- approx(x = seq_along(actuals), y = actuals, xout = seq(1, length(actuals), length.out = length(data[[response_var]])))$"y"

      # Compute correlation
      crossval_interp_similarity_Pearson <- cor(data[[response_var]], y_crossval_interp, method = "pearson")
      crossval_interp_similarity_Kendall <- cor(data[[response_var]], y_crossval_interp, method = "kendall")
      crossval_interp_similarity_Spearman <- cor(data[[response_var]], y_crossval_interp, method = "spearman")


      # Store the R-squared value
      r_squared_values[i] <- r_squared
      r_squared_values_V2[i] <- r_squared_V2
      r_squared_values_V3[i] <- r_squared_V3
      target_differences_perc_means[i] <- thisDiff
      crossval_dtw_similarities[i] <- crossval_dtw_similarity
      crossval_interp_similarities_Pearson[i] <- crossval_interp_similarity_Pearson
      crossval_interp_similarities_Kendall[i] <- crossval_interp_similarity_Kendall
      crossval_interp_similarities_Spearman[i] <- crossval_interp_similarity_Spearman

      # cat("r_squared = ", r_squared, "\t r_squared_V2 =  ", r_squared_V2, "\n", sep="")

    }

    # Calculate average R-squared across all folds
    # average_r_squared <- mean(r_squared_values)
    # Print the average R-squared
    # cat("Average R-squared V1 over", k, "folds:", average_r_squared, "\n")
    # average_r_squared_V2 <- mean(r_squared_values_V2)
    # cat("Average R-squared V2 over", k, "folds:", average_r_squared_V2, "\n")
    average_r_squared_V3 <- mean(r_squared_values_V3)
    sd_r_squared_V3 <- sd(r_squared_values_V3)
    if(VERBOSE) cat("Average R-squared V3 over ",  k, " folds: ", average_r_squared_V3, " ± ", sd_r_squared_V3, "  ", sep="")
    if(VERBOSE) cat(" in the [", dec_five(min(r_squared_values_V3)), ",", dec_five(max(r_squared_values_V3)), "] interval\n", sep="")

    # cat("(a=", a, ") absolute average absolute percentage difference between iteration target and original target = ", dec_three(mean(target_differences_perc_means)), "% ± ", dec_three(sd(target_differences_perc_means)), "%\n", sep="")
    # cat("(a=", a, ") mean DTW similarity between fold targets and global targets = ", dec_five(mean(crossval_dtw_similarities)), " ± ", dec_five(sd(crossval_dtw_similarities)), "\n", sep="")
    # cat("(a=", a, ") mean crossval_interpolation similarity between fold targets and global targets = ", dec_five(mean(crossval_interp_similarities_Pearson)), " ± ", dec_five(sd(crossval_interp_similarities_Pearson)), "\n", sep="")

    crossval_mean_differences_between_global_targets_and_fold_targets[a] <- mean(target_differences_perc_means)
    cross_dtw_mean_similarities[a] <- mean(crossval_dtw_similarities)
    crossval_interp_mean_similarities_Pearson[a] <- mean(crossval_interp_similarities_Pearson)
    crossval_interp_mean_similarities_Kendall[a] <- mean(crossval_interp_similarities_Kendall)
    crossval_interp_mean_similarities_Spearman[a] <- mean(crossval_interp_similarities_Spearman)

}


cat("\n[mean] global absolute average absolute percentage difference between iteration target and original target = ", dec_three(mean(crossval_mean_differences_between_global_targets_and_fold_targets)), "% ± ", dec_three(sd(crossval_mean_differences_between_global_targets_and_fold_targets)), "% ↑↑↑\n", sep="")
cat("[DTW] global DTW similarities between iteration target and original target = ", dec_five(mean(cross_dtw_mean_similarities)), " ± ", dec_five(sd(cross_dtw_mean_similarities)), " ↓↓↓\n", sep="")
cat("[Pearson crossval_interpolation] global crossval_interpolated similarities between iteration target and original target = ", dec_five(mean(crossval_interp_mean_similarities_Pearson)), " ± ", dec_five(sd(crossval_interp_mean_similarities_Pearson)), " ↑↑↑\n", sep="")
cat("[Kendal crossval_interpolation] global crossval_interpolated similarities between iteration target and original target = ", dec_five(mean(crossval_interp_mean_similarities_Kendall)), " ± ", dec_five(sd(crossval_interp_mean_similarities_Kendall)), " ↑↑↑\n", sep="")
cat("[Spearman crossval_interpolation] global crossval_interpolated similarities between iteration target and original target = ", dec_five(mean(crossval_interp_mean_similarities_Spearman)), " ± ", dec_five(sd(crossval_interp_mean_similarities_Spearman)), " ↑↑↑\n", sep="")

cat(k, "-fold cross-validation\t Random Forests regression analysis\n", sep="")

computeExecutionTime()
