# Load necessary libraries

library("pacman")
p_load("randomForest", "caret", "MLmetrics", "dtw")
set.seed(123)  # For reproducibility

source("utils.r")

VERBOSE <- FALSE

# Load the data
# dataFile <- "../data/CKD_CVD_journal.pone.0199920.s002_EDITED_IMPUTED_v2.csv"
# response_var <- "TimeToEventMonths"

# dataFile <- "../data/sepsis_severity_dataset_col_norm_edited_target-SOFA-score.csv"
# response_var <- "SOFA.score"  # Replace with your response variable name

# dataFile <- "../data/EHRs_heart_failure_S1Data_EDITED_v3.csv"
# response_var <- "TIME_DAYS"  # Replace with your response variable name

dataFile <- "../data/ObesityDataSet_raw_and_data_sinthetic_EDITED.csv"
response_var <- "obesity_class"  # Replace with your response variable name

# dataFile <- "../data/journal.pone.0216416_Takashi2019_diabetes_type1_dataset_preprocessed.csv"
# response_var <- "duration.of.diabetes"

data <- read.csv(dataFile, header=T)
cat(dataFile,"\n")

data %>% dim() %>% print()


data <- data[sample(nrow(data)), ]

# Ensure the target variable is a factor
# data$SOFA.score <- as.factor(data$"SOFA.score")
predictor_vars <- setdiff(names(data), response_var)

# Set parameters for repeated hold-out validation
num_iterations <- 10
train_fraction <- 0.9
cat("training set = ", dec_five(100*train_fraction), "%\n", sep="")
r_squared_values <- numeric(num_iterations)
r_squared_values_V2 <- numeric(num_iterations)
r_squared_values_V3 <- numeric(num_iterations)
target_differences_perc_means <- numeric(num_iterations)

cat(num_iterations, "-times repeated hold-out validation\t Random Forests regression analysis\n", sep="")

cat("global target mean ", mean(data[[response_var]]),  " +- ", sd(data[[response_var]]), "\n", sep="")

mean_differences_between_global_targets_and_fold_targets <- c()
dtw_mean_similarities <- 0
interp_mean_similarities <- 0


global_interations <- 10
cat("global_interations = ", global_interations, "\n", sep="")

for(a in seq(1:global_interations)) {

    dtw_similarities <- c()
    interp_similarities <- c()

    for (i in 1:num_iterations) {

      shuffled_df <- data[sample(nrow(data)), ]
      data <- shuffled_df

      # Split data into training and testing sets
      train_indices <- sample(1:nrow(data), size = floor(train_fraction * nrow(data)))
      train_data <- data[train_indices, ]
      test_data <- data[-train_indices, ]

      # Fit Random Forest model
      rf_model <- randomForest(as.formula(paste(response_var, paste(predictor_vars, collapse = "+"), sep = " ~ ")),
                              data = train_data)

      # Make predictions on the test set
      predictions <- predict(rf_model, newdata = test_data)

      # Calculate R-squared
      actuals <- test_data[[response_var]]

      if(VERBOSE) cat("a=", a, ", i=", i, ") target mean ", mean(actuals),  " +- ", sd(actuals), "\t", sep="")

      if(VERBOSE) cat(" mean diff with original = ")
      thisDiff <- abs(computeDiffPerc(mean(data[[response_var]]),mean(actuals)))
      if(VERBOSE) cat(dec_five(thisDiff), "%\n", sep="")

      dtw_alignment <- dtw(data[[response_var]], actuals)
      dtw_similarity <- 1 / (1 + dtw_alignment$"distance")  # Convert distance to similarity
      # DTW similarity (common transforms): values between 0 (very dissimilar) and 1 (identical)

      # Interpolate y to length of x
      y_interp <- approx(x = seq_along(actuals), y = actuals, xout = seq(1, length(actuals), length.out = length(data[[response_var]])))$"y"

      # Compute correlation
      interp_similarity <- cor(data[[response_var]], y_interp)

      r_squared <- 1 - (sum((actuals - predictions) ^ 2) / sum((actuals - mean(actuals)) ^ 2))
      # cat("(i=", i,") R-squared = ", r_squared, "\n", sep="" )
      r_squared_V2 <- cor(actuals, predictions)^2
      # cat("(i=", i,") R-squared V2:", r_squared_V2, "\n", sep="")
      r_squared_V3 <- R2_Score(predictions, actuals)

      # Store the R-squared value
      r_squared_values[i] <- r_squared
      r_squared_values_V2[i] <- r_squared_V2
      r_squared_values_V3[i] <- r_squared_V3
      target_differences_perc_means[i] <- thisDiff
      dtw_similarities[i] <- dtw_similarity
      interp_similarities[i] <- interp_similarity

    }

    # Calculate average R-squared across all folds
    # average_r_squared <- mean(r_squared_values)
    # Print the average R-squared
    # cat("Average R-squared V1 over", k, "folds:", average_r_squared, "\n")
    # average_r_squared_V2 <- mean(r_squared_values_V2)
    # cat("Average R-squared V2 over", k, "folds:", average_r_squared_V2, "\n")
    average_r_squared_V3 <- mean(r_squared_values_V3)
    sd_r_squared_V3 <- sd(r_squared_values_V3)
    if(VERBOSE) cat("Average R-squared V3 over ", num_iterations, " iterations: ", average_r_squared_V3, " +- ", sd_r_squared_V3, "  ", sep="")
    if(VERBOSE) cat(" in the [", dec_five(min(r_squared_values_V3)), ",", dec_five(max(r_squared_values_V3)), "] interval\n", sep="")

    # cat("(a=", a, ") absolute average absolute percentage difference between iteration target and original target = ", dec_five(mean(target_differences_perc_means)), "% ± ", dec_five(sd(target_differences_perc_means)), "%\n", sep="")
    # cat("(a=", a, ") mean DTW similarity between fold targets and global targets = ", dec_five(mean(dtw_similarities)), " ± ", dec_five(sd(dtw_similarities)), "\n", sep="")
    # cat("(a=", a, ") mean interpolation similarity between fold targets and global targets = ", dec_five(mean(interp_similarities)), " ± ", dec_five(sd(interp_similarities)), "\n", sep="")

    mean_differences_between_global_targets_and_fold_targets[a] <- mean(target_differences_perc_means)
    dtw_mean_similarities[a] <- mean(dtw_similarities)
    interp_mean_similarities[a] <- mean(interp_similarities)

}

cat("\nglobal absolute average absolute percentage difference between iteration target and original target = ", dec_five(mean(mean_differences_between_global_targets_and_fold_targets)), "% ± ", dec_five(sd(mean_differences_between_global_targets_and_fold_targets)), "%\n", sep="")
# cat("\nglobal DTW similarities between iteration target and original target = ", dec_five(mean(dtw_mean_similarities)), " ± ", dec_five(sd(dtw_mean_similarities)), "\n", sep="")
cat("\nglobal interpolation similarities between iteration target and original target = ", dec_five(mean(interp_mean_similarities)), " ± ", dec_five(sd(interp_mean_similarities)), "\n", sep="")


cat(num_iterations, "-times repeated hold-out validation\t Random Forests regression analysis\n", sep="")


computeExecutionTime()
