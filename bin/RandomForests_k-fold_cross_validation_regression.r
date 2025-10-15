# Load necessary libraries

library("pacman")
p_load("randomForest", "caret", "MLmetrics")
# set.seed(11)  # For reproducibility

source("utils.r")

k <- 10

cat(k, "-fold cross-validation\t Random Forests regression analysis\n", sep="")

# Load the data
# dataFile <- "../data/CKD_CVD_journal.pone.0199920.s002_EDITED_IMPUTED_v2.csv"
# response_var <- "TimeToEventMonths"

# dataFile <- "../data/sepsis_severity_dataset_col_norm_edited_target-SOFA-score.csv"
# response_var <- "SOFA.score"  # Replace with your response variable name

# dataFile <- "../data/EHRs_heart_failure_S1Data_EDITED_v3.csv"
# response_var <- "TIME_DAYS"  # Replace with your response variable name

# dataFile <- "../data/ObesityDataSet_raw_and_data_sinthetic_EDITED.csv"
# response_var <- "obesity_class"  # Replace with your response variable name

# dataFile <- "../data/journal.pone.0216416_Takashi2019_diabetes_type1_dataset_preprocessed.csv"
# response_var <- "duration.of.diabetes"

data <- read.csv(dataFile, header=T)
cat(dataFile,"\n")

data <- data[sample(nrow(data)), ]

# Ensure the target variable is a factor
# data$SOFA.score <- as.factor(data$"SOFA.score")
predictor_vars <- setdiff(names(data), response_var)


cat("global target mean ", mean(data[[response_var]]),  " +- ", sd(data[[response_var]]), "\n", sep="")

# Set parameters for repeated hold-out validation
set.seed(123)  # For reproducibility
train_fraction <- 1-1/k
cat("k = ", k,  " somehow like training set = ", dec_two(100*train_fraction), "%\n", sep="")
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

  cat("(i=", i, " of ", k, ") ", sep="")
  cat("target mean ", mean(actuals),  " +- ", sd(actuals), " ", sep="")

  cat(" mean diff with original = ")
  thisDiff <- computeDiffPerc(mean(data[[response_var]]),mean(actuals))
  cat(dec_two(thisDiff), "%\n", sep="")


  # Store the R-squared value
  r_squared_values[i] <- r_squared
  r_squared_values_V2[i] <- r_squared_V2
  r_squared_values_V3[i] <- r_squared_V3
  target_differences_perc_means[i] <- abs(thisDiff)

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
cat("Average R-squared V3 over ",  k, " folds: ", average_r_squared_V3, " +- ", sd_r_squared_V3, "  ", sep="")
cat(" in the [", dec_two(min(r_squared_values_V3)), ",", dec_two(max(r_squared_values_V3)), "] interval\n", sep="")

cat("average absolute percentage difference between iteration target and original target = ", dec_two(mean(target_differences_perc_means)), "%\n", sep="")

cat(k, "-fold cross-validation\t Random Forests regression analysis\n", sep="")

computeExecutionTime()
