# Load necessary libraries

library("pacman")
p_load("randomForest", "caret")
set.seed(11)  # For reproducibility

source("utils.r")

k <- 10

cat(k, "-fold cross-validation\t Random Forests\n", sep="")

# Load the data
# dataFile <- "../../clustering_internal_metrics/data/real_datasets/journal_pone_0161784_preterm_sepsis_edited_dataset_IMPUTED.csv"
dataFile <- "../data/EHRs_heart_failure_S1Data_EDITED_v3.csv"
data <- read.csv(dataFile)
cat(dataFile,"\n")


# Ensure the target variable is a factor
data[, ncol(data)] <- as.factor(data[, ncol(data)])

cat("dataset dimensions: ", nrow(data), " rows and ", ncol(data), " columns\n")

# Set up k-fold cross-validation (e.g., k = 10)
control <- trainControl(method = "cv", number = k)

# Train the Random Forest model
model <- train(data[, -ncol(data)], data[, ncol(data)],
               method = "rf",
               trControl = control)

# Make predictions using cross-validated folds
predictions <- predict(model, data[, -ncol(data)])

# Create confusion matrix
confusion <- confusionMatrix(predictions, data[, ncol(data)])

# Calculate Matthews correlation coefficient (MCC)
mcc <- (confusion$table[1, 1] * confusion$table[2, 2] -
         confusion$table[1, 2] * confusion$table[2, 1]) /
        sqrt((confusion$table[1, 1] + confusion$table[1, 2]) *
             (confusion$table[1, 1] + confusion$table[2, 1]) *
             (confusion$table[2, 2] + confusion$table[1, 2]) *
             (confusion$table[2, 2] + confusion$table[2, 1]))

# Print the Matthews correlation coefficient
cat("Matthews Correlation Coefficient (MCC):", mcc, "\n")

computeExecutionTime()
