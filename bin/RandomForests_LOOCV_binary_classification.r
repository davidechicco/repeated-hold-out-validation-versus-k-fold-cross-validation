# Load necessary libraries

library("pacman")
p_load("randomForest", "caret")

cat("LOOCV\t Random Forests\n")
cat("neuroblastoma Hiyama2010 GSE16237 gene expression\n")

# Load the data
data <- read.csv("../../clustering_internal_metrics/data/real_datasets/microarray_gene_expression_Hiyama2010_GSE16237_neuroblastoma_DCangelosi_signature_probesets_dataset_3791.csv")


# Ensure the target variable is a factor
data[, ncol(data)] <- as.factor(data[, ncol(data)])

# Set up leave-one-out cross-validation
control <- trainControl(method = "LOOCV")

# Train the Random Forest model
set.seed(123)  # For reproducibility
model <- train(data[, -ncol(data)], data[, ncol(data)],
               method = "rf",
               trControl = control)

# Print model details
# print(model)

# Make predictions
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
