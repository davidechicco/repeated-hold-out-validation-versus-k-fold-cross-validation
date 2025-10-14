# Load necessary libraries

library("pacman")
p_load("randomForest", "caret")
# set.seed(11)  # For reproducibility

source("utils.r")

k <- 15

cat(k, "-fold cross-validation\t Random Forests regression analysis\n", sep="")

# Load the data
# dataFile <- "../../clustering_internal_metrics/data/real_datasets/journal_pone_0161784_preterm_sepsis_edited_dataset_IMPUTED.csv"
dataFile <- "../data/sepsis_severity_dataset_col_norm_edited_target-SOFA-score.csv"
data <- read.csv(dataFile, header=T)
cat(dataFile,"\n")

# colnames(data[,ncol(data)]) <- "y"

shuffled_df <- data[sample(nrow(data)), ]
data <- shuffled_df

cat("dataset dimensions: ", nrow(data), " rows and ", ncol(data), " columns\n", sep="")

missing_counts <- colSums(is.na(data))

# cat("missing data: ", missing_counts, "\n", sep="")

# Check the structure of the data
#str(data)

# Separate features and target variable
# Assuming the last column is the target variable
target_variable <- data[, ncol(data)]
features <- data[, -ncol(data)]

# Create a trainControl object for k-fold cross-validation
train_control <- trainControl(method = "cv", number = k)

# Fit the Random Forest model using caret
model <- train(SOFA.score ~ ., data = data, method = "rf", trControl = train_control)

# Make predictions on the training set
predictions <- predict(model, data)

# Calculate R-squared
r_squared <- cor(target_variable, predictions)^2
cat("R-squared:", r_squared, "\n")

computeExecutionTime()
