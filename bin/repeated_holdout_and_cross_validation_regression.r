setwd(".")
options(stringsAsFactors = FALSE)
cat("\014")
set.seed(123)
options(repos = list(CRAN="http://cran.rstudio.com/"))

list.of.packages <- c("pacman")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library("pacman")
p_load("randomForest", "caret", "MLmetrics", "dtw", "crayon")

source("utils.r")

VERBOSE <- FALSE

cat(" =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  = \n")


# rep_holdout_global_interations <- 5000

args <- commandArgs(trailingOnly = TRUE)  # get only the user-supplied arguments

if (length(args) < 2) {
  stop("Please provide two arguments")
}

rep_holdout_global_interations <- args[1]
dataset_name <- args[2]

cat("argument 1, rep_holdout_global_interations:", rep_holdout_global_interations, "\n")
cat("argument 2, dataset:", dataset_name, "\n")

# Load the data

if(dataset_name == "chronic_kidney_disease") {
  dataFile <- "../data/CKD_CVD_journal.pone.0199920.s002_EDITED_IMPUTED_v2.csv"
  response_var <- "TimeToEventMonths"
} else if(dataset_name == "sepsis") {
  dataFile <- "../data/sepsis_severity_dataset_col_norm_edited_target-SOFA-score.csv"
  response_var <- "SOFA.score"  # Replace with your response variable name
} else if(dataset_name == "heart_failure") {
  dataFile <- "../data/EHRs_heart_failure_S1Data_EDITED_v3.csv"
  response_var <- "TIME_DAYS"  # Replace with your response variable name
} else if(dataset_name == "obesity") {
  dataFile <- "../data/ObesityDataSet_raw_and_data_sinthetic_EDITED.csv"
  response_var <- "obesity_class"  # Replace with your response variable name
} else if(dataset_name == "diabetes_type_one") {
  dataFile <- "../data/journal.pone.0216416_Takashi2019_diabetes_type1_dataset_preprocessed.csv"
  response_var <- "duration.of.diabetes"
} else {
  cat("Error: the dataset chosen should be one among chronic_kidney_disease, sepsis, heart_failure, obesity, diabetes_type_one. The program stops here\n")
  quit(save = "no", status = 0)
}


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
rep_holdout_target_differences_perc_means <- numeric(num_iterations)

cat(num_iterations, "-times repeated hold-out validation\t Random Forests regression analysis\n", sep="")

cat("global target mean ", mean(data[[response_var]]),  " ± ", sd(data[[response_var]]), "\n", sep="")

rep_holdout_mean_differences_between_global_targets_and_fold_targets <- c()
rep_holdout_dtw_mean_similarities <- 0
rep_holdout_interp_mean_similarities_Pearson <- 0
rep_holdout_interp_mean_similarities_Kendall <- 0
rep_holdout_interp_mean_similarities_Spearman <- 0

cat("rep_holdout_global_interations = ", rep_holdout_global_interations, "\n", sep="")

for(a in seq(1:rep_holdout_global_interations)) {

    rep_holdout_dtw_similarities <- c()
    rep_holdout_interp_similarities_Pearson <- c()
    rep_holdout_interp_similarities_Kendall <- c()
    rep_holdout_interp_similarities_Spearman <- c()

    for (i in 1:num_iterations) {

      shuffled_df <- data[sample(nrow(data)), ]
      data <- shuffled_df

      # Split data into training and testing sets
      train_indices <- sample(1:nrow(data), size = floor(train_fraction * nrow(data)))
      train_data <- data[train_indices, ]
      test_data <- data[-train_indices, ]

      # Calculate R-squared
      actuals <- test_data[[response_var]]

      if(VERBOSE) cat("a=", a, ", i=", i, ") target mean ", mean(actuals),  " ± ", sd(actuals), "\t", sep="")

      if(VERBOSE) cat(" mean diff with original = ")
      thisDiff <- abs(computeDiffPerc(mean(data[[response_var]]),mean(actuals)))
      if(VERBOSE) cat(dec_five(thisDiff), "%\n", sep="")

      rep_holdout_dtw_alignment <- dtw(data[[response_var]], actuals)
      rep_holdout_dtw_similarity <- 1 / (1 + rep_holdout_dtw_alignment$"distance")  # Convert distance to similarity
      # DTW similarity (common transforms): values between 0 (very dissimilar) and 1 (identical)

      # rep_holdout_interpolate y to length of x
      y_rep_holdout_interp <- approx(x = seq_along(actuals), y = actuals, xout = seq(1, length(actuals), length.out = length(data[[response_var]])))$"y"

      # Compute correlation
      rep_holdout_interp_similarity_Pearson <- cor(data[[response_var]], y_rep_holdout_interp, method = "pearson")
      rep_holdout_interp_similarity_Kendall <- cor(data[[response_var]], y_rep_holdout_interp, method = "kendall")
      rep_holdout_interp_similarity_Spearman <- cor(data[[response_var]], y_rep_holdout_interp, method = "spearman")

      rep_holdout_target_differences_perc_means[i] <- thisDiff
      rep_holdout_dtw_similarities[i] <- rep_holdout_dtw_similarity
      rep_holdout_interp_similarities_Pearson[i] <- rep_holdout_interp_similarity_Pearson
      rep_holdout_interp_similarities_Kendall[i] <- rep_holdout_interp_similarity_Kendall
      rep_holdout_interp_similarities_Spearman[i] <- rep_holdout_interp_similarity_Spearman

    }

    # Calculate average R-squared across all folds
    # average_r_squared <- mean(r_squared_values)
    # Print the average R-squared
    # cat("Average R-squared V1 over", k, "folds:", average_r_squared, "\n")
    # average_r_squared_V2 <- mean(r_squared_values_V2)
    # cat("Average R-squared V2 over", k, "folds:", average_r_squared_V2, "\n")
#     average_r_squared_V3 <- mean(r_squared_values_V3)
#     sd_r_squared_V3 <- sd(r_squared_values_V3)
#     if(VERBOSE) cat("Average R-squared V3 over ", num_iterations, " iterations: ", average_r_squared_V3, " ± ", sd_r_squared_V3, "  ", sep="")
#     if(VERBOSE) cat(" in the [", dec_five(min(r_squared_values_V3)), ",", dec_five(max(r_squared_values_V3)), "] interval\n", sep="")

    # cat("(a=", a, ") absolute average absolute percentage difference between iteration target and original target = ", dec_five(mean(rep_holdout_target_differences_perc_means)), "% ± ", dec_five(sd(rep_holdout_target_differences_perc_means)), "%\n", sep="")
    # cat("(a=", a, ") mean DTW similarity between fold targets and global targets = ", dec_five(mean(rep_holdout_dtw_similarities)), " ± ", dec_five(sd(rep_holdout_dtw_similarities)), "\n", sep="")
    # cat("(a=", a, ") mean rep_holdout_interpolation similarity between fold targets and global targets = ", dec_five(mean(rep_holdout_interp_similarities)), " ± ", dec_five(sd(rep_holdout_interp_similarities)), "\n", sep="")

    rep_holdout_mean_differences_between_global_targets_and_fold_targets[a] <- mean(rep_holdout_target_differences_perc_means)
    rep_holdout_dtw_mean_similarities[a] <- mean(rep_holdout_dtw_similarities)
    rep_holdout_interp_mean_similarities_Pearson[a] <- mean(rep_holdout_interp_similarities_Pearson)
    rep_holdout_interp_mean_similarities_Kendall[a] <- mean(rep_holdout_interp_similarities_Kendall)
    rep_holdout_interp_mean_similarities_Spearman[a] <- mean(rep_holdout_interp_similarities_Spearman)

}
cat("\n[mean] global absolute average absolute percentage difference between iteration target and original target = ", dec_five(mean(rep_holdout_mean_differences_between_global_targets_and_fold_targets)), "% ± ", dec_five(sd(rep_holdout_mean_differences_between_global_targets_and_fold_targets)), "% ↑↑↑\n", sep="")
cat("[DTW] global DTW similarities between iteration target and original target = ", dec_five(mean(rep_holdout_dtw_mean_similarities)), " ± ", dec_five(sd(rep_holdout_dtw_mean_similarities)), " ↓↓↓\n", sep="")
cat("[Pearson rep_holdout_interpolation] global rep_holdout_interpolated similarities between iteration target and original target = ", dec_five(mean(rep_holdout_interp_mean_similarities_Pearson)), " ± ", dec_five(sd(rep_holdout_interp_mean_similarities_Pearson)), " ↑↑↑\n", sep="")
cat("[Kendal rep_holdout_interpolation] global rep_holdout_interpolated similarities between iteration target and original target = ", dec_five(mean(rep_holdout_interp_mean_similarities_Kendall)), " ± ", dec_five(sd(rep_holdout_interp_mean_similarities_Kendall)), " ↑↑↑\n", sep="")
cat("[Spearman rep_holdout_interpolation] global rep_holdout_interpolated similarities between iteration target and original target = ", dec_five(mean(rep_holdout_interp_mean_similarities_Spearman)), " ± ", dec_five(sd(rep_holdout_interp_mean_similarities_Spearman)), " ↑↑↑\n", sep="")

cat(num_iterations, "-times repeated hold-out validation\t Random Forests regression analysis\n", sep="")
cat(" =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  = \n")

k <- 10

crossval_mean_differences_between_global_targets_and_fold_targets <- c()

crossval_global_interations <- rep_holdout_global_interations
cat("crossval_global_interations = ", crossval_global_interations, "\n", sep="")

crossval_dtw_mean_similarities <- 0
crossval_interp_mean_similarities_Pearson <- 0
crossval_interp_mean_similarities_Kendall <- 0
crossval_interp_mean_similarities_Spearman <- 0


for(a in seq(1:crossval_global_interations)) {

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

    crossval_target_differences_perc_means <- numeric(k)

    folds <- sample(1:k, nrow(data), replace = TRUE)

    # Perform k-fold cross-validation
    for (i in 1:k) {
      # Split data into training and testing sets
      train_data <- data[folds != i, ]
      test_data <- data[folds == i, ]

      actuals <- test_data[[response_var]]

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


      crossval_target_differences_perc_means[i] <- thisDiff
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
#     average_r_squared_V3 <- mean(r_squared_values_V3)
#     sd_r_squared_V3 <- sd(r_squared_values_V3)
#     if(VERBOSE) cat("Average R-squared V3 over ",  k, " folds: ", average_r_squared_V3, " ± ", sd_r_squared_V3, "  ", sep="")
#     if(VERBOSE) cat(" in the [", dec_five(min(r_squared_values_V3)), ",", dec_five(max(r_squared_values_V3)), "] interval\n", sep="")

    # cat("(a=", a, ") absolute average absolute percentage difference between iteration target and original target = ", dec_five(mean(crossval_target_differences_perc_means)), "% ± ", dec_five(sd(crossval_target_differences_perc_means)), "%\n", sep="")
    # cat("(a=", a, ") mean DTW similarity between fold targets and global targets = ", dec_five(mean(crossval_dtw_similarities)), " ± ", dec_five(sd(crossval_dtw_similarities)), "\n", sep="")
    # cat("(a=", a, ") mean crossval_interpolation similarity between fold targets and global targets = ", dec_five(mean(crossval_interp_similarities_Pearson)), " ± ", dec_five(sd(crossval_interp_similarities_Pearson)), "\n", sep="")

    crossval_mean_differences_between_global_targets_and_fold_targets[a] <- mean(crossval_target_differences_perc_means)
    crossval_dtw_mean_similarities[a] <- mean(crossval_dtw_similarities)
    crossval_interp_mean_similarities_Pearson[a] <- mean(crossval_interp_similarities_Pearson)
    crossval_interp_mean_similarities_Kendall[a] <- mean(crossval_interp_similarities_Kendall)
    crossval_interp_mean_similarities_Spearman[a] <- mean(crossval_interp_similarities_Spearman)

}


cat("\n[mean] global absolute average absolute percentage difference between iteration target and original target = ", dec_five(mean(crossval_mean_differences_between_global_targets_and_fold_targets)), "% ± ", dec_five(sd(crossval_mean_differences_between_global_targets_and_fold_targets)), "% ↑↑↑\n", sep="")
cat("[DTW] global DTW similarities between iteration target and original target = ", dec_five(mean(crossval_dtw_mean_similarities)), " ± ", dec_five(sd(crossval_dtw_mean_similarities)), " ↓↓↓\n", sep="")
cat("[Pearson crossval interpolation] global crossval_interpolated similarities between iteration target and original target = ", dec_five(mean(crossval_interp_mean_similarities_Pearson)), " ± ", dec_five(sd(crossval_interp_mean_similarities_Pearson)), " ↑↑↑\n", sep="")
cat("[Kendal crossval interpolation] global crossval_interpolated similarities between iteration target and original target = ", dec_five(mean(crossval_interp_mean_similarities_Kendall)), " ± ", dec_five(sd(crossval_interp_mean_similarities_Kendall)), " ↑↑↑\n", sep="")
cat("[Spearman crossval interpolation] global crossval_interpolated similarities between iteration target and original target = ", dec_five(mean(crossval_interp_mean_similarities_Spearman)), " ± ", dec_five(sd(crossval_interp_mean_similarities_Spearman)), " ↑↑↑\n", sep="")

cat(k, "-fold cross-validation\t Random Forests regression analysis\n", sep="")
cat(" =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  = \n")


cat("\n[mean of the means] repeated holdout ", dec_five(mean(rep_holdout_mean_differences_between_global_targets_and_fold_targets)), "% > ", dec_five(mean(crossval_mean_differences_between_global_targets_and_fold_targets)), "%  cross validation? ↑↑↑ ", sep="")
if(mean(rep_holdout_mean_differences_between_global_targets_and_fold_targets) > mean(crossval_mean_differences_between_global_targets_and_fold_targets)) {
  cat("TRUE, repeated holdout wins ", green("\u2714"), " ", sep="")
} else cat("FALSE, cross-validation wins ", red("\u2716"), " ", sep="")
cat("diff = ", dec_two(computeDiffPerc(mean(rep_holdout_mean_differences_between_global_targets_and_fold_targets),mean(crossval_mean_differences_between_global_targets_and_fold_targets))), "%", sep="")

cat("\n[DTW comparison] repeated holdout ", dec_five(mean(rep_holdout_dtw_mean_similarities)), " < ", dec_five(mean(crossval_dtw_mean_similarities)), " cross validation? ↓↓↓ ", sep="")
if(mean(rep_holdout_dtw_mean_similarities) > mean(crossval_dtw_mean_similarities)) {
  cat("TRUE, repeated holdout wins ", green("\u2714"), " ", sep="")
} else cat("FALSE, cross-validation wins ", red("\u2716"), " ", sep="")
cat("diff = ", dec_two(computeDiffPerc(mean(rep_holdout_dtw_mean_similarities),mean(crossval_dtw_mean_similarities))), "%", sep="")

cat("\n[Pearson interpolation] repeated holdout ", dec_five(mean(rep_holdout_interp_mean_similarities_Pearson)), " > ", dec_five(mean(crossval_interp_mean_similarities_Pearson)), " cross validation? ↑↑↑ ", sep="")
if(mean(rep_holdout_interp_mean_similarities_Pearson) > mean(crossval_interp_mean_similarities_Pearson)) {
  cat("TRUE, repeated holdout wins ", green("\u2714"), " ", sep="")
} else cat("FALSE, cross-validation wins ", red("\u2716"), " ", sep="")
cat("diff = ", dec_two(computeDiffPerc(mean(rep_holdout_interp_mean_similarities_Pearson),mean(crossval_interp_mean_similarities_Pearson))), "%", sep="")

cat("\n[Kendall interpolation] repeated holdout ", dec_five(mean(rep_holdout_interp_mean_similarities_Kendall)), " > ", dec_five(mean(crossval_interp_mean_similarities_Kendall)), " cross validation? ↑↑↑ ", sep="")
if(mean(rep_holdout_interp_mean_similarities_Kendall) > mean(crossval_interp_mean_similarities_Kendall)) {
  cat("TRUE, repeated holdout wins ", green("\u2714"), " ", sep="")
} else cat("FALSE, cross-validation wins ", red("\u2716"), " ", sep="")
cat("diff = ", dec_two(computeDiffPerc(mean(rep_holdout_interp_mean_similarities_Kendall),mean(crossval_interp_mean_similarities_Kendall))), "%", sep="")

cat("\n[Spearman interpolation] repeated holdout ", dec_five(mean(rep_holdout_interp_mean_similarities_Spearman)), " > ", dec_five(mean(crossval_interp_mean_similarities_Spearman)), " cross validation? ↑↑↑ ", sep="")
if(mean(rep_holdout_interp_mean_similarities_Spearman) > mean(crossval_interp_mean_similarities_Spearman)) {
  cat("TRUE, repeated holdout wins ", green("\u2714"), " ", sep="")
} else cat("FALSE, cross-validation wins ", red("\u2716"), " ", sep="")
cat("diff = ", dec_two(computeDiffPerc(mean(rep_holdout_interp_mean_similarities_Spearman),mean(crossval_interp_mean_similarities_Spearman))), "%", sep="")

cat("\n")
cat("number of global interactions = ", rep_holdout_global_interations, "\n", sep="")
cat("dataFile: ", dataFile, "\n", sep="")

cat("\n")
computeExecutionTime()

cat("~ : ~ : ~ : ~ : ~ The end ~ : ~ : ~ : ~ : ~\n")
