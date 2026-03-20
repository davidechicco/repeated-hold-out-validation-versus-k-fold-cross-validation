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

if (length(args) < 4) {
  stop("Please provide four arguments")
}

rep_holdout_global_interations <- as.numeric(args[1])
dataset_name <- args[2]
k <- as.numeric(args[3])
type <- toString(args[4])


cat("argument 1, rep_holdout_global_interations:", rep_holdout_global_interations, "\n")
cat("argument 2, dataset:", dataset_name, "\n")
cat("argument 3, k:", k, "\n")
cat("argument 4, type:", type, "\n")

BINARY_CLASS <- NULL
REGRESSION <- NULL

if(type=="binary") {
  BINARY_CLASS <- TRUE
  REGRESSION <- !(BINARY_CLASS)
} else if(type=="regression") {
  REGRESSION <- TRUE
  BINARY_CLASS <- !(REGRESSION)
}

# Load the data
if(dataset_name == "chronic_kidney_disease") {
  dataFile <- "../data/CKD_CVD_journal.pone.0199920.s002_EDITED_IMPUTED_v2.csv"
  if(REGRESSION) response_var <- "TimeToEventMonths"
  if(BINARY_CLASS) response_var <- "EventCKD35"
} else if(dataset_name == "sepsis") {
  dataFile <- "../data/sepsis_severity_dataset_col_norm_edited_target-SOFA-score.csv"
  if(REGRESSION) response_var <- "SOFA.score"  # Replace with your response variable name
  if(BINARY_CLASS) response_var <- "cancer"
} else if(dataset_name == "heart_failure") {
  dataFile <- "../data/EHRs_heart_failure_S1Data_EDITED_v3.csv"
  if(REGRESSION) response_var <- "TIME_DAYS"  # Replace with your response variable name
  if(BINARY_CLASS) response_var <- "DEATH_EVENT"
} else if(dataset_name == "obesity") {
  dataFile <- "../data/ObesityDataSet_raw_and_data_sinthetic_EDITED.csv"
  if(REGRESSION) response_var <- "obesity_class"  # Replace with your response variable name
  if(BINARY_CLASS) response_var <- "family_history_with_overweight"
} else if(dataset_name == "diabetes_type_one") {
  dataFile <- "../data/journal.pone.0216416_Takashi2019_diabetes_type1_dataset_preprocessed.csv"
  if(REGRESSION) response_var <- "duration.of.diabetes"
  if(BINARY_CLASS) response_var <- "insulin_regimen_binary"
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

cat("predictor_vars:\n")
#cat("response_var: ", data[, c(response_var)], "\n", sep="")

print((data[, c(response_var)] %>% table())*100/nrow(data))


# Set parameters for repeated hold-out validation
num_iterations <- k
cat("here p = k = ", k, "\n", sep="")

train_fraction <- 1-(1/k)

cat("training set = ", dec_five(100*train_fraction), "%\n", sep="")
rep_holdout_target_differences_perc_means <- numeric(num_iterations)

cat(num_iterations, "-times repeated hold-out validation\t ", sep="")
if(BINARY_CLASS) cat(" binary classification\n")
if(REGRESSION) cat(" regression analysis\n")

cat("global target mean ", mean(data[[response_var]]),  " ± ", sd(data[[response_var]]), "\n", sep="")

rep_holdout_mean_differences_between_global_targets_and_fold_targets <- c()
rep_holdout_dtw_mean_similarities <- 0
rep_holdout_interp_mean_similarities_Pearson <- 0
rep_holdout_interp_mean_similarities_Kendall <- 0
rep_holdout_interp_mean_similarities_Spearman <- 0

cat("rep_holdout_global_interations = ", rep_holdout_global_interations, "\n", sep="")

countTimesRepetedHoldoutBetterThanKfold <- 0

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

      rep_holdout_dtw_similarities[i] <- rep_holdout_dtw_similarity
    }

        rep_holdout_dtw_mean_similarities[a] <- mean(rep_holdout_dtw_similarities)
}

cat("\n")
cat("[DTW] global DTW similarities between iteration target and original target = ", dec_five(mean(rep_holdout_dtw_mean_similarities)), " ± ", dec_five(sd(rep_holdout_dtw_mean_similarities)), " ↑↑↑\n", sep="")
cat(num_iterations, "-times repeated hold-out validation\t ", sep="")
cat("\n =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  = \n")


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

    crossval_target_differences_perc_means <- numeric(k)

    folds <- sample(rep(1:k, length.out = nrow(data)))

    # Perform k-fold cross-validation
    for (i in 1:k) {
      # Split data into training and testing sets
      train_data <- data[folds != i, ]
      test_data <- data[folds == i, ]

      actuals <- test_data[[response_var]]

      if(VERBOSE) cat("[k=", k, ",a=", a, ",i=", i, " of ", k, "] ", sep="")

      if(VERBOSE) cat("#rows of this training set = ", nrow(train_data), " ", sep="")
      if(VERBOSE) cat("#rows of this test set = ", nrow(test_data), "\n", sep="")
      if(VERBOSE) cat("target mean ", mean(actuals),  " ± ", sd(actuals), "\n", sep="")

      if(VERBOSE) cat(" mean diff with original = ")
      thisDiff <- abs(computeDiffPerc(mean(data[[response_var]]),mean(actuals)))
      if(VERBOSE) cat(dec_five(thisDiff), "%\n", sep="")

      crossval_dtw_alignment <- dtw(data[[response_var]], actuals)
      crossval_dtw_similarity <- 1 / (1 + crossval_dtw_alignment$"distance")  # Convert distance to similarity
      # DTW similarity (common transforms): values between 0 (very dissimilar) and 1 (identical)

         crossval_dtw_similarities[i] <- crossval_dtw_similarity
    }

    crossval_dtw_mean_similarities[a] <- mean(crossval_dtw_similarities)

    if(VERBOSE) cat("\n[corresponding iteration DTW comparison] repeated holdout ", dec_five(rep_holdout_dtw_mean_similarities[a]), " > ", dec_five(crossval_dtw_mean_similarities[a]), " cross validation? ↑↑↑ ", sep="")
    if(rep_holdout_dtw_mean_similarities[a] > crossval_dtw_mean_similarities[a]) {

    if(VERBOSE)  cat("TRUE, repeated holdout wins ", green("\u2714"), " ", sep="")
      countTimesRepetedHoldoutBetterThanKfold <- countTimesRepetedHoldoutBetterThanKfold + 1

    } else {

      if(VERBOSE) cat("FALSE, cross-validation wins ", red("\u2716"), " ", sep="")

    }

     if(VERBOSE) cat("diff = ", dec_two(computeDiffPerc(mean(rep_holdout_dtw_mean_similarities),crossval_dtw_mean_similarities[a])), "%", sep="")

}

cat("\n")

cat("[DTW] global DTW similarities between iteration target and original target = ", dec_five(mean(crossval_dtw_mean_similarities)), " ± ", dec_five(sd(crossval_dtw_mean_similarities)), " ↑↑↑\n", sep="")

cat(k, "-fold cross-validation\t ", sep="")
if(BINARY_CLASS) cat(" binary classification\n")
if(REGRESSION) cat(" regression analysis\n")
cat(" =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  = \n")


percTimesRepetedHoldoutBetterThanKfold <- (countTimesRepetedHoldoutBetterThanKfold * 100) / crossval_global_interations
cat("\n[DTW] number of times repeated holdout was better than kfold: ", countTimesRepetedHoldoutBetterThanKfold, "/", crossval_global_interations,  " = ", dec_two(percTimesRepetedHoldoutBetterThanKfold),"% ", sep="")

percMajorityTimes <- 50
if(percTimesRepetedHoldoutBetterThanKfold > percMajorityTimes) {
  cat("TRUE, repeated holdout wins in most of the times ", green("\u2714"), " ", sep="")
} else {
  cat("FALSE, cross-validation wins in most of the times ", red("\u2716"), " ", sep="")
  }

cat("\n")
cat("[mean DTW comparison] repeated holdout ", dec_five(mean(rep_holdout_dtw_mean_similarities)), " > ", dec_five(mean(crossval_dtw_mean_similarities)), " cross validation? ↑↑↑ ", sep="")
if(mean(rep_holdout_dtw_mean_similarities) > mean(crossval_dtw_mean_similarities)) {
  cat("TRUE, repeated holdout wins ", green("\u2714"), " ", sep="")
} else cat("FALSE, cross-validation wins ", red("\u2716"), " ", sep="")
cat("diff = ", dec_two(computeDiffPerc(mean(rep_holdout_dtw_mean_similarities),mean(crossval_dtw_mean_similarities))), "%", sep="")

cat("\n")
cat("number of global interactions = ", rep_holdout_global_interations, "\n", sep="")
cat("dataFile: ", dataFile, "\n", sep="")

cat("\n")
computeExecutionTime()

cat("~ : ~ : ~ : ~ : ~ The end ~ : ~ : ~ : ~ : ~\n")
