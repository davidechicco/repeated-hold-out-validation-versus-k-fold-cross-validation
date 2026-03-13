#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#
set -o nounset -o pipefail -o errexit
set -o xtrace

current_date_and_time=$(date +"%Y-%m-%d %H:%M:%S" | sed 's/ /_h/; s/:/m/; s/:/s/')
date=$(date +"%Y-%m-%d")

# Define the list as an array
datasets=("chronic_kidney_disease" "diabetes_type_one" "heart_failure" "obesity" "sepsis")
number_of_global_interations="100 200 300 500 1000 2000 5000"
#number_of_global_interations="5000"

thisK=3

outputFolder="../results"
mkdir -p $outputFolder

# Loop through the array and print each condition
for dataset in "${datasets[@]}"; do
  echo "$dataset"
  for global_iteration in $number_of_global_interations; do
    echo "$global_iteration"
    outputFile=$outputFolder"/output_"$thisK"k"_$dataset"_iterations"$global_iteration"_"$current_date_and_time".txt"
    Rscript repeated_holdout_and_cross_validation_regression.r $global_iteration $dataset $thisK >> $outputFile 2>> $outputFile
  done
done



