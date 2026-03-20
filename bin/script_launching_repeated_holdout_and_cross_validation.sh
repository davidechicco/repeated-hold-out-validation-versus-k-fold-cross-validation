#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#
set -o nounset -o pipefail -o errexit
# set -o xtrace

current_date_and_time=$(date +"%Y-%m-%d %H:%M:%S" | sed 's/ /_h/; s/:/m/; s/:/s/')
date=$(date +"%Y-%m-%d")

# Define the list as an array
datasets=("diabetes_type_one" "heart_failure" "sepsis" "chronic_kidney_disease" "obesity")
number_of_global_interations="100 200 300 500 1000 2000 5000"
# number_of_global_interations="1000 2000 5000"
length_of_number_of_global_interations=($number_of_global_interations)

# datasets=($(printf "%s\n" "${datasets_list[@]}" | sort))

k_list="3 5 7 10 12 20"

analysis_type="binary"
# analysis_type="regression"

outputFolder="../results"
mkdir -p $outputFolder

test_number=1
length_datasets=${#datasets[@]}
length_arr=${#length_of_number_of_global_interations[@]}
length_k_list_temp=($k_list)
length_k_list=${#length_k_list_temp[@]}
number_of_tests_foreseen=$(( length_datasets * length_arr * length_k_list))
echo -e "We're going to do "$number_of_tests_foreseen" tests, fasten your seatbelts"

# Loop through the array and print each condition
for thisK in $k_list; do
  for dataset in "${datasets[@]}"; do
    for global_iteration in $number_of_global_interations; do

    echo -n "[test "$test_number"/"$number_of_tests_foreseen"]"
    echo -n " iterations="$global_iteration
    echo -n " dataset chosen: "$dataset
    echo -n " k="$thisK" "

    outputFile=$outputFolder"/output_"$analysis_type"_"$thisK"k"_$dataset"_iterations"$global_iteration"_"$current_date_and_time".txt"
    Rscript repeated_holdout_and_cross_validation_regression.r $global_iteration $dataset $thisK $analysis_type >> $outputFile 2>> $outputFile
    echo -e $outputFile
    test_number=$((test_number+1))
     done
  done
done



