#!/usr/bin/bash

FAIL="\033[0;31m"
PASS="\033[0;32m"
NC="\033[0m"

cpus=4

# metabolome search parameters
min_met=170
max_met=250
met_steps=10
# microbiome search parameters
min_mic=70
max_mic=150
mic_steps=10

# paths
mic_data_file="Data/16S_data.RData"
base_results="Results"

echo "Starting analysis..."

err_file="$base_results/model_fitting.err"
if Rscript multi_omics.R "$mic_data_file" "$base_results" \
 --min-metabolome "$min_met" --max-metabolome "$max_met" \
 --metabolome-step-size "$met_steps" \
 --min-microbiome "$min_mic" --max-microbiome "$max_mic" \
 --microbiome-step-size "$mic_steps" \
 --cpus=$cpus \
 2> "$err_file"; then
	echo -e "${PASS}Finished model fitting!${NC}\n=====\n"
else
	echo -e "${FAIL}Model fitting failed. Logs can be found in $err_file. ${NC}\n=====\n"
fi
