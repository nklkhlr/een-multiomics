#!/usr/bin/bash

FAIL="\033[0;31m"
PASS="\033[0;32m"
NC="\033[0m"

# taxonomic levels to cover
base_results="Results"

x_thresh=1.0
y_thresh=0.7
greater="EEN"
dim="x"


echo "Starting downstream analysis..."

err_file="$base_results/downstream.err"
if Rscript downstream_analysis.R "$base_results" \
	--x-threshold "$x_thresh" --y-threshold "$y_thresh" \
	--association-greater "$greater" --association-dim "$dim" \
	  2> "$err_file"; then
	echo -e "${PASS}Finished downstream analysis!${NC}\n=====\n"
else
	echo -e "${FAIL}Downstream analysis failed. Please see logs in $err_file. ${NC}\n=====\n"
fi

