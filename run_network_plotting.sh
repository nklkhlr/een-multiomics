#!/usr/bin/bash


FAIL="\033[0;31m"
PASS="\033[0;32m"
NC="\033[0m"

function check_status() {
  if [ $1 -eq 0 ]; then
    echo -e  "${PASS}Successfully finished $2!${NC}"
  else
    echo -e  "${FAIL}Failed $2! Please check logs in $3.${NC}"
  fi
}

base_results="Results"
plot_results="community_network_plotting"


echo "Starting plotting of network..."
plot_path="$plot_results"

if [ ! -d "$plot_path" ]; then
  mkdir "$plot_path"
fi

Rscript community_network_plotting/community_network_data.R \
  "$base_results" "$plot_path" \
  2> "$plot_path/data_extraction.err"
status=$?
check_status $status "data extraction" "$plot_path/data_extraction.err"

if [ $status -eq 0 ]; then
	python community_network_plotting/community_network_plotting.py \
	  --data-path "$plot_path" 2> "$plot_path/network_plotting.err"
	status=$?
	check_status $status "network plotting" "$plot_path/network_plotting.err"
fi

echo -e "=====\n"
