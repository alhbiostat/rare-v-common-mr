#!/bin/bash

while IFS= read -r line; do 
    filename=$(echo $line | sed 's/ /_/g')
    awk 'NR==1' UKBB_470k_deeprvat_results.csv > deeprvat_${filename}.csv
    grep -F "$line" UKBB_470k_deeprvat_results.csv >> deeprvat_${filename}.csv 
done < deeprvat_phenotypes.txt