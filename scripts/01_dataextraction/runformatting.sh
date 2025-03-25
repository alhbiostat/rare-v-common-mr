#!/bin/bash

source "../../config.env"

#studies=${project_dir}/studypairings.csv
studies=${project_dir}/test.csv
scripts_dir=${project_dir}/scripts/01_dataextraction
out_dir=${data_dir}/harmonised
log_dir=${out_dir}/logs

mkdir -p $log_dir

cd ${project_dir}

# Loop over exposure and outcome study pairs and format according to study source and instrument class

tail -n +2 $studies | while IFS=, read -r exposure_study outcome_study exposure_name outcome_name source class
do
  exposure_study=$(echo $exposure_study | tr -d '"')
  outcome_study=$(echo $outcome_study | tr -d '"')

  log_file="${log_dir}/${source}_${class}_${exposure_name}_${outcome_name}.log"
  out_file="${out_dir}/${source}_${class}_${exposure_name}_${outcome_name}.rda"

  if [[ $source == "genebass" && ! -f $out_file ]]; then

    echo "Formatting studies...
    Exposure: $exposure_name ($exposure_study)
    Outcome: $outcome_name ($outcome_study)
    Source: $source
    Class: $class" | tee -a $log_file

    Rscript ${scripts_dir}/003_formatgenebassstudies.R $exposure_study $outcome_study $exposure_name $outcome_name $source $class >> $log_file 2>&1 

    if [ $? -ne 0 ]; then
          echo "Error: Harmonisation failed for $exposure_name and $outcome_name. Skipping..."
          continue
    fi

  else
    echo "File $out_file already exists. Skipping..."
    continue
  fi
done