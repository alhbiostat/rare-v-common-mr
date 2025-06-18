#!/bin/bash
#To run: source scripts/02_dataextraction_proteins/runformatting_proteins.sh

source "config.env"

#studies=${project_dir}/allstudypairings_proteins.csv
studies=${project_dir}/allstudypairings_proteins_sunout.csv
scripts_dir=${project_dir}/scripts/02_dataextraction_proteins
out_dir=${data_dir}/harmonised/proteomics
log_dir=${out_dir}/logs

mkdir -p $log_dir

cd ${project_dir}

# Loop over exposure and outcome study pairs and format according to study source and instrument class

tail -n +2 $studies | while IFS=, read -r exposure_study outcome_study exposure_name outcome_name exposure_source outcome_source class
do
  log_file="${log_dir}/${exposure_source}_${class}_${exposure_name}_${outcome_name}.log"
  out_file="${out_dir}/${exposure_source}_${class}_${exposure_name}_${outcome_name}.rda"

  if [[ -f $out_file ]]; then
    echo "File $out_file already exists. Skipping..."
    continue
  fi
  
  if [[ $exposure_source == "sun_gwas" ]]; then

    echo "Formatting studies...
    Exposure: $exposure_name ($exposure_study)
    Outcome: $outcome_name ($outcome_study)
    Source: $exposure_source
    Class: $class" | tee -a $log_file

    Rscript ${scripts_dir}/003_formatcommonstudies.R $exposure_study $outcome_study $exposure_name $outcome_name $exposure_source $outcome_source $class >> $log_file 2>&1 

    if [ $? -ne 0 ]; then
          echo "Error: Harmonisation failed for $exposure_name and $outcome_name. Skipping..."
          continue
    fi

  elif [[ $exposure_source == "dhindsa_exwas" ]]; then

    echo "Formatting studies...
    Exposure: $exposure_name ($exposure_study)
    Outcome: $outcome_name ($outcome_study)
    Source: $exposure_source
    Class: $class" | tee -a $log_file

    Rscript ${scripts_dir}/004_formatrarestudies.R $exposure_study $outcome_study $exposure_name $outcome_name $exposure_source $outcome_source $class >> $log_file 2>&1 

    if [ $? -ne 0 ]; then
          echo "Error: Harmonisation failed for $exposure_name and $outcome_name. Skipping..."
          continue
    fi
  fi
done