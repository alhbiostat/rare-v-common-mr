#!/bin/bash

source ../../config.env

phenocodes="${data_dir}/sumstats/deeprvat/deeprvat_ukbdatafields_continuous.txt"
binary="${data_dir}/sumstats/deeprvat/deeprvat_binary.txt"
genebass_studies="${data_dir}/sumstats/genebass/genebass_wes_studymetadata.tsv"

# Continuous traits
awk -F '\t' 'NR==FNR {vals[$1]; next} $7 in vals' ${phenocodes} ${genebass_studies} > ${data_dir}/sumstats/genebass/genebass_wes_requiredstudies.tsv

# Binary traits (requires subsequent manual filtering)
awk -F '\t' 'NR==FNR {vals[$1]; next} {for (i in vals) if (tolower($15) ~ tolower(i)) print $0}' ${binary} ${genebass_studies} >> ${data_dir}/sumstats/genebass/genebass_wes_requiredstudies.tsv