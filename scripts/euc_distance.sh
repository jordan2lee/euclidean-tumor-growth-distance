#!/usr/bin/bash
set -e
cancer=${1}
prep_path=${2}
# Run all cohorts
if [[ ${cancer} == 'ALL' ]]; then
  declare -a StringArray=('BLCA' 'BRCA' 'COADREAD' 'ESO' 'KIRCKICH' 'LGGGBM' 'LIHCCHOL' 'LUNG' 'OV' 'PAAD' 'SARC' 'SKCM' 'UCEC')
  for cancer in ${StringArray[@]}; do
    mkdir -p data/distance_metric/${cancer}_GEXP
    python scripts/dist_model-tumor.py --cancer PAAD --prep_path ${prep_path}
    done
# Run single cohort
else
  mkdir -p data/distance_metric/${cancer}_GEXP
  python scripts/dist_model-tumor.py --cancer PAAD --prep_path ${prep_path}
fi