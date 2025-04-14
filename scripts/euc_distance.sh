#!/usr/bin/bash
set -e
cancer=${1}
prep_path=${2}

mkdir -p data/distance_metric/${cancer}_GEXP
python scripts/dist_model-tumor.py --cancer PAAD --prep_path ${prep_path}
