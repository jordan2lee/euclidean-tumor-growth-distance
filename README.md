# tumor-growth-distance
Biological distance between tumors and their lab grown models

# set up

```bash
python3 -m venv venv; source venv/bin/activate
pip install -r requirements.txt 
```

# files not included in repo but are required are:
data/distance_metric/src/missing_ohsu_euclidean_distance_all.02.12.2024.csv

# Calculate distances and visualize
Set up
```
mkdir -p data/distance_metric/${cancer}_GEXP
```

Run
```
scripts/dist_model-tumor.ipynb
```
> creates file data/distance_metric/main_results/important/mdelta.tsv that contains all mdelta values (match and mismatchs, within same cohort and across cohorts)
> Requres input files, including info on which samples are correct matches `src/model_tumor_pairings_Pulled_analysis_tracker_OHSU_tab_on_05.03.24.txt`


If want to get pvalues, continue with (WIP)
```
scripts/fdr_correction_only.Rmd

```
then can generate plots of the distances here
```
scripts/visualize_distances.ipynb
violin plots visualize_califano.ipynb
```

### Notes
The source file `data/distance_metric/src/model-tumor-relations.json` was created by running this code `python scripts/create_model-tumor-relations.py` and is used as a look up dictionary for samples when taking other methods distance file (all vs all samples) and converting to p-values.
