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
```bash
python3 scripts/dist_model-tumor.py
```