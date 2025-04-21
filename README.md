<h1 align="center">Euclidean Distance Method</h1>
<h4 align="center">Biological distance between tumors and their lab grown models</h4>


## Table of contents
- [Quickstart Guide](#quickstart-guide)
- [Download Data from Manifest File Using the GDC Client](#download-data-from-manifest-file-using-the-gdc-client)
- [Run Processing Pipeline](#run-processing-pipeline)
- [Calculate Euclidean Distances](#calculate-euclidean-distances)



## Quickstart Guide

### Setup

Install requirements - detailed instructions are found on the [Requirements page](doc/requirements.md):

1. Install Python 3+
2. Install GDC Data Transfer Tool Client

Ensure that steps are completed on the [Requirements page](doc/requirements.md) - *(includes creating working environment, fetching classifier submodule, and manually downloading required data)*

## Download Data from Manifest File Using the GDC Client
Download gene expression data. 

Specify a single cancer cohort (see bullets below on options, example `PAAD` for pancreatic cancer cohort) or use `ALL` to run all cancer cohorts.

```bash
bash scripts/gdc_download.sh PAAD
```

This will create subfolders in `data-raw/<CANCER>_GEXP_<TYPE>` and place GDC molecular matrices here.

> Options for cancer cohort includes `ALL`, `BLCA`, `BRCA`, `COADREAD`, `ESO`, `HNSC`, `KID`, `LGGGBM`, `LIHCCHOL`, `LUNG`, `OV`, `PAAD`, `SARC`, `SKCM`, `UCEC`

For more details on each cancer cohort option see [Cohort Options Page](doc/cohort_options.md)


## Pre-process Data
Run pre-processing pipeline.

Specify a single cancer cohort (see bullets below on options, example `PAAD` for pancreatic cancer cohort) or use `ALL` to run all cancer cohorts.
```bash
bash scripts/process.sh PAAD data/prep
```

> Creates file `data/prep/<CANCER>_GEXP/<CANCER>_GEXP_prep2_<TYPE>.tsv` that is prepped for distance calculations

> Options for cancer cohort includes `ALL`, `BLCA`, `BRCA`, `COADREAD`, `ESO`, `HNSC`, `KID`, `LGGGBM`, `LIHCCHOL`, `LUNG`, `OV`, `PAAD`, `SARC`, `SKCM`, `UCEC`

For more details on each cancer cohort option see [Cohort Options Page](doc/cohort_options.md)


## Determine Biological Distances Between Tumor and Derived Model Pairs
Calculate euclidean distances.

Specify a single cancer cohort (see bullets below on options, example `PAAD` for pancreatic cancer cohort) or use `ALL` to run all cancer cohorts.
```bash
bash scripts/euc_distance.sh PAAD data/prep
```

> Options for cancer cohort includes `ALL`, `BLCA`, `BRCA`, `COADREAD`, `ESO`, `HNSC`, `KID`, `LGGGBM`, `LIHCCHOL`, `LUNG`, `OV`, `PAAD`, `SARC`, `SKCM`, `UCEC`

For more details on each cancer cohort option see [Cohort Options Page](doc/cohort_options.md)


Final results are found in `data/distance_metric/main_results/`: 

+ Euclidean distances and z-scores: `euclidean_distances_HCMITumor.Model_PAAD.tsv`
+ Outlier information: `euclidean_distances_outlier_samples_HCMITumor.Model_PAAD.tsv`
