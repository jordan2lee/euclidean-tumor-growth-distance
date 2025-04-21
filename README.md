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

> *Note on the selection of Cancer Type: These HCMI selections of Cancer Type were made and grouped together for each corresponding TCGA cancer cohort:*
>
> + 'BLCA':['Bladder Cancer']
> + 'BRCA':['Breast Cancer','Intraductal Carcinoma']
> + 'COADREAD': ['Colorectal Cancer', 'Rare Cancer','Serrated Adenoma','Tubular Adenoma','Tubulovilluos Adenoma']
> + 'ESO': ['Esophageal Cancer','Stomach Cancer']
> + 'HNSC':['Head and Neck']
> + 'KID':['Kidney Cancer']
> + 'LGGGBM':['Glioblastoma']
> + 'LIHCCHOL':['Extrahepatic Cholangiocarcinoma','Ampulla of Vater','Hepatocellular carcinoma','Intrahepatic Cholangiocarcinoma']
> + 'LUNG':['Lung Cancer']
> + 'OV':['Ovarian Cancer']
> + 'PAAD': ['Pancreatic Cancer']
> + 'SARC': ['Bone Cancer','Clear Cell Sarcoma','Desmoid Tumors','Epithelial Sarcoma','Ewing Sarcoma','Intimal Sarcoma','Leiomyosarcoma','Rhabdomyosarcoma','Spindle Cell Sarcoma','Undifferentiated Pleomorphic Sarcoma']
> + 'SKCM':['Melanoma']
> + 'UCEC':['Endometrial Cancer'],


## Pre-process Data
Run pre-processing pipeline.

Specify a single cancer cohort (see bullets below on options, example `PAAD` for pancreatic cancer cohort) or use `ALL` to run all cancer cohorts.
```bash
bash scripts/process.sh PAAD data/prep
```

> Creates file `data/prep/<CANCER>_GEXP/<CANCER>_GEXP_prep2_<TYPE>.tsv` that is prepped for distance calculations

## Determine Biological Distances Between Tumor and Derived Model Pairs
Calculate euclidean distances.

Specify a single cancer cohort (see bullets below on options, example `PAAD` for pancreatic cancer cohort) or use `ALL` to run all cancer cohorts.
```bash
bash scripts/euc_distance.sh PAAD data/prep
```

Final results are found in `data/distance_metric/main_results/`: 

+ Euclidean distances and z-scores: `euclidean_distances_HCMITumor.Model_PAAD.tsv`
+ Outlier information: `euclidean_distances_outlier_samples_HCMITumor.Model_PAAD.tsv`
