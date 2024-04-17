# pLGG-Immune-Clinicoradiomics
This code repository includes the data and source codes used in the manuscript "Multiparametric MRI Along with Machine Learning Informs on Molecular Underpinnings, Prognosis, and Treatment Response In Pediatric Low-Grade Glioma"

## Software Requirements
- CaPTk, v1.8.1 (https://cbica.github.io/CaPTk/)
- Python3 
- R v4.3
- MATLAB 2023A (v23.2)


## MRI Pre-processing and Tumor Segmentation:
- Required MRI sequences: T1, T1CE, T2, FLAIR (ADC optional)
- Pre-processing using BraTS Pre-processing Pipeline, details explained in: https://cbica.github.io/CaPTk/preprocessing_brats.html
- Tumor Segmentation, all details provided in: https://github.com/d3b-center/peds-brain-auto-seg-public
- Skull-stripping to generate a brain mask: https://github.com/d3b-center/peds-brain-auto-skull-strip
- Image normalization:
- Whole tumor generation:
- Radiomic feature extraction, using CaPTk v1.8.1: https://cbica.github.io/CaPTk/ht_FeatureExtraction.html; parameter file for radiomic feature extraction:


## Immune Profiling:
- analyses/lgg_xcell_analyses

## Radioimmunomic Analysis:
- analyses/Radioimmunomic_Signature

## Clinicoradiomic Risk Stratification:
- analyses/Clinicoradiomics

## Assessment of transcriptomic pathways associated with clinicoradiomic risk:
- analyses/lgg_risk_analysis
