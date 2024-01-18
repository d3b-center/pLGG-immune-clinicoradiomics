
### Authors: Komal S. Rathi, Adam Kraya

### Purpose

Given a cox model for prediction of risk of progression in LGG; the module identifies what pathways are associated with the predicted risk.

### Data

[Release 20230826](https://cavatica.sbgenomics.com/u/d3b-bixu-ops/monthly-release-data/files/#q?path=20230826_release) for annotated histologies and TPM data.

201 unique cohort participant identifiers in pLGG risk group file. 

### 1. ssGSEA

`01-ssgsea.R`:  161 biospecimen identifiers corresponding to 201 cohort participant identifiers in the risk group file were pulled from the `20230826_release` histology file. GSVA using `ssgsea` method was performed on `REACTOME` gene sets with minimum set size of 10 and maximum set size of 500.

#### Run script

```
Rscript --vanilla 01-ssgsea.R --help

Usage: 01-ssgsea.R [options]

Options:
	--mat=MAT
		input RNA-seq gene expression matrix (preferably TPM)

	--histology_file=HISTOLOGY_FILE
		histology file for all OpenPedCan samples (.tsv)

	--risk_file=RISK_FILE
		file with risk scores

	--gtf_file=GTF_FILE
		gencode GTF file

	--output_dir=OUTPUT_DIR
		output directory

	-h, --help
		Show this help message and exit
```

#### Output

Matrix of 1292 gene sets across 161 biospecimens was saved under `ssgsea_matrix.rds`.

```
results
└── ssgsea_matrix.rds
```

### 2. Elastic NET regression

`02-elastic_net_risk.R`: Run elastic net regression on pathways obtained in step 1.

#### Run script

```
Rscript --vanilla 02-elastic_net_risk.R --help

Usage: 02-elastic_net_risk.R [options]

Options:
	--histology_file=HISTOLOGY_FILE
		histology file for all OpenPedCan samples (.tsv)

	--risk_file=RISK_FILE
		file with risk scores

	--output_dir=OUTPUT_DIR
		output directory

	--plots_dir=PLOTS_DIR
		plots directory

	-h, --help
		Show this help message and exit
```

#### Output

Coefficients and outputs:

```
results
├── EL_testing_prediction_summary.txt
├── EL_training_prediction_summary.txt
├── estimated_testing_error_cv.txt
└── full_coefficient_table.txt
```

Associated plots:

```
plots
├── EL_metrics.pdf
└── coef_EL_LGGrisk.pdf
```
