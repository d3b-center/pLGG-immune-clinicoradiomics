### Author: Komal Rathi

### Purpose

We have built a Cox model for prediction of risk of progression in LGG; we would like to know what pathways are associated with the predicted risk.

### Data

[Release 20230826](https://cavatica.sbgenomics.com/u/d3b-bixu-ops/monthly-release-data/files/#q?path=20230826_release) for annotated histologies and TPM data.

201 unique cohort participant identifiers in Anahita's risk group file (LGG). 

### Methods

161 biospecimen identifiers corresponding to 201 cohort participant identifiers in the risk group file were pulled from the `20230826_release` histology file. GSVA using `ssgsea` method was performed on `REACTOME` gene sets with minimum set size of 10 and maximum set size of 500.

### Results

Matrix of 1292 gene sets across 161 biospecimens was saved under `ssgsea_matrix.rds`.
