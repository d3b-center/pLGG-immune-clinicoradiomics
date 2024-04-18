# Code Developed for Clinicoradiomic Risk Stratification in pLGG

```
Developers:

- Meen Chul Kim, PhD (kimm11@chop.edu)
- Anahita Fathi Kazerooni, PhD (fathikazea@chop.edu)
```

This repository contains code that processes data files divided into Discovery and Replication sets, with the following features:

- Input features
  - Normalized clinical variables (`X_clinical.pkl`):
    - Sex (categorical)
    - Age at Diagnosis (continuous)
    - Tumor Location (categorical)
    - NF1 status (categorical)
    - Extent of Tumor Resection (continuous)
    - Chemotherapy (categorical)
    - Radiation treatment (categorical)
  - Normalized radiomic variables (`X_radiomic.pkl`): A total of 439 variables selected through statistical filter methods.

- Output features (`y.pkl`):
  - Progression-free survival (continuous)
  - Censoring status (boolean)

- Clinicoradiomic risk stratification pipeline (`Risk Stratification.ipynb`)

Running the risk stratification pipeline produces the following outputs:

- For the Discovery set:
  - Risk scores and categories (`rs_discovery.csv`)
  - Kaplan-Meier survival curves (`km_discovery.png`)
  - A table of survival data (`st_discovery.csv`)
  - A forest plot visualizing the impact of variables (`fp_discovery.png`)

- For the Replication set:
  - Risk scores and categories (`rs_replicate.csv`)
  - Kaplan-Meier survival curves (`km_replicate.png`)
  - A table of survival data (`st_replicate.csv`)
