# Code Developed for Clinicoradiomic Risk Stratification in pLGG
## By: Meen Chul Kim and Anahita Fathi Kazerooni

This code inputs the files including normalized radiomic and clinical variables with data splits (into Discovery/Replication sets), and the progression-free survival and censoring information:
- 02_cox_radiomic.ipynb: Clinicoradiomic Risk stratification
- Inputs:
  - Clinical Variables: X_clinical.pkl
  - Radiomic Features: X_radiomic.pkl
  - Progression-Free Survival and Censoring Information: y.pkl
- Ouputs:
  - Risk scores and groups for the Disocvery set: st_radiomic_discovery.csv
  - Risk scores and groups for the Replication set: rs_radiomic_replicate.csv
