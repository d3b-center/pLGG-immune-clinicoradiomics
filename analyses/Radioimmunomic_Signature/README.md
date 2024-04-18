## MATLAB codes for radioimmunomic signature
```
Author: Anahita Fathi Kazerooni

```

### Software Requirements
- MATLAB 2023A
  - Parallel Computing Toolbox
  - Statistics and Machine Learning Toolbox

### Hardware Requirements
multicore processors, GPUs, or computer clusters


### Code and Data Description
- `dnafinder-roc-b84f535.zip` : ROC analysis package. This file needs to be unzipped before running the scripts
- `ModelTraining.m` : Script for ML model training for radioimmunomic signature
  - Input: `Features_ModelTraining.mat`
- `PredictImmuneClass.m`
  - Inputs: `Features_ModelTraining.mat`, `SelectedModel.mat`


