clear; clc; close all;
addpath('dnafinder-roc-b84f535');

%% load the data

load('Features_ModelTraining.mat')
%% Discovery vs Replication
load('SelectedModel.mat')

%
All_Features_discovery2 = All_Features_discovery(:,SelectedModel.selected_ranking);
All_Features_replicate2 = All_Features_replicate(:,SelectedModel.selected_ranking);



SVMModel = SelectedModel.SVMModel;
ScoreSVMModel = fitPosterior(SVMModel,All_Features_discovery2,class_discovery);

[~,score_rep] = predict(ScoreSVMModel,All_Features_replicate2);


sscore_rep = score_rep(:,2);


ROCout_rep = roc_dnafinder([sscore_rep,class_replicate],0,0.05,0);
auc_rep = ROCout_rep.AUC;


acc_rep = (1-ROCout_rep.xr + ROCout_rep.yr)/2;
[i1,i2] = max(acc_rep);
bacc_rep = i1;
sens_rep = ROCout_rep.yr(i2);
spec_rep = 1-ROCout_rep.xr(i2);


