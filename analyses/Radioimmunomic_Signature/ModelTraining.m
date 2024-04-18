%%%%% Author: Anahita Fathi Kazerooni %%%%%%
%%%%% Email: anahita.fathi@gmail.com
%%%%% GitHub handle: anahitafk
%%%%% Citation: ""


%% add paths to the functions
clear; clc; close all;
addpath('dnafinder-roc-b84f535');


%% Load Data
load('Features_ModelTraining.mat', 'All_Features_discovery', 'class_discovery');


%% Settings
rng(1)
cv = 10;
kernel_type = 'linear';

X = All_Features_discovery;
y = class_discovery;


CV1 = cvpartition(y,'KFold',cv,'Stratify',true); % outer CV

selected_ranking_outer = struct;
for idx1 = 1:CV1.NumTestSets
    disp(num2str(idx1))
    clearvars -except idx1 CV1  y X selected_ranking_outer ...
         All_Features_discovery selection_method kernel_type ...
         All_Features_replicate class_discovery class_replicate ...
        Updated_FeaturesList Gene_Name cv
    X_train = double( X(CV1.training(idx1),:) );
    Y_train = double(logical( y(CV1.training(idx1)) ));
    X_test = double( X(CV1.test(idx1),:) );
    Y_test = double(logical( y(CV1.test(idx1)) ));
    
    CV2 = cvpartition(Y_train, 'KFold', 5, 'Stratify',true); %inner CV
    ranking_mat = [];

    [ranking,scores] = fscmrmr(X_train, Y_train);
    
    selected_ranking_inner = struct;
    auc_inner_cv_save = [];
    auc_inner_train_save = [];
    for idx2 = 1:CV2.NumTestSets
        clearvars -except idx1 CV1 y X selected_ranking_outer All_Features_discovery ...
            cv selection_method X_train Y_train X_test Y_test  CV2 numF ranking ...
            All_Features_replicate class_discovery class_replicate ...
            selected_ranking_inner auc_inner_cv_save  auc_inner_train_save ...
            idx2 kernel_type Updated_FeaturesList Gene_Name
        
        X_CV_training = X_train(CV2.training(idx2),:);
        Y_CV_training = Y_train(CV2.training(idx2),:);
        X_CV_test = X_train(CV2.test(idx2),:);
        Y_CV_test = Y_train(CV2.test(idx2),:);
        
        ranking_fold = struct;
        auc_cv_save = [];
        auc_train_save = [];
        k_save = [];
        c = cvpartition(Y_CV_training,'KFold',3, 'Stratify', true);
        jj=5; 
        for  k = jj:1:50 % select the first n features
            
            cdata = X_CV_training(:,ranking(1:k));
            grp = Y_CV_training;
            % Use a linear or kernel support vector machine classifier
            sigma = optimizableVariable('sigma',[1,1e3],'Transform','log');
            box = optimizableVariable('box',[1,1e3],'Transform','log');
            
            minfn = @(z)kfoldLoss(fitcsvm(cdata,grp,'CVPartition',c,...
                'KernelFunction',kernel_type,'KernelScale',z.sigma));
            results = bayesopt(minfn,sigma,'IsObjectiveDeterministic',true,...
                'AcquisitionFunctionName','expected-improvement-plus',...
                'Verbose',0, 'PlotFcn', [], 'UseParallel',true);
            z(1) = results.XAtMinObjective.sigma;
            SVMModel = fitcsvm(cdata,grp,'KernelFunction',kernel_type,...
                'KernelScale',z(1));
                    
                       
            
           [label, score] = predict(SVMModel,X_CV_test(:,ranking(1:k)));
            if score(1,1)<0
                sscore_cv = score(:,2);
            else
                sscore_cv = score(:,1);
            end
            
            y_cv = zeros(size(Y_CV_test));
            y_cv(find(Y_CV_test==1)) = 1;
            ROCout_CV = roc_dnafinder([sscore_cv,y_cv],0,0.05,0);
            auc_cv = ROCout_CV.AUC;
            
            
            [label, score] = predict(SVMModel,X_CV_training(:,ranking(1:k)));
            if score(1,1)<0
                sscore_train = score(:,2);
            else
                sscore_train = score(:,1);
            end
            y_train = zeros(size(Y_CV_training));
            y_train(find(Y_CV_training==1)) = 1;
            ROCout_train = roc_dnafinder([sscore_train,y_train],0,0.05,0);
            auc_train = ROCout_train.AUC;
            
            ranking_fold(k-jj+1).ranking = ranking(1:k);
            ranking_fold(k-jj+1).SVM_HypParam = results;
            ranking_fold(k-jj+1).SVMModel = SVMModel;
            ranking_fold(k-jj+1).auc_cv = auc_cv;
            ranking_fold(k-jj+1).auc_train = auc_train;
            ranking_fold(k-jj+1).FNames = Updated_FeaturesList(ranking(1:k));
            
            auc_cv_save = [auc_cv_save; auc_cv];
            auc_train_save = [auc_train_save; auc_train];
            k_save = [k_save; k];
            
            
        end
        [~,ind] = max(auc_cv_save);
        if length(ind) ==1
            k_max = k_save(ind); 
        else
            k_max = k_save(ind(1));
        end
        
        selected_ranking = ranking_fold(k_max-jj+1).ranking;
        
       
        
        cdata = X_CV_training(:, selected_ranking);
        grp = Y_CV_training;
       
        SVMModel = ranking_fold(k_max-jj+1).SVMModel;
           
        [label, score] = predict(SVMModel,X_CV_training(:,selected_ranking));
        if score(1,1)<0
            sscore_tr = score(:,2);
        else
            sscore_tr = score(:,1);
        end
        
        y_tr = zeros(size(Y_CV_training));
        y_tr(find(Y_CV_training==1)) = 1;
        ROCout_tr = roc_dnafinder([sscore_tr,y_tr],0,0.05,0);
        auc_tr = ROCout_tr.AUC;
        auc_inner_train_save = [auc_inner_train_save; auc_tr];
        
        
        
        [label, score] = predict(SVMModel,X_CV_test(:,selected_ranking));
        if score(1,1)<0
            sscore_cv = score(:,2);
        else
            sscore_cv = score(:,1);
        end
        
        y_cv = zeros(size(Y_CV_test));
        y_cv(find(Y_CV_test==1)) = 1;
        ROCout_CV = roc_dnafinder([sscore_cv,y_cv],0,0.05,0);
        auc_cv = ROCout_CV.AUC;
        auc_inner_cv_save = [auc_inner_cv_save; auc_cv];
        %
        %
        selected_ranking_inner(idx2).ranking = selected_ranking;
        selected_ranking_inner(idx2).auc_cv = auc_cv;
        selected_ranking_inner(idx2).auc_tr = auc_tr;
        selected_ranking_inner(idx2).FNames = Updated_FeaturesList(selected_ranking);
        %
    end
    %
    % find optimum hyperparameter options for the selected features

    [~, ind] = max(auc_inner_cv_save);
    selected_ranking2 = selected_ranking_inner(ind).ranking;

    
    
    % we have K(number of folds) models
    c2 = cvpartition(Y_train,'KFold',3, 'Stratify', true);
    cdata = X_train(:, selected_ranking2);
    grp = Y_train;
    sigma = optimizableVariable('sigma',[1,1e3],'Transform','log');
    box = optimizableVariable('box',[1,1e3],'Transform','log');
    
       
    minfn = @(z)kfoldLoss(fitcsvm(cdata,grp,'CVPartition',c2,...
        'KernelFunction',kernel_type,'KernelScale',z.sigma));
    results = bayesopt(minfn,sigma,'IsObjectiveDeterministic',true,...
        'AcquisitionFunctionName','expected-improvement-plus',...
        'Verbose',0, 'PlotFcn', [], 'UseParallel',true);
    z(1) = results.XAtMinObjective.sigma;
    SVMModel = fitcsvm(cdata,grp,'KernelFunction',kernel_type,...
        'KernelScale',z(1));


    
    [label, score] = predict(SVMModel,X_train(:,selected_ranking2));
    if score(1,1)<0
        sscore_tr = score(:,2);
    else
        sscore_tr = score(:,1);
    end
    
    y_tr = zeros(size(Y_train));
    y_tr(find(Y_train==1)) = 1;
    ROCout_tr = roc_dnafinder([sscore_tr,y_tr],0,0.05,0);
    auc_train = ROCout_tr.AUC;
    label_train = label;
    
    [label_test, score] = predict(SVMModel,X_test(:,selected_ranking2));
    if score(1,1)<0
        sscore_test = score(:,2);
    else
        sscore_test = score(:,1);
    end
    
    y_test = zeros(size(Y_test));
    y_test(find(Y_test==1)) = 1;
    ROCout_test = roc_dnafinder([sscore_test,y_test],0,0.05,0);
    auc_test = ROCout_test.AUC;
    
    selected_ranking_outer(idx1).selected_ranking = selected_ranking2;
    selected_ranking_outer(idx1).selected_Hyp = results;
    selected_ranking_outer(idx1).SVMModel = SVMModel;
    selected_ranking_outer(idx1).TestPredictions = label_test;
    selected_ranking_outer(idx1).TrainPredictions = label_train;
    selected_ranking_outer(idx1).X_train = X_train(:,selected_ranking2);
    selected_ranking_outer(idx1).Y_train = Y_train;
    selected_ranking_outer(idx1).ROCout_train = ROCout_tr;
    selected_ranking_outer(idx1).AUC_train = auc_train;
    selected_ranking_outer(idx1).X_test = X_test(:,selected_ranking2);
    selected_ranking_outer(idx1).Y_test = Y_test;
    selected_ranking_outer(idx1).ROCout_test = ROCout_test;
    selected_ranking_outer(idx1).AUC_test = auc_test;
    selected_ranking_outer(idx1).FList = Updated_FeaturesList(selected_ranking2);
    selected_ranking_outer(idx1).inner = selected_ranking_inner;



    save ('Radioimmunomic_adc_mrmr_discovery.mat')
end




