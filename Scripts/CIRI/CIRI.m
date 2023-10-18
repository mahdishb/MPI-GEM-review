clear
clc

% reading GS
[~,RXNs,~] = xlsread('CompetitiveInhibitors.xlsx','F2:F312');
[~,REGs,~] = xlsread('CompetitiveInhibitors.xlsx','I1:I189');

% reading fingerprints for molecules
[~,REGs_CHEBI,~] = xlsread('CompetitiveInhibitors.xlsx','Regulators','A1:A454');
REGs_CID = xlsread('CompetitiveInhibitors.xlsx','Regulators','B1:B454');

FP_Subs = xlsread('CompetitiveInhibitors.xlsx','FP_Substrates','A1:A134');
FP_Regs = xlsread('CompetitiveInhibitors.xlsx','FP_Regulators','A1:A408');

FP_Subs_Vals = xlsread('CompetitiveInhibitors.xlsx','FP_Substrates','B1:AMK134');
FP_Regs_Vals = xlsread('CompetitiveInhibitors.xlsx','FP_Regulators','B1:AMK408');

emptyCells_rxn = cellfun(@isempty,RXNs);
emptyCells_reg = cellfun(@isempty,REGs);
RXNs(find(emptyCells_reg + emptyCells_rxn == 1)) = [];
REGs(find(emptyCells_reg + emptyCells_rxn == 1)) = [];

% the metabolic model
model = readCbModel;

RXNs_unique = unique(RXNs);
REGs_unique = unique(REGs);

RXNs_index_in_unique = zeros(length(RXNs),1);
for r = 1:length(RXNs)
    for r_u = 1:length(RXNs_unique)
        if strcmp(RXNs_unique{r_u},RXNs{r})
            RXNs_index_in_unique(r) = r_u;
            break
        end
    end
    
end

REGs_index_in_unique = zeros(length(REGs),1);
for r = 1:length(REGs)
    for r_u = 1:length(REGs_unique)
        if strcmp(REGs_unique{r_u},REGs{r})
            REGs_index_in_unique(r) = r_u;
            break
        end
    end
end

% Constructing the gold standard network
GS = zeros(length(RXNs_unique),length(REGs_unique));
for r = 1:length(RXNs_index_in_unique)
    GS(RXNs_index_in_unique(r),REGs_index_in_unique(r)) = 1;
end


% Positive features construction
inx_positive = 0;

[r_pos,c_pos] = find(GS == 1);

for r = 1:length(r_pos)
    rxn = RXNs_unique(r_pos(r));
    SUBs_in_model = [];
    for rr = 1:length(model.rxns)
        if strcmpi(rxn,model.rxns{rr})
            SUBs_in_model = find(model.S(:,rr) < 0);
            if model.lb(rr) < 0
                SUBs_in_model = union(SUBs_in_model,(find(model.S(:,rr) > 0)));
            end
            break
        end
    end
    
    SUBs_in_FP = [];
    ss = 0;
    for s = 1:length(SUBs_in_model)
        if ~isempty(find(FP_Subs == SUBs_in_model(s)))
            ss = ss + 1;
            SUBs_in_FP{ss} = find(FP_Subs == SUBs_in_model(s));
        end
    end
    
    if length(SUBs_in_FP) > 0
        Inx_r = strfind(REGs_CHEBI,REGs_unique{c_pos(r)});
        Index_r = find(not(cellfun('isempty',Inx_r)));
        if ~isempty(Index_r)
            if ~isempty(find(FP_Regs == REGs_CID(Index_r)))
                
                TM_temp = zeros(length(SUBs_in_FP),1);
                for i = 1:length(SUBs_in_FP)
                    FP_s = FP_Subs_Vals(SUBs_in_FP{i},:);
                    FP_r = FP_Regs_Vals(find(FP_Regs == REGs_CID(Index_r)),:);
                    
                    N_s = sum(FP_s);
                    N_r = sum(FP_r);
                    N_sr = sum(FP_s .* FP_r);
                    
                    
                    TM_temp(i) = 1 - (N_sr/(N_s + N_r - N_sr + eps));
                end
                [~,ii] = min(TM_temp);
                FP_s = FP_Subs_Vals(SUBs_in_FP{ii},:);
                FP_r = FP_Regs_Vals(find(FP_Regs == REGs_CID(Index_r)),:);
                
                inx_positive = inx_positive + 1;
                Features_positive(inx_positive,:) = [FP_s,FP_r];
                
                
            end
        end
    end
end

% Negative features construction
inx_negative = 0;

[r_neg,c_neg] = find(GS == 0);

for r = 1:length(r_neg)
    rxn = RXNs_unique(r_neg(r));
    SUBs_in_model = [];
    for rr = 1:length(model.rxns)
        if strcmpi(rxn,model.rxns{rr})
            SUBs_in_model = find(model.S(:,rr) < 0);
            if model.lb(rr) < 0
                SUBs_in_model = union(SUBs_in_model,(find(model.S(:,rr) > 0)));
            end
            break
        end
    end
    
    SUBs_in_FP = [];
    ss = 0;
    for s = 1:length(SUBs_in_model)
        if ~isempty(find(FP_Subs == SUBs_in_model(s)))
            ss = ss + 1;
            SUBs_in_FP{ss} = find(FP_Subs == SUBs_in_model(s));
        end
    end
    
    if length(SUBs_in_FP) > 0
        Inx_r = strfind(REGs_CHEBI,REGs_unique{c_neg(r)});
        Index_r = find(not(cellfun('isempty',Inx_r)));
        if ~isempty(Index_r)
            if ~isempty(find(FP_Regs == REGs_CID(Index_r)))
                for i = 1:length(SUBs_in_FP)
                    FP_s = FP_Subs_Vals(SUBs_in_FP{i},:);
                    FP_r = FP_Regs_Vals(find(FP_Regs == REGs_CID(Index_r)),:);
                    
                    inx_negative = inx_negative + 1;
                    
                    Features_negative(inx_negative,:) = [FP_s,FP_r];
                end
            end
        end
    end
end


%% several SVMs for predicting negative labels
% The labels for negatives are already in the folder labels.mat, but here
% is where the procedure for choice of negatives.

% disp('Prediction of negative samples...     (takes a while)');
% 
% rand_num_Pos = randperm(size(Features_positive,1));
% rand_num_Neg = randperm(size(Features_negative,1));
% 
% STR = [num2str(floor(length(rand_num_Neg)/length(rand_num_Pos))),' SVMs are going to be trained for predicting negative samples.'];
% disp(STR);
% 
% for i = 1: floor(length(rand_num_Neg)/length(rand_num_Pos))
%     Start_train_neg(i) = 1 + length(rand_num_Pos) * (i - 1);
%     end_train_neg(i) = length(rand_num_Pos) * i;
% end

% Labels = zeros(length(rand_num_Neg),1);
%
% for i = 1:length(Start_train_neg)
%
%     STR = ['SVM number ', num2str(i),' for predicting negative samples.'];
%     disp(STR);
%
%     Data_Train = [Features_positive(rand_num_Pos(1:end),:) ; Features_negative(rand_num_Neg(Start_train_neg(i):end_train_neg(i)),:)];
%
%     Group_Train = [ones(length(rand_num_Pos),1) ; zeros(length(rand_num_Pos),1)];
%
%    SVMModel_h = fitcsvm(Data_Train,Group_Train,'Standardize',true,'KernelFunction','rbf',...
%        'KernelScale','auto','OptimizeHyperparameters','all', ...
%        'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName', ...
%        'expected-improvement-plus','Holdout',0.1,'ShowPlots',false));
%
%     range_test_neg = setdiff(1:length(rand_num_Neg),rand_num_Neg(Start_train_neg(i):end_train_neg(i)));
%     Data_Test_Neg = [Features_negative(range_test_neg,:)];
%
%     [label,score] = predict(SVMModel_h,Data_Test_Neg);
%
%     Labels(range_test_neg) = label + Labels(range_test_neg);
% end

load Labels
Neg_labels = find(Labels == 0);

%% Final svm
Train_number = floor(9 * size(Features_positive,1)/10);
Test_number = size(Features_positive,1) - Train_number;
num_of_training = 10;
accuracy_h = zeros(num_of_training,1);
AUCsvm_h = zeros(num_of_training,1);
PRsvm_h = zeros(num_of_training,1);
X_PR = zeros(2*Test_number,num_of_training);
Y_PR = zeros(2*Test_number,num_of_training);
X_AUC = zeros(2*Test_number,num_of_training);
Y_AUC = zeros(2*Test_number,num_of_training);

i = 1;
while i <= num_of_training
    i
    rand_num_Pos = randperm(size(Features_positive,1));
    rand_Neg = randperm(size(Neg_labels,1));
    
    %  Training data is selected randomly:the training set is balanced- having the same number of positives
    %  and negatives
    Data_Train = [Features_positive(rand_num_Pos(1:Train_number),:) ; Features_negative(Neg_labels(rand_Neg(1:Train_number)),:)];
    
    Group_Train = [ones(Train_number,1) ; zeros(Train_number,1)];
    
    SVMModel_h = fitcsvm(Data_Train,Group_Train,'Standardize',true,'KernelFunction','rbf',...
        'KernelScale','auto','OptimizeHyperparameters','all', ...
        'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName', ...
        'expected-improvement-plus','Holdout',0.1,'ShowPlots',false));
    
    % Test data
    Data_Test = [Features_positive(rand_num_Pos(Train_number+1:end),:) ; Features_negative(Neg_labels(rand_Neg(Train_number+1:size(Features_positive,1))),:)];
    Test_number = size(Features_positive,1) - Train_number;
    Group_Test = [ones(Test_number,1) ; zeros(Test_number,1)];
        
    [label,score] = predict(SVMModel_h,Data_Test);
    
    if length(perfcurve(logical(Group_Test),score(:,logical(SVMModel_h.ClassNames)),'true')) == 2*Test_number
        accuracy_h(i,1) = sum(predict(SVMModel_h, Data_Test) == Group_Test)/length(Group_Test)*100;

        [X_AUC(:,i),Y_AUC(:,i),~,AUCsvm_h(i,1)] = perfcurve(logical(Group_Test),score(:,logical(SVMModel_h.ClassNames)),'true');
        [X_PR(:,i),Y_PR(:,i),~,PRsvm_h(i,1)] = perfcurve(logical(Group_Test),score(:,logical(SVMModel_h.ClassNames)),'true','xCrit', 'reca', 'yCrit', 'prec');
        i = i + 1;
    end
    
end
%% Display results for the classifier with best AUC
[~,inx] = max(AUCsvm_h)
STR1 = ['The area under ROC is : ', num2str(AUCsvm_h(inx))];
disp(STR1);

STR2 = ['The area under PR curve is : ', num2str(PRsvm_h(inx))];
disp(STR2);   

figure(1)
hold on
t = linspace(realmin ( 'single' ),1);
plot(t,t,'--','Color',[0,0.45,0.74]);
plot(X_AUC(:,inx),Y_AUC(:,inx),'Color',[0.50,0.50,0.50]);

xlabel('False positive rate')
ylabel('True positive rate')
title('ROC')

figure(2)
hold on
t = linspace(realmin ( 'single' ),1);
plot(t,1-t,'--','Color',[0,0.45,0.74]);
plot(X_PR,Y_PR,'Color',[0.50,0.50,0.50]);

xlabel('Recall')
ylabel('Precision')
title('PR')

% Learning curves
% For checking the learning curves you can uncomment this part

m = 20;
errCV = zeros(m,30);
errTrain = zeros(m,30);
for r = 1:m
    r
    rand_num_Pos = randperm(size(Features_positive,1));
    rand_Neg = randperm(size(Neg_labels,1));
    for i = 1:30
        Data_Train_l = [Features_positive(rand_num_Pos(1:Train_number - (i-1) * 5),:) ; Features_negative(Neg_labels(rand_Neg(1:Train_number - (i-1) * 5)),:)];
        Group_Train_l = [ones(Train_number - (i-1) * 5,1) ; zeros(Train_number - (i-1) * 5,1)];
        
    SVMModel_h = fitcsvm(Data_Train,Group_Train,'Standardize',true,'KernelFunction','rbf',...
        'KernelScale','auto','OptimizeHyperparameters','all', ...
        'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName', ...
        'expected-improvement-plus','Holdout',0.1,'ShowPlots',false));

         errCV(r,i) = (1 - sum(str2num(cell2mat(predict(SVMModel_h, Data_Test))) == Group_Test)/length(Group_Test))*100;
        
        errTrain(r,i) = (1 - sum(str2num(cell2mat(predict(SVMModel_h, Data_Train_l))) == Group_Train_l)/length(Group_Train_l))*100;
        
    end
end

errCV_plot = zeros(m,30);
errTrain_plot = zeros(m,30);

for j = 1:size(errCV,1)
    for i = 1:size(errCV,2)
        errCV_plot(j,i) = errCV(j,length(errCV)-i+1);
        errTrain_plot(j,i) = errTrain(j,length(errTrain)-i+1);
    end
end
%
figure(1)
hold on

plot(20:10:310,mean(errCV_plot),'b-','linewidth',2);
plot(20:10:310,mean(errTrain_plot),'r-','linewidth',2);

y_train_1 = mean(errTrain_plot) - std(errTrain_plot)/2;
y_train_2 = mean(errTrain_plot) + std(errTrain_plot)/2;
y_CV_1 = mean(errCV_plot) - std(errCV_plot)/2;
y_CV_2 = mean(errCV_plot) + std(errCV_plot)/2;

x = 20:10:310;
patch ([x fliplr(x)], [y_train_1 fliplr(y_train_2)], [1 0 0]*0.5, 'EdgeColor','none', 'FaceAlpha',0.2 )
patch ([x fliplr(x)], [y_CV_1 fliplr(y_CV_2)],  [0 0 1]*0.5, 'EdgeColor','none', 'FaceAlpha',0.2)

plot(20:10:310,y_train_1,'r--','linewidth',0.5);
plot(20:10:310,y_train_2,'r--','linewidth',0.5);
plot(20:10:310,y_CV_1,'b--','linewidth',0.5);
plot(20:10:310,y_CV_2,'b--','linewidth',0.5);
hold off

