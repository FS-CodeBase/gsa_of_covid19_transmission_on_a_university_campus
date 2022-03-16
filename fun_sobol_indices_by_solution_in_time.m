function [FirstOrderIdx,TotalOrderIdx,Stats] = ...
                            fun_sobol_indices_by_solution_in_time(...
                                y_sol_idx,...
                                    sub_sample_size,...
                                        n_resamples,...
                                            YA_cell,...
                                                YB_cell,...
                                                    YA_Bj_cell,...
                                                        prms_bool_info)
% FUNCTION: fun_sobol_indices_by_solution_in_time
% AUTHOR: Fabian Santiago
% EMAIL: fsantiago3@ucmerced.edu
% DATE: 8/27/2021
%     Computes the FIRST ORDER and TOTAL ORDER sensitivity indices of 
%     parameters COV-19 Model indices following Saltelli et al. 2008
% clear
% close all
% 
% % Load Y_A_cell, YB_cell, YC_cell, prms_info, and prms_srt
% load('model_sols/model_solsN4096matYAYB.mat');
% load('model_sols/model_solsN4096matYA_Bj_cell.mat');
% load('model_sols/model_solsN4096matYB_Aj_cell.mat');

% Compute first order indicies in time
% YA_cell = {yA1;yA2;...;yAN};
% yAi = [Su:1,   Sd:2,   Sg:3,   Sf:4
%        Eu:5,   Ed:6,   Eg:7,   Ef:8,
%        Iau:9,  Isu:10, Iad:11, Isd:12,
%        Iag:13, Isg:14 ,Iaf:14, Isf:16,
% ] 
% F =@(x) mean(x,1);
% MEAS =@(x) min(x,[],1);
% y_sol_idx = 3;
% SubSampN  = N; 
% ReplicasN = 100;
% prm_idices = find(prms_info(:,1));
N = numel(YA_cell); % Number of parameters sampled uniformly
Ntpts = size(YA_cell{1},1);
numPrms = sum(prms_bool_info);
si_prms = find(prms_bool_info)';
Scell   = cell(1,n_resamples);
STicell = cell(1,n_resamples);
Varcell = cell(1,n_resamples);
Mucell  = cell(1,n_resamples);
% tic
parfor R = 1:n_resamples

    % RESAMPLE WITH REPLACEMENT
    resmpl_idx = randsample(y_sol_idx:1:(1*N),sub_sample_size,true);

    YA = cell2mat(YA_cell);
    YA = YA(:,resmpl_idx);
    YB = cell2mat(YB_cell);
    YB = YB(:,resmpl_idx);
    
    % COMPUTE MU AND VAR    
    Mu = mean(YA,2);
    F0sq = Mu.^2;
    yAyA = mean(YA.^2,2); 
    VarY = yAyA-F0sq;
    
    % Save Mu and Var
    Varcell{R} = VarY;
    Mucell{R}  = Mu;
    S   = zeros(Ntpts,numPrms);
    STi = zeros(Ntpts,numPrms);
    si_idx = 1;
    % YA_Bj_cell is a broadcast variable 
    for prm_idx = si_prms
        YA_B = cell2mat(YA_Bj_cell{prm_idx}); % Agree in Jth column
%         YB_A = cell2mat(YB_Aj_cell{prm_idx}); % Disagree in Jth column
        YA_Bj = YA_B(:,resmpl_idx);
%         YB_Aj = MEAS(YB_A(:,resmpl_idx));
        
        %%% FIRST ORDER INDICES
%         S(si_idx) = (mean(YA1.*YB_Aj)-F0sq)/VarY; % Sobol 1993
%         S(si_idx) = mean(YA1.*(YB_Aj-YB1))/VarY; % Saltelli 2002
         S(:,si_idx) = mean(YB.*(YA_Bj-YA),2)./VarY; % Saltelli 2010
%         S(si_idx) = (VarY-1/2*mean((YB1-YA_Bj).^2))/VarY;

        %%% TOTAL ORDER INDICES 
%        STi(si_idx) = (VarY-mean(YA1.*YA_Bj)+F0sq)/VarY; % Homma 1996
%         STi(si_idx) = mean(YA1.*(YA1-YA_Bj))/VarY; % Sobol 2007
        STi(:,si_idx) = 1/2*mean((YA-YA_Bj).^2,2)./VarY; % Jansen 1999
%         STi(si_idx) = 1-mean(YB1.*YB_Aj-F0sq)/VarY; % Saltelli 2008 

        % INCREASE Si INDEX
        si_idx = si_idx + 1;
    end
    Scell{R}   = S;
    STicell{R} = STi;
end

% SENSITIVITIES TO MATRICES
SIMAT = cell2mat(Scell);
STMAT = cell2mat(STicell);

% MEAN AND VARIANCE MATS
MU = mean(cell2mat(Mucell),2)';
VAR = mean(cell2mat(Varcell),2)';

% MEAN and ERROR of FIRST ORDER AND TOTAL ORDER SENSITIVITIES
FO_idx = zeros(Ntpts,numPrms);
FO_err = zeros(Ntpts,numPrms);

TO_idx = zeros(Ntpts,numPrms);
TO_err = zeros(Ntpts,numPrms);

Pvals = zeros(Ntpts,numPrms);
CIs = zeros(Ntpts,numPrms,2);
for P = 1:numPrms
    for T = 1:Ntpts
        % FIRST ORDER SENSITIVITY
        S = sort(SIMAT(T,P:numPrms:end),'ascend');
        ES = S(ceil(n_resamples*0.95))-mean(S);
        FO_idx(T,P) = mean(S);
        FO_err(T,P) = ES;
        
        % TOTAL ORDER SENSITIVITY
        Ti = sort(STMAT(T,P:numPrms:end),'ascend');
        ET = Ti(ceil(n_resamples*0.95))-mean(Ti);
        TO_idx(T,P) = mean(Ti);
        TO_err(T,P) = ET;
        
        % TWO SAMPLE t-test for difference in mean between FO and TO
        [~,pval,CI] = ttest2(SIMAT(T,P:numPrms:end),STMAT(T,P:numPrms:end),'Vartype','unequal');
        Pvals(T,P) = pval; 
        CIs(T,P,:) = CI;
    end
end
FirstOrderIdx.indices = FO_idx; FirstOrderIdx.error = FO_err;
TotalOrderIdx.indices = TO_idx; TotalOrderIdx.error = TO_err;
Stats.mu = MU; Stats.var = VAR;
Stats.pvals = Pvals; Stats.CI = CIs;