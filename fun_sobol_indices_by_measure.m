function [FirstOrderIdx,TotalOrderIdx,Stats] = ...
            fun_sobol_indices_by_measure(...
                y_sol_idx,...
                    sub_sample_size,...
                        n_resamples,...
                            YA_cell,...
                                YB_cell,...
                                    YA_Bj_cell,...
                                        prms_bool_info)
% % FUNCTION: FUN_SOBOL_INDICES_BY_MEASURE
% % AUTHOR: Fabian Santiago
% % EMAIL: fsantiago3@ucmerced.edu
% % DATE: 11/16/2020
% %     Computes the FIRST ORDER and TOTAL ORDER sensitivity indices of 
% %     parameters COV-19 Model indices following Saltelli et al. 2008
% %     notation, using the fun_measure 

n_solutions = numel(YA_cell);
sol_dim = size(YA_cell{1},2);
n_paramaters = sum(prms_bool_info);
Si_prm_idx = find(prms_bool_info)';
Scell   = cell(n_resamples,1);
STicell = cell(n_resamples,1);

DT_Var = ones(1,n_resamples);
DT_Mu  = ones(1,n_resamples);

CI_Var = ones(1,n_resamples);
CI_Mu  = ones(1,n_resamples);

% tic
parfor resamp_idx = 1:n_resamples
    % Resample with replacement (true flag)
    resmpl_idx = randsample(y_sol_idx:sol_dim:(sol_dim*n_solutions),...
                                                    sub_sample_size,true);
    YA = cell2mat(YA_cell);  YA = (YA(:,resmpl_idx));
    YB = cell2mat(YB_cell);  YB = (YB(:,resmpl_idx));

    % Compute mean and std for Cumulative Infections (CIs)
    mu_ci = mean(YA(end,:),2);
    F0sq = mu_ci.^2;
    yCIyCI = mean(YA(end,:).^2,2); 
    var_ci = yCIyCI-F0sq;
    
    % Save mean and std for CIs
    CI_Mu(resamp_idx) = mu_ci; CI_Var(resamp_idx) = var_ci;
    
    % Compute Doubling Times (DT)
    YB = fun_doubling_time(YB);
    YA = fun_doubling_time(YA);    

    % Compute mu and std for DT
    mu_dt  = mean(YA,2);
    F0sq   = mu_dt.^2;
    yAyA   = mean(YA.^2,2); 
    var_dt = yAyA-F0sq;
    
    % Save mu and std for DT
    DT_Mu(resamp_idx) = mu_dt; DT_Var(resamp_idx) = var_dt;
        
    % SENSITIVITY MATRIX
    S   = zeros(1,n_paramaters);
    STi = zeros(1,n_paramaters);
    
    si_idx = 1;
    for param_idx = Si_prm_idx
        YA_B = cell2mat(YA_Bj_cell{param_idx}); % Agree in Jth column
%         YB_A = cell2mat(YB_Aj_cell{prm_idx}); % Disagree in Jth column
        YA_Bj = fun_doubling_time(YA_B(:,resmpl_idx));
%         YB_Aj = MEAS(YB_A(:,resmpl_idx));
        
        %%% FIRST ORDER INDICES
%         S(si_idx) = (mean(YA1.*YB_Aj)-F0sq)/VarY; % Sobol 1993
%         S(si_idx) = mean(YA1.*(YB_Aj-YB1))/VarY; % Saltelli 2002
         S(si_idx) = mean(YB.*(YA_Bj-YA))./var_dt; % Saltelli 2010
%         S(si_idx) = (VarY-1/2*mean((YB1-YA_Bj).^2))/VarY;

        %%% TOTAL ORDER INDICES 
%        STi(si_idx) = (VarY-mean(YA1.*YA_Bj)+F0sq)/VarY; % Homma 1996
%         STi(si_idx) = mean(YA1.*(YA1-YA_Bj))/VarY; % Sobol 2007
        STi(si_idx) = 1/2*mean((YA-YA_Bj).^2)./var_dt; % Jansen 1999
%         STi(si_idx) = 1-mean(YB1.*YB_Aj-F0sq)/VarY; % Saltelli 2008 

        % Increase i index in Si
        si_idx = si_idx + 1;
    end
    Scell{resamp_idx}   = S;
    STicell{resamp_idx} = STi;
end

% Allocate space for sensitivity matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sensitivity matrices
Si_mat  = cell2mat(Scell);
STi_mat = cell2mat(STicell);

% First order Sobol Indices
FO_idx = zeros(1,n_paramaters);
FO_err = zeros(1,n_paramaters);

% Total order Sobol Indices
TO_idx = zeros(1,n_paramaters);
TO_err = zeros(1,n_paramaters);

Pvals = zeros(1,n_paramaters);
CIs = zeros(n_paramaters,2);
for P = 1:n_paramaters
    % FIRST ORDER SENSITIVITY
    S = sort(Si_mat(:,P),'ascend');
    MuS = mean(S);
    ES = S(ceil(n_resamples*0.95))-MuS;
    FO_idx(P) = MuS;
    FO_err(P) = ES;

    % TOTAL ORDER SENSITIVITY
    Ti = sort(STi_mat(:,P),'ascend');
    MuTi = mean(Ti);
    ET = Ti(ceil(n_resamples*0.95))-MuTi;
    TO_idx(P) = MuTi;
    TO_err(P) = ET;
    
    % TWO SAMPLE Whelch's t-test for difference in mean between FOI and TOI
    [~,pval,CI] = ttest2(STi_mat(:,P),Si_mat(:,P),'Vartype','unequal');
    Pvals(P) = pval; 
    CIs(P,:) = CI;
end

% Store first order indices, total order indices, and stats
FirstOrderIdx.indices = FO_idx; FirstOrderIdx.error = FO_err;
TotalOrderIdx.indices = TO_idx; TotalOrderIdx.error = TO_err;
Stats.mu.CI = CI_Mu; Stats.var.CI = CI_Var;
Stats.mu.DT = DT_Mu; Stats.var.DT = DT_Var;
Stats.pvals = Pvals; Stats.CI = CIs;
