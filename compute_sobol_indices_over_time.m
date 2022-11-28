% SCRIPT: compute_sobol_indices_doubling_time
% AUTHOR: Fabian Santiago
% EMAIL: fsantiago3@ucmerced.edu
% DATE: 8/27/2021

% Make a folder to store Sobol indices if one does not exist
if ~exist('sobol_indices_in_time', 'dir')
       mkdir('sobol_indices_in_time')
end

% Load parameter information for global sensitivity analysis
parameters_info = fun_model_parameter_ranges;

% Load model solutions for each contact scenario
model_sols_file_id = dir('model_sols/');
model_sols_file_id(1:2) = [];

% Set sub-sampling size for bootstrapping and number of resamples
SubSampN  = 1500*(sum(parameters_info(:,1))+2); % Sub-samples to compute Si from resampled solutions
n_resamples = 2000;                                    % Number of times to perform re-sampling
y_sol_idx = 1; 

% Compute Sensitivity Indices
for case_solution_idx = 1:numel(model_sols_file_id)
tic
% Load model solutions YA, YB, YA_Bj, and
% parameter information prms_info and prms_str
load(['model_sols/',model_sols_file_id(case_solution_idx).name])

% Display contact scenario being considered
disp(['contact scenario: ',num2str(class_cap)]);

% Compute first and total-order Sobol indices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
[FirstOrderIdx,TotalOrderIdx,Stats] = ...
    fun_sobol_indices_by_solution_over_time(...
        y_sol_idx,...
            SubSampN,...
                n_resamples,...
                    YA_cell,...
                        YB_cell,...
                            YA_Bj_cell,...
                                parameters_info(:,1));
                                                    
save(['sobol_indices_in_time/SobolT_',...
        'sol',num2str(y_sol_idx),...
            'sub',num2str(SubSampN),...
                'rep',num2str(n_resamples),...
                    'CC_',num2str(class_cap),...
                        'Vu',num2str(Vu),'Vd',num2str(Vd),...
                            'Vg',num2str(Vg),'Vf',num2str(Vf),'.mat'],...
                    'FirstOrderIdx','TotalOrderIdx','Stats','class_cap');
disp(['cs: ',num2str(class_cap),', time: ',num2str(toc/60),'min']);
end
