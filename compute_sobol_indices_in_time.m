% SCRIPT: compute_sobol_indices_doubling_time
% AUTHOR: Fabian Santiago
% EMAIL: fabiansantiago707@gmail.com
% DATE: 8/27/2021
% c = parpool('local',16); 

% Load parameter information for global sensitivity analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameters_info = fun_model_parameter_ranges;

% % Load model solutions for each contact scenario
model_sols_file_id = dir('model_solutions/');
model_sols_file_id(1:2) = [];

% Set sub-sampling size for bootstrapping and number of resamples
sub_sample_size  = 1200*(sum(parameters_info(:,1))+2); % Sub-samples to compute Si from resampled solutions
n_resamples = 2000; % Number of times to perform re-sampling
y_sol_idx = 1; 

% Compute Sensitivity Indices
for case_solution_idx = 1:numel(model_sols_file_id)
tic
% Load model solutions YA, YB, YA_Bj, and
% parameter information prms_info and prms_str
load(['model_solutions/',model_sols_file_id(case_solution_idx).name])

% Display contact scenario being considered
disp(['contact scenario: ',num2str(contact_scenario)]);

% Compute first and total-order Sobol indices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
[FirstOrderIdx,TotalOrderIdx,Stats] = ...
    fun_sobol_indices_by_solution_in_time(...
        y_sol_idx,...
            sub_sample_size,...
                n_resamples,...
                    YA_cell,...
                        YB_cell,...
                            YA_Bj_cell,...
                                parameters_info(:,1));
                                                    
save(['sobol_indices_in_time/SobolT_',...
        'sol',num2str(y_sol_idx),...
            'sub',num2str(sub_sample_size),...
                'rep',num2str(n_resamples),...
                    'CS_',num2str(contact_scenario),...
                        'Vu',num2str(Vu),'Vd',num2str(Vd),...
                            'Vg',num2str(Vg),'Vf',num2str(Vf),'.mat'],...
                    'FirstOrderIdx','TotalOrderIdx','Stats','contact_scenario');
disp(['cs: ',num2str(contact_scenario),', time: ',num2str(toc/60),'min']);
end