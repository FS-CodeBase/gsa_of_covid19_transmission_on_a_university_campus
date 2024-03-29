% SCRIPT: RUN_ME                                    
% AUTHOR: Fabian Santiago                                             
% EMAIL: fsantiago3@ucmerced.edu                                  
% DATE: 11/27/2022  
%
% DESCRIPTION:
%     Script shows how to use the code to run simulations and perform a
%     global sensitivity analysis.

% Define cases to be numerically computed. These are the cases considered in the
% main text.
% CASES = [50000,  0,  0;... % [50000 = No class cap  | No vaccination] (Baseline)
%            100,  0,  0;... % [100 student class cap | No vaccination]
%             50,  0,  0;... % [50 student class cap  | No vaccination]
%          50000, 40, 50;... % [50000 = No class cap  | vU=vD=40 | vG=vF=50 ]
%            100, 40, 50;... % [100 student class cap | vU=vD=40 | vG=vF=50 ]
%             50, 40, 50;... % [50 student class cap  | vU=vD=40 | vG=vF=50 ]
%          50000, 80,100;... % [50000 = No class cap  | vU=vD=80 | vG=vF=100]
%            100, 80,100;... % [100 student class cap | vU=vD=80 | vG=vF=100]
%             50, 80,100];   % [50 student class cap  | vU=vD=80 | vG=vF=100]

% Used to test code. First row is the baseline case and is necessary for adding
% the reduction in cumulative infections (RECI) to plots. 
CASES = [50000,  0, 0;... % [50000 = No class cap  | No vaccination] (Baseline)
            50, 40, 50];  % [50 student class cap  | vU=vD=40 | vG=vF=50 ]

% Get parameter information
prms_info = fun_model_parameter_ranges;

% Number of samples from parameter space using Latin hypercube sampling
SubSampN  = 1500*(sum(prms_info(:,1))+2);

% Reduced infection for 'far' contacts. 
reduced_infect = 25; % Options: 25 (for 25%) or 10 (for 10%)

% Compute model solutions
fun_compute_model_solutions(CASES,reduced_infect,SubSampN)

% Compute Sobol indices for 
% doubling time using all solutions in model_sols folder
compute_sobol_indices_doubling_time

% Compute Sobol indices 
% over time using all solutions in model_sols folder
compute_sobol_indices_over_time

% Plot sensitivity of doubling time
class_scenario = 50; % Class cap scenario
undg_vax = 40; % Undergraduate vaccination
fac_grad_vax = 50; % Faculty vaccination
param_set = 'inf and con'; % Options: 'inf and con' OR 'init'
print_contact_info = true;
print_ylabel = true;
font_size = 16;
plot_sobol_indices_dt(class_scenario,...
                            undg_vax,fac_grad_vax,...
                                param_set,...
                                    print_contact_info,...
                                        print_ylabel,...
                                            SubSampN,...
                                                font_size)

% Plot sensitivity of cumulative infections over time
class_scenario = 50; % Class cap scenario
undg_vax = 40; % Undergraduate vaccination
fac_grad_vax = 50; % Faculty vaccination
param_set = 'init'; % Options: 'inf', 'con', OR 'init'
print_contact_info = true;
print_ylabel = true;
font_size = 16;
Sobol_order = 'FO'; % OPTIONS: 'TE' OR 'FO'
plot_sobol_indices_time(class_scenario,...
                            undg_vax,fac_grad_vax,...
                                param_set,...
                                    print_contact_info,...
                                        print_ylabel,...
                                            Sobol_order,...
                                                SubSampN,...
                                                    font_size)
