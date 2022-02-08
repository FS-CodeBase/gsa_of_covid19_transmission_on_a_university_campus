% fun_compute_model_solutions
% SCRIPT: compute_model_solutions                                    
% AUTHOR: Fabian Santiago                                             
% EMAIL: fabiansantiago707@gmail.com                                  
% DATE: 11/27/2021                                                    
%     Computes needed model solutions for computation of the FIRST    
%     ORDER and TOTAL ORDER sensitivity indices of parameters COV-19  
%     Model indices following (Saltelli et al. 2008 notation)         

% Clear workspace and close all figures
% clear
% close all

% Start parpool for parfor
% parpool('local',20);

% Load parameters of interest and their ranges for Latin Hypercube Sampling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameters_info = fun_model_parameter_ranges;

% Load contact matrices and number of students by classification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ConMat,Con_livingMat,...
    nUVec,nDVec,nGVec,nFVec] = fun_initialize_contact_matrices;

% Model solution length in days
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DaysRefined = 30; % Number of days model is solved refined
MaxDays = 7*15; % End of semester (15 week semester)

% Refine solutions during the first 30 days, then every 3 days
tspan = 24*[0:DaysRefined (DaysRefined+3):3:MaxDays];
converttohrs = 24; % 24 hrs in a day

% Particular cases
%       case  Vu   Vg
cases = [24   0    0;...
         24   40   50;...
         24   80   100;...
         23   0    0;...
         23   40   50;...
         23   80   100;...
         22   0    0;...
         22   40   50;...
         22   80   100];

% Number of Monte Carlo (MC) base parameter samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_mc_samples = 1200*(sum(parameters_info(:,1))+2);

% Generate MC Parameter Value Matrices: A, B
A = fun_lhs_sampling_of_parameters(parameters_info,n_mc_samples);
B = fun_lhs_sampling_of_parameters(parameters_info,n_mc_samples);

% Pre-allocate space for matrices B_Aj and A_Bj. Notations indicates 
% matrix A or B from above but the jth column is replaced with that of the 
% other matrix.
B_Aj_cell = cell(size(parameters_info,1),1); 
A_Bj_cell = cell(size(parameters_info,1),1);

% Define a matrix Ci formed by all columns of B except the ith column,
% which is taken from A
for param = 1:size(parameters_info,1)
    % Only parameters whose sensitivity is being estimated
    if logical(parameters_info(param,1))
        ABjtmp = A;
        BAjtmp = B;
       
        % Matrix, where column J comes from matrix B and 
        % all other k − 1 columns come from matrix A
        ABjtmp(:,param) = B(:,param); 
        
        % Matrix, where column J comes from matrix A and 
        % all other k − 1 columns come from matrix B
        BAjtmp(:,param) = A(:,param); 
        
        % Save A_Bj and B_Aj
        A_Bj_cell{param} = ABjtmp;
        B_Aj_cell{param} = BAjtmp;
    end
end

% Compute model solutions
%%%%%%%%%%%%%%%%%%%%%%%%%
for faculty_vax = [0 50 100] % Faculty Vaccination
for contact_scenario = 22:24 % Class cap: 50,100,200
    disp(['CONTACT SCENARIO: ',num2str(contact_scenario)])
    % Percent Faculty and Graduate Students Vaccinated 
    Vf = faculty_vax; Vg = faculty_vax; 
    for undergraduate_vax = [0 40 80] %(2*vax_case:(2*vax_case+1))*10 
        if any(ismember(cases,[contact_scenario undergraduate_vax faculty_vax],'rows'))

        % Percent of undergraduates vaccinated
        Vu = undergraduate_vax; Vd = undergraduate_vax; 
        
        % Load a particular contact scenario
        [Con,Con_living,nu,nd,ng,nf] = deal(...
                            ConMat(:,:,contact_scenario),...
                                Con_livingMat(:,:,contact_scenario),...
                                    nUVec(contact_scenario),...
                                        nDVec(contact_scenario),...
                                            nGVec(contact_scenario),...
                                                nFVec(contact_scenario)...
                                                    );                                              
        % Compute YA
        YA_cell = solve_model(A,tspan,converttohrs,Con,Con_living,...
                                nu,nd,ng,nf,Vu,Vd,Vg,Vf);

        % Compute YB
        YB_cell = solve_model(B,tspan,converttohrs,Con,Con_living,...
                                nu,nd,ng,nf,Vu,Vd,Vg,Vf);

        % Compute YA_Bj
        YA_Bj_cell = cell(1,numel(parameters_info(:,1)));
        for param = find(parameters_info(:,1))'
            A_Bj_tmp = A_Bj_cell{param};
            YA_Bj_tmp = solve_model(...
                                    A_Bj_tmp,tspan,converttohrs,...
                                        Con,Con_living,nu,nd,ng,nf,...
                                            Vu,Vd,Vg,Vf...
                                                );
            YA_Bj_cell{param} = YA_Bj_tmp;
        end

        % Save the solutions using case specific information
        save(... % File name:
                ['./model_sols/model_sols_dailyMCn',...
                    num2str(n_mc_samples),...
                        'Vu',num2str(Vu),'Vd',num2str(Vd),...
                            'Vg',num2str(Vg),'Vf',num2str(Vf),...
                                'CS_',num2str(contact_scenario),'.mat'],...
                                    'YA_cell','YB_cell','YA_Bj_cell',...
            ... % Variables:
                'parameters_info',...
                    'n_mc_samples',...
                        'Vu','Vd','Vg','Vf',...
                            'contact_scenario',...
                                '-v7.3');
        end
    end % END FOR-Loop: CONTACT SCENARIO
end
end