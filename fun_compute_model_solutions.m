function fun_compute_model_solutions(CASES,reduced_infect,SubSampN)
% FUNCTION: FUN_COMPUTE_MODEL_SOLUTIONS                                    
% AUTHOR: Fabian Santiago                                             
% EMAIL: fsantiago3@ucmerced.edu                                  
% DATE: 11/27/2021   
%
% DESCRIPTION:
%     Computes model solutions for computation of the FIRST ORDER 
%     and TOTAL ORDER sensitivity indices of parameters COV-19  
%     Model indices following (Saltelli et al. 2008 notation). 
% 
% INPUTS: 
%     CASES: A vector of class caps: 25, 50, 100 or 50000 (no cap)
%     reduced_infect: Undergraduate vaccination scenarios vector (U and D)

% Check for invalid class-caps
for ccap_idx = unique(CASES(:,1)')
    if not(any(ccap_idx==[25,50,100,50000]))
        error(['Acceptable class caps: 25, 50, 100, or 50000',...
            ' (no class cap)'])
    end
end

% Check for invalid vacc entry (values outside of the interval [0,100])
if any([[CASES(:,2)' CASES(:,3)']<0 [CASES(:,2)' CASES(:,3)']>100])
    error('Percentage of initally vaccinated must be a value in [0,100]')
end

% Make a folder to store model solutions if one does not exist
if ~exist('model_sols','dir')
       mkdir('model_sols')
end

% Load parameters of interest and their ranges for Latin Hypercube Sampling
parameters_info = fun_model_parameter_ranges;

% Load contact matrices and number of students by classification
[ConMat,Con_livingMat,...
    nUVec,nDVec,nGVec,nFVec] = fun_initialize_contact_matrices(reduced_infect);

% Model solution length in days
DaysRefined = 30; % Number of days model is solved daily
MaxDays = 7*15; % End of semester (15 week semester)

% Refine solutions during the first 30 days, then every 3 days
tspan = 24*[0:DaysRefined (DaysRefined+3):3:MaxDays];
converttohrs = 24; % 24 hrs/day

% Number of times to sample the parameter space (ps)
% SubSampN (user input)

% Generate Parameter Value Matrices: A, B
A = fun_lhs_sampling_of_parameters(parameters_info,SubSampN);
B = fun_lhs_sampling_of_parameters(parameters_info,SubSampN);

% Pre-allocate space for matrices B_Aj and A_Bj. Notations indicates 
% matrix A or B from above but the jth column is replaced with that of the 
% other matrix.
B_Aj_cell = cell(size(parameters_info,1),1); 
A_Bj_cell = cell(size(parameters_info,1),1);

% Define a matrix C_i formed by all columns of B except the ith column,
% which is taken from A
for param = 1:size(parameters_info,1)
    % Only parameters whose sensitivity is being computed
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
tic
for ccap_idx = 1:size(CASES,1)
    TOC = toc;
    class_cap = CASES(ccap_idx,1);
    switch class_cap
        case 25
            con_idx = 1;
        case 50
            con_idx = 2;
        case 100
            con_idx = 3;
        case 50000
            con_idx = 4;
    end
    if class_cap ~= 50000
        disp(['Class Cap: ',num2str(class_cap)])
    else
        disp('Class Cap: None')
    end
    % Percent Faculty and Graduate Students Vaccinated  
    Vu = CASES(ccap_idx,2); Vd = Vu;
    Vf = CASES(ccap_idx,3); Vg = Vf; 
    
    disp(['und vax: ',num2str(Vu),' & fac vax: ',num2str(Vf)])
    % Percent of undergraduates vaccinated
    
    % Load a particular contact scenario
    [Con,Con_living,nu,nd,ng,nf] = deal(...
                        ConMat(:,:,con_idx),...
                            Con_livingMat,...
                                nUVec(con_idx),...
                                    nDVec(con_idx),...
                                        nGVec(con_idx),...
                                            nFVec(con_idx)...
                                                );                                              
    % Compute YA
    YA_cell = fun_solve_covid_model(A,tspan,converttohrs,Con,Con_living,...
                            nu,nd,ng,nf,Vu,Vd,Vg,Vf);
    
    % Compute YB
    YB_cell = fun_solve_covid_model(B,tspan,converttohrs,Con,Con_living,...
                            nu,nd,ng,nf,Vu,Vd,Vg,Vf);
    
    
    % Compute YA_Bj
    YA_Bj_cell = cell(1,numel(parameters_info(:,1)));
    for param = find(parameters_info(:,1))'
        A_Bj_tmp = A_Bj_cell{param};
        YA_Bj_tmp = fun_solve_covid_model(...
                                A_Bj_tmp,tspan,converttohrs,...
                                    Con,Con_living,nu,nd,ng,nf,...
                                        Vu,Vd,Vg,Vf...
                                            );
        YA_Bj_cell{param} = YA_Bj_tmp;
    end
    
    % Save the solutions using case specific information
    save(... % File name:
        ['./model_sols/model_sols_dailyMCn',...
                num2str(SubSampN),...
                    'Vu',num2str(Vu),'Vd',num2str(Vd),...
                        'Vg',num2str(Vg),'Vf',num2str(Vf),...
                            'CC_',num2str(class_cap),'.mat'],...
                                'YA_cell','YB_cell','YA_Bj_cell',...
        ... % Variables:
            'parameters_info',...
                    'Vu','Vd','Vg','Vf',...
                        'class_cap',...
                            '-v7.3');
    
    disp(['und vax: ',num2str(Vu),...
    ' & fac vax: ',num2str(Vf),...
     ' & elapsed Time: ',num2str(round(toc/60,2)),' mins',...
     ' & simulation Time: ',num2str(round((toc-TOC)/60,2)),' mins'])
end
