function [Con,Con_living,nU,nD,nG,nF] ...
                = fun_initialize_contact_matrices(reduced_infect)
% % FUNCTION: FUN_INITIALIZE_CONTACT_MATRICES
% % AUTHOR: Fabian Santiago
% % EMAIL: fsantiago3@ucmerced.edu
% % DATE: 8/27/2021
% % DESCRIPTION: Parse contact matrix informatiom in contact_matrices folder. 
% % ABBREVIATIONS: 
% %         U: undergraduates living off campus
% %         D: undergraduates living in the doorms
% %         G: graduate students
% %         F: faculty and staff
 
% Structure of contact matrices:
%  	C[1,1]*: Undergraduates to Undergraduates 
%   C[1,2]*: Undergraduates to Dorm-Living Students
%  	C[1,3]: Undergraduates to Graduates 
%  	C[1,4]: Undergraduates to Faculty
% 
%  	C[2,1]*: Dorm-Living to Dorm-Living 
%  	C[2,2]*: Dorm-Living Students to Dorm-Living Students 
%  	C[2,3]: Dorm-Living Students to Graduates 
%  	C[2,4]: Dorm-Living Students to Faculty 
%  	
%  	C[3,1]: Graduates to Undergraduates 
%  	C[3,2]: Graduates to Dorm-Living Students 
%  	C[3,3]: Graduate to Graduates 
%  	C[3,4]: Graduate to Faculty 
%  	
%  	C[4,1]: Faculty to Undergraduates 
%  	C[4,2]: Faculty to Dorm-Living Students 
%  	C[4,3]: Faculty to Graduates 
%  	C[4,4]: Faculty to Faculty 
%
%   The contact matrix is being read in from UCM data processed by Erica
%   Rutter. 

% Error detection: 
if reduced_infect~=25 && reduced_infect~=10
    error(['ERROR: in fun_initialize_contact_matrices, ',...
            'input reduced_infect must be 25 (25%) or 10 (10%).']);
end

% Number of class caps considered
num_class_caps = 4; 

% Initialize Contact matrix
Con = zeros(4,5,num_class_caps);

% Pre-allocate space for sub-populations
nU = zeros(1,num_class_caps); nD = zeros(1,num_class_caps);
nG = zeros(1,num_class_caps); nF = zeros(1,num_class_caps);

% Hour conversion
converttohrs = 24; 

% Load Classroom Contact Scenarios
% Class caps: 25, 50, 100 students, and no-class cap (label = 50000)
class_cap = [25 50 100 50000];
for contact_idx = 1:4
    % Load contact matrices
    load(['./contact_matrices/scaled0.',num2str(reduced_infect),...
            '_reciprocal_sens_avg_contact_matrix_classcap',...
            num2str(class_cap(contact_idx)),'.mat'],'contact_matrix');
    
    % Add population information
    nU(contact_idx) = contact_matrix(1,5); 
    nD(contact_idx) = contact_matrix(2,5);
    nG(contact_idx) = contact_matrix(3,5); 
    nF(contact_idx) = contact_matrix(4,5);
    
    % Add contact matrices to contact matrix
    Con(:,:,contact_idx) = contact_matrix./converttohrs; 
end

% Load living contact scenarios
load('contact_matrices/reciprocal_sens_living_contact_matrix.mat',...
                                                        'living_matrix');
Con_living = living_matrix./converttohrs;