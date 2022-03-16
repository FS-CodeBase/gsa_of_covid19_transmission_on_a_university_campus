function [Con,Con_living,nU,nD,nG,nF] = fun_initialize_contact_matrices
% FUNCTION: fun_initialize_contact_matrices
% AUTHOR: Fabian Santiago
% EMAIL: fabiansantiago707@gmail.com
% DATE: 8/27/2021
% DESCRIPTION: 

% Subpopulations 
%%%%%%%%%%%%%%%%%
% Undergraduate Students: u 
% Students who live in the Dorms: d
% Graduate Students: g
% Faculty/Staff: f
% 
% Structure:
% 1) 16 Classes:  S (4 classes), E (4 classes), 
%    I(s/a) (4*2 = 8 classes), R (1 class)
%  	(a) u,d,g,f  = undergraduates, dorm-living students, graduates, faculty
%  	(b) a vs s = asymptomatic vs  symptomatic
%   (c) H vs N = Health Services vs not 
%  	
%  2) Campus Population: (https://www.ucmerced.edu/fast-facts)
%  	(a) Undergraduates: population being computed from contact data 
%  	(b) Graduates = 696
%  	(c) Faculty/Staff = 250 Ladder Rank + 157 Non Ladder Rank =  407
%  	
%  3) Contact Rates (# contacts per hour) 
%     Note we often measure hours of contact per day and divide by 24 
%  	  Infectious Contact Rate = # Hourly Contacts*Transm. Prob Per Contact
%  	
%  	C[1,1]: Undergraduates to Undergraduates 
%   C[1,2]: Undergraduates to Dorm-Living Students
%  	C[1,3]: Undergraduates to Graduates 
%  	C[1,4]: Undergraduates to Faculty
% 
%  	C[2,1]: Graduates to Undergraduates 
%  	C[2,2]: Dorm-Living Students to Dorm-Living Students 
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
%  	
%  4) Infectious Period
%  	(a) Assume Latent Period lasts 5 days;
%  	(b) Assume period of infectivity is 14 days; 
%  	(c) Average time of infectivity for unvaccinated asymptomatic? 1.9 day
%
%  
%  5) Symptomatic vs Asymptomatic 
%  	(a) Assumed 50% of infected would be asymptomatic for undergrads 
%       \[Phi]u (will have different rates)
%   (b) For now assume that symptomatic people behave the same as 
%       asymptomatic people (i.e., there's no change in behavior)
%
% Assumptions and Problems:
% - Heterogeneous information is  being averaged out. Faculty, Staff, Age
%   Distributions.
% - The model is HIGHLY sensitive towards probability of trans. [Beta] 
% - We are not modeling a change in behavior.

% Contact Matrices 
% Con = #Contacts/Hour
% Note beta will change based on which contact scenario you are considering
% so if you change the matrices - also need to change beta below 
con_files = dir('./contact_matrices/');
con_files(1:2) = [];
nCon = numel(con_files)/2; % Number of contact matrix scenarios to consider 

% Initialize Contact matrices
Con        = zeros(4,5,nCon); % Initialization Housing Contact Matrix 
Con_living = zeros(4,5,nCon); % Initialize Living Contact Matrix

%Population can change with students NOT coming to Merced.
nU = zeros(1,nCon); nD = zeros(1,nCon);
nG = zeros(1,nCon); nF = zeros(1,nCon);

% Hour conversion parameter
converttohrs = 24; 

% Load Contact Scenarios
for con = 1:nCon
    % Load contact matrices
    load(['./contact_matrices/sens_avg_contact_matrix_',...
            num2str(con),'.mat'],'contact_matrix');
    load(['./contact_matrices/sens_living_contact_matrix_',...
            num2str(con),'.mat'],'living_matrix');
    
    % Add population information
    nU(con) = contact_matrix(1,5); nD(con) = contact_matrix(2,5);
    nG(con) = contact_matrix(3,5); nF(con) = contact_matrix(4,5);
    
    % Add contact matrices to contact matrix
    Con(:,:,con) = contact_matrix./converttohrs;
    Con_living(:,:,con) = living_matrix./converttohrs;
end