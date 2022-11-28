function [parameters_info, prms_str] = fun_model_parameter_ranges
% % FUNCTION: MODEL_PARAMETER_RANGES
% % AUTHOR: Fabian Santiago
% % EMAIL: fsantiago3@ucmerced.edu
% % UPDATE: 10/27/2022
% % DESCRIPTION: Parameter values and their ranges. See Table 1 in main
% %              manuscript. Each ROW represents ONE parameter and has the
% %              following form: [flag1 lower_bound upper_bound flag2]
% %              flag1 := Compute sensitivity of parameter (true/false)
% %              lower_bound := Parameter lower bound
% %              upper_bound := Parameter upper bound
% %              flag2 := Convert from days to hourly rate (true/false)
% %        NOTE: If flag1 is set to false, then the lower_bound value is 
% %              used as the model parameter value.
% % 
% % ABBREVIATIONS: 
% %         U: undergraduates living off campus
% %         D: undergraduates living in the doorms
% %         G: graduate students
% %         F: faculty and staff

parameters_info = [... % [SA_Bool a  b Days_Bool]
% CONTACT MULTIPLIERS: multiplies contacts with:  
   false 1.00  1.00    false;... % [~] U
   false 1.00  1.00    false;... % [~] D
   true  1.00  10.00   false;... % [c] Community
% PERCENTAGE OF MERCED COMMUNITY INFECTIONS: (M = pS + pA)
   true  0.00  0.005   false;... % [pS] Symptomatic
   false 0.01  0.005   false;... % [pA] Asymptomatic
% SOCIAL CONTACT PARAMETERS
   false 0.75  0.75    false;... % [~] Frac of time socializing
   true  0.50  1.50    false;... % [p] Social percentage
   true  1.00  10.00   false;... % [w] Weekend multiplier
% INFECTION PARAMETERS
   true  1.50  6.50    false;... % [beta] Transmission rate
   true  0.40  1.00    false;... % [m] 1-m reduction in beta by wearing masks
   false 0.50  0.50    false;... % [~] Infection transmission probability of asymptomatics
% Expected Time in Exposed State (varied once and assumed equal for each subpopulation)
   true  1.00  7.50    true;...  % [sigma] U
   false 1.00  7.50    true;...  % [sigma] D
   false 1.00  7.50    true;...  % [sigma] G
   false 1.00  7.50    true;...  % [sigma] F
% Probability of Becoming Asymptomatic
   true  0.20  0.80    false;... % [phi] U
   false 0.50  0.50    false;... % [phi] U
   false 0.50  0.50    false;... % [phi] G
   false 0.50  0.50    false;... % [phi] F
% Expected time in the symptomatic state 
% before voluntary self-isolation
   false 2.00  2.00    true;...  % [~] U
   false 2.00  2.00    true;...  % [~] D
   false 2.00  2.00    true;...  % [~] G
   false 2.00  2.00    true;...  % [~] F
% Expected duration in the asymptomatic infectious state
   false 14.00 14.00   true;...  % [~] U
   false 14.00 14.00   true;...  % [~] D
   false 14.00 14.00   true;...  % [~] G
   false 14.00 14.00   true;...  % [~] F
% Probability of a symptomatic individual will self-isolate 
   true  0.00  1.00    false;... % [alpha] U
   false 0.50  1.00    false;... % [alpha] D
   false 0.50  1.00    false;... % [alpha] G
   false 0.50  1.00    false;... % [alpha] F
% INITIAL CONDITIONS/INITIALLY INFECTED 
   true  0.00  0.50    false;... % [v_{us,0}] symptomatic infected U
   true  0.00  0.50    false;... % [v_{ua,0}] asymptomatic infected U
   true  0.00  0.50    false;... % [v_{da,0}] symptomatic infected D
   true  0.00  0.50    false;... % [v_{ds,0}] asymptomatic infected D
% Graduate students
   true  0.00  0.50    false;... % [v_{gs,0}]  symptomatic infected G
   true  0.00  0.50    false;... % [v_{ga,0}]  asymptomatic infected G
% Faculty and staff
   true  0.00  0.50    false;... % [v_{fs,0}]  symptomatic infected F
   true  0.00  0.50    false;... % [v_{fa,0}]  asymptomatic infected F
% IMMUNITY DURATION (ADDITIONAL INFECTION PARAMETERS)
   true  90    270     true;...  % [tR] Length of natural immunity
   true  180   270     true;...  % [tV] Length of vaccine immunity
];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% PARAMETER LABELS FOR FIGURES IN PAPER %%%%%%%%%%%%%%%%%%%%
% The notation (a:b) indicate the index number for their position in the 
% prms_info object defined above.
prms_str = {... % CONTACT MULTIPLIERS 
            '$l_m$','$d_m$','$c$',... (1:3) [living, dorm, community]
            ... % COMMUNITY INFECTIONS 
            '$M$','$M$',... (4:5) Fraction: [Symp (M) and Asymp]
            ... % SOCIAL CONTACT PARAMETERS
             '$f_s$','$p$','$w$',... (6:8) [social,party size,weekend x]
            ... % INFECTION PARAMETERS 
             '$\beta$','$m$',... % (9:10) [transmission, mask usage]
             '$aS$',...% (11) Probability of asymptomatic infectiousness
             '$\sigma$','$\sigma$','$\sigma$','$\sigma$',... % (12:15) SIGMA  = 1/(Expected Time in Exposed State)
             '$\phi$','$\phi$','$\phi$','$\phi$',... % (16:19) PROB of BECOMING ASYMPTIOMATIC
             '$\gamma_{su}$','$\gamma_{sd}$','$\gamma_{sg}$','$\gamma_{sf}$',... % (20:23) GAMMA: 1/(Expected Time in Symptomatic State BEFORE you go to health services)
             '$\gamma_{au}$','$\gamma_{ad}$','$\gamma_{ag}$','$\gamma_{af}$',... % (32:35) GAMMA: 1/(Expected Duration in Asymptomatic Infectious State)
             '$\alpha$','$\alpha$','$\alpha$','$\alpha$',... %(28:31) ALPHA: Prob of symptomatic who go to Hospital HS
             ... % INITIAL CONDITIONS: UNDERGRADUATE STUDENTS
             '$I_{u0}^s$','$I_{u0}^a$','$I_{d0}^s$','$I_{d0}^a$',... %(32:35)
             ... % INITIAL CONDITIONS: GRADUATE STUDENTS
             '$I_{g0}^s$','$I_{g0}^a$',... %(36:37)
             ... % INITIAL CONDITIONS: FACULTY
             '$I_{f0}^s$','$I_{f0}^a$',... %(38:39)
             ... % IMMUNITY DURATION PARAMETERS (natural,vax) (40:41)
             '$t_R$','$t_V$',... % Natural t_R and Vaccination t_V
             }; 