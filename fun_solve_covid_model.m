function Y_cell = fun_solve_covid_model(params_mat,tspan,tohrs,Con,Con_living,...
                                nU,nD,nG,nF,vU,vD,vG,vF)
% % FUNCTION: FUN_SOLVE_COVID_MODEL
% % AUTHOR: Fabian Santiago
% % EMAIL: fsantiago3@ucmerced.edu
% % DATE: 10/27/2022
% % 
% % DESCRIPTION: Solve covid_model_ode for each of the N parameter sets.
% % 
% % ABBREVIATIONS: 
% %         nU: number of undergraduates living off campus
% %         nD: number of undergraduates living in the doorms
% %         nG: number of graduate students
% %         nF: number of faculty and staff

% Determine number of parameters and allocate space to store solutions.
N = size(params_mat,1);
Y_cell = cell(1,N);

% Solve model for each set of LHS sampled parameters
parfor n = 1:N
    PARAMS = params_mat(n,:);

    % CONTACT MULTIPLIERS
    livingM = PARAMS(1); dormM = PARAMS(2); communityM = PARAMS(3);

    % SOCIAL CONTACT
    p_soc_con_x  = PARAMS(7);

    % INFECTION PARAMETERS
    maskeffect = PARAMS(10); 
    
    % DURATION IN ASY/SYMPTOMATIC INFECTIOUS STATE
    % (SYMPTOMATIC)
    gammasu = PARAMS(20); gammasd = PARAMS(21); 
    gammasg = PARAMS(22); gammasf = PARAMS(23);
    % (ASYMPTOMATIC)
    gammaau = PARAMS(24); gammaad = PARAMS(25); 
    gammaag = PARAMS(26); gammaaf = PARAMS(27);

    % INITIAL CONDITIONS: Undergraduates
    pI0Us = PARAMS(32); pI0Ua = PARAMS(33); 
    pI0Ds = PARAMS(34); pI0Da = PARAMS(35); 
    
    % INITIAL CONDITIONS: Graduate + Faculty
    pI0Gs = PARAMS(36); pI0Ga = PARAMS(37); 
    pI0Fs = PARAMS(38); pI0Fa = PARAMS(39);

    % EXPECTED TIME IN SELF-ISOLATION
    hu = 1/(1/gammaau - 1/gammasu);    hd = 1/(1/gammaad - 1/gammasd);
    hg = 1/(1/gammaag - 1/gammasg);    hf = 1/(1/gammaaf - 1/gammasf);
    
    % EXPECTED TIME IN SYMPTOMATIC STATE WITHOUT SELF-ISOLATION
    deltau = hu; deltad = hd; deltag = hg; deltaf = hf;

    % SAVE PARAMETERS FOR ODE MODEL
    PARAMS(42:49) = [hu,hd,hg,hf,deltau,deltad,deltag,deltaf];
    
    % SET CONTACT MATRICES
    Con_living_tmp = livingM*Con_living./maskeffect;

    % APPLY DORM CONTACT MULTIPLIER
    Con_living_tmp(2,2) = dormM*Con_living_tmp(2,2);

    % UPDATE CONTACT MATRIX WITH LIVING CONTACTS
    Con_tmp = Con + Con_living_tmp;
    
    % UPDATE COMMUNITY CONTACTS: This also accounts for the fact that we
    % assume students do not use masks off campus.

    % Undergrad to merced community
    Con_tmp(1,5) = (communityM*5/tohrs)/maskeffect;
    % Dorm to merced community 
    Con_tmp(2,5) = (communityM*1/tohrs)/maskeffect;	
    % Grad to merced community 
    Con_tmp(3,5) = (communityM*5/tohrs)/maskeffect;	
    % Faculty to merced community 
    Con_tmp(4,5) = (communityM*15/tohrs)/maskeffect;	 
    

    % SET SOCIAL CONTACT MATRIX
    SocCon = zeros(size(Con_tmp));
    SocCon(1:2,1:2) = p_soc_con_x*Con(1:2,1:2);

    % SET INITIAL CONDITIONS
    y0     = zeros(34,1); 
    y0(1)  = (nU*(1 - (pI0Ua + pI0Us)/100))*(1 - vU/100); % Su 
    y0(2)  = (nD*(1 - (pI0Da + pI0Ds)/100))*(1 - vD/100); % Sd 
    y0(3)  = (nG*(1 - (pI0Ga + pI0Gs)/100))*(1 - vG/100); % Sg
    y0(4)  = (nF*(1 - (pI0Fa + pI0Fs)/100))*(1 - vF/100); % Sf 
    y0(9)  = nU*(pI0Ua)/100; % Iau (asym-Und-off-campus)
    y0(10) = nU*(pI0Us)/100; % Isu (sym-Und-off-campus)
    y0(11) = nD*(pI0Da)/100; % Iad (asym-Und-on-campus)
    y0(12) = nD*(pI0Ds)/100; % Isd (sym-Und-on-campus)
    y0(13) = nG*(pI0Ga)/100; % Iag (asym-Grad)
    y0(14) = nG*(pI0Gs)/100; % Isg (sym-Grad)
    y0(15) = nF*(pI0Fa)/100; % Iaf (asym-Fac)
    y0(16) = nF*(pI0Fs)/100; % Isf (sym-Fac)

    % Cumulative infections
    y0(17) = sum(y0(9:16));

    % Vaccinated: U,D,G,F
    % 27: Vu; 28: Vd; 29: Vg; 30: Vf 
    y0(31) = (nU - y0(9)  - y0(10))*(vU/100); % U
    y0(32) = (nD - y0(11) - y0(12))*(vD/100); % D
    y0(33) = (nG - y0(13) - y0(14))*(vG/100); % G
    y0(34) = (nF - y0(15) - y0(16))*(vF/100); % F

    % Solve covid model with ode45
    [~,Ysols] = ode45(@(t,y) fun_covid_model_ode(t,y,Con_tmp,SocCon,...
                                                nU,nD,nG,nF,PARAMS),...
                                                                 tspan,y0);
% ONLY SAVE CUMULATIVE INFECTIONS  
%     Y_cell{n} = [sum(Ysols(:,1:4),2)  sum(Ysols(:,5:8),2) 
%                  sum(Ysols(:,9:16),2) Ysols(:,27) Ysols(:,17)];
%     Solutions: Susceptible = sum(Ysols(:,1:4),2)
%                Exposed     = sum(Ysols(:,5:8),2)
%                Infectious  = sum(Ysols(:,9:16),2)
%                Recovered  = Ysols(:,27)
%                Cumulative Infections = Ysols(:,17)
    Y_cell{n} = Ysols(:,17);
end