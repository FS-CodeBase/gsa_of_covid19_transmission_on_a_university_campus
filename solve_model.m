function Y_cell = solve_model(params_mat,tspan,tohrs,Con,Con_living,...
                                nU,nD,nG,nF,vU,vD,vG,vF)
% FUNCTION: solve_model
% AUTHOR: Fabian Santiago
% EMAIL: fsantiago@math.arizona.edu
% DATE: 8/27/2021
N = size(params_mat,1);
Y_cell = cell(1,N);

% Solve model for each set of MC sampled parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parfor n = 1:N
    PARAMS = params_mat(n,:);
    % CONTACT MULTIPLIERS
    livingM = PARAMS(1); dormM = PARAMS(2); communityM = PARAMS(3);
    % SOCIAL CONTACT
    frac_social = PARAMS(6); partysize = PARAMS(7);
    % INFECTION PARAMETERS
%     R0 = PARAMS(13); 
    maskeffect = PARAMS(14); 
    % INITIAL CONDITIONS: Undergraduates
    I0Us = PARAMS(9); I0Ua = PARAMS(10); 
    I0Ds = PARAMS(11); I0Da = PARAMS(12); 
    % INITIAL CONDITIONS: Graduate + Faculty
    I0Gs = PARAMS(36); I0Ga = PARAMS(37); 
    I0Fs = PARAMS(38); I0Fa = PARAMS(39); 
    
    % ALPHA -> Prob of going to Hospital
    % gammaa: Prob of symptomatic who go to Hospital HS
    gammaau = PARAMS(24); gammaad = PARAMS(25); 
    gammaag = PARAMS(26); gammaaf = PARAMS(27);
%     % alpha: Prob of symptomatic who go to HS 
% 	alphau = PARAMS(28); alphad = PARAMS(29); 
%   alphag = PARAMS(30); alphaf = PARAMS(31);
    % gammas: 1/(Expected Time in Symptomatic State 
    %               BEFORE you go to health services)
    gammasu = PARAMS(32); gammasd = PARAMS(33); 
    gammasg = PARAMS(34); gammasf = PARAMS(35);
    
    % hi -> 1/Expected Time in the health services/Quarantine
    hu = 1/(1/gammaau - 1/gammasu);    hd = 1/(1/gammaad - 1/gammasd);
    hg = 1/(1/gammaag - 1/gammasg);    hf = 1/(1/gammaaf - 1/gammasf);
    % deltai -> 1/(Expected Time in Symptomatic State 
    %                 IF you do NOT go to health services)
    deltau = hu; deltad = hd; deltag = hg; deltaf = hf;
    PARAMS(36:43) = [hu,hd,hg,hf,deltau,deltad,deltag,deltaf];
    
    % Update CONTACT MATRICES
    Con_living_tmp = livingM*Con_living./maskeffect; %(MODIFIED)*
    % Modify the living contact rate for dorm to dorm 
    % students by multiplier dormM 
    Con_living_tmp(2,2) = dormM*Con_living_tmp(2,2);
    % Add living matrix to the course contact matrices
    Con_tmp = Con + Con_living_tmp;
    
    % Community contacts added into the matrices note communityM 
    % multiplies this effect assumed masks are not worn in these contacts
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Undergrad to merced community
    Con_tmp(1,5) = (communityM*5/tohrs)/maskeffect;
    % Dorm to merced community 
    Con_tmp(2,5) = (communityM*1/tohrs)/maskeffect;	
    % Grad to merced community 
    Con_tmp(3,5) = (communityM*5/tohrs)/maskeffect;	
    % Faculty to merced community
    Con_tmp(4,5) = (communityM*15/tohrs)/maskeffect;	 
    
%     % For cases 23 - 35 remove graduate students and faculty/staff
%     if ConScenario >= 23 && ConScenario <= 35
%         Con_tmp(1,3:5) = 0; Con_tmp(2,3:5) = 0;
%         Con_tmp(3,:)   = 0; Con_tmp(4,:)   = 0;
%     end
    % Social contacts
    % partysize: daily contact hours for daily unscheduled contact
    %            (varied from 1hr to 10hrs)
    SocCon = zeros(size(Con_tmp));
    SocCon(1,1) = partysize * (frac_social/tohrs);
    SocCon(2,2) = partysize * (frac_social/tohrs);
    SocCon(1,2) = (partysize * (1-frac_social))/tohrs;
    SocCon(2,1) = (partysize * (1-frac_social))/tohrs;
    
    % EVALUATE MODEL
    %Initial conditions 
    y0     = zeros(31,1); 
    y0(1)  = (nU - I0Ua - I0Us)*(1-vU/100); %Su 
    y0(2)  = (nD - I0Da - I0Ds)*(1-vD/100); %Sd 
% %     y0(1)  = (nU)*(1-vU/100) - I0Us - I0Ua; %Su 
% %     y0(2)  = (nD)*(1-vD/100) - I0Ds - I0Da; %Sd 
    y0(3)  = (nG - I0Ga - I0Gs)*(1-vG/100); %Sg
    y0(4)  = (nF - I0Fa - I0Fs)*(1-vF/100); %Sf 
    y0(9)  = I0Ua; %Iau (asym-Und-off-campus)
    y0(10) = I0Us; %Isu (sym-Und-off-campus)
    y0(11) = I0Da; %Iad (asym-Und-on-campus)
    y0(12) = I0Ds; %Isd (sym-Und-on-campus)
    y0(13) = I0Ga; %Iag (asym-Grad)
    y0(14) = I0Gs; %Isg (sym-Grad)
    y0(15) = I0Fa; %Iaf (asym-Fac)
    y0(16) = I0Fs; %Isf (sym-Fac)
    y0(17) = I0Ua + I0Us + I0Da + I0Ds + I0Ga + I0Gs + I0Fa + I0Fs;   % CI
    y0(28) = I0Ua + I0Us;
    y0(29) = I0Da + I0Ds;
    y0(30) = I0Ga + I0Gs;
    y0(31) = I0Fa + I0Fs;
    y0(27) = (nU - I0Ua - I0Us)*(vU/100)+(nD - I0Da - I0Ds)*(vU/100)+...
             (nG - I0Ga - I0Gs)*(vG/100)+(nF - I0Fa - I0Fs)*(vF/100);
    % Solve covid model with ode45
    [~,Ysols] = ode45(@(t,y) covid_model_ode(t,y,Con_tmp,SocCon,...
                                                nU,nD,nG,nF,PARAMS),...
                        tspan,y0);
%   SAVE: [S E I R C]
%     Y_cell{n} = [sum(Ysols(:,1:4),2)  sum(Ysols(:,5:8),2) 
%                  sum(Ysols(:,9:16),2) Ysols(:,27) Ysols(:,17)];
    Y_cell{n} = Ysols(:,17); % Only save C (cumulative cases)
end