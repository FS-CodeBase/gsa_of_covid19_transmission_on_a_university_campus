% FUNCTION: covid_model_ode
% AUTHOR: Fabian Santiago and Shilpa Khatri
% EMAIL: fabiansantiago707@gmail.com
% DATE: 8/27/2021
function dy = covid_model_ode(t,y,ConMat,SocCon,nU,nD,nG,nF,PARAMS)
    % CONTACT MULTIPLIERS
    %%%%%%%%%%%%%%%%%%%%%
    %     livingM = PARAMS(1); dormM = PARAMS(2); communityM = PARAMS(3); 
    pM = PARAMS(4); pA = PARAMS(5);
    % SOCIAL CONTACT
    %%%%%%%%%%%%%%%%
    %     frac_social = PARAMS(6); partysize = PARAMS(7); 
    weekend_multiplier = PARAMS(8);
    % INITIAL CONDITIONS
    %%%%%%%%%%%%%%%%%%%%
    %     I0Us = PARAMS(9); I0Ua = PARAMS(10); I0Ds = PARAMS(11); I0Da = PARAMS(12); 
    % INFECTION PARAMETERS
    %%%%%%%%%%%%%%%%%%%%%%
    R0 = PARAMS(13); maskeffect = PARAMS(14); 
    %beta --> 2614.552209348986*beta = R0
    beta = R0/2614.552209348986;

    aS = PARAMS(15);
    % SIGMA -> 1/(Expected Time in Exposed State)
    sigmau = PARAMS(16); sigmad = PARAMS(17);
    sigmag = PARAMS(18); sigmaf = PARAMS(19);
    % PHI -> PROB of BECOMING ASYMPTIOMATIC
    phiu = PARAMS(20); phid = PARAMS(21); 
    phig = PARAMS(22); phif = PARAMS(23); 
    % ALPHA -> Prob of going to Hospital
    % gammaa: Prob of symptomatic who go to Hospital HS
    gammaau = PARAMS(24); gammaad = PARAMS(25); 
    gammaag = PARAMS(26); gammaaf = PARAMS(27);
    % alpha: Prob of symptomatic who go to HS 
    alphau = PARAMS(28); alphad = PARAMS(29); 
    alphag = PARAMS(30); alphaf = PARAMS(31);
    % gammas: 1/(Expected Time in Symptomatic State BEFORE you go to health services)
    gammasu = PARAMS(32); gammasd = PARAMS(33); 
    gammasg = PARAMS(34); gammasf = PARAMS(35);
    % hi -> 1/Expected Time in the health services/Quarantine
    hu = PARAMS(36); hd = PARAMS(37); hg = PARAMS(38); hf = PARAMS(39);
    % deltai -> 1/(Expected Time in Symptomatic State IF you do NOT go to health services)
    deltau = PARAMS(40); deltad = PARAMS(41); 
    deltag = PARAMS(42); deltaf = PARAMS(43);
    
    % Weekend Social Contacts
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Weekends are assumed to go from Friday at 5pm to Sunday at 5pm 
    if (rem(t,168)<=168 && rem(t,168)>=120)
        % party on the weekend
        SocCon = weekend_multiplier * SocCon;
    end
    ConMat = ConMat + SocCon/maskeffect; % not wearing masks at parties;

% Key for ODE below 
% 1: Su
% 2: Sd
% 3: Sg
% 4: Sf 
% 5: Eu
% 6: Ed
% 7: Eg
% 8: Ef 
% 9: Iau 
% 10: Isu
% 11: Iad 
% 12: Isd 
% 13: Iag 
% 14: Isg 
% 15: Iaf 
% 16: Isf 
% 17: CI 
% 18: Hu
% 19: Hd
% 20: Hg
% 21: Hf
% 22: CH 
% 23: Nu
% 24: Nd
% 25: Ng 
% 26: Nf
% 27: R 
% 28: Cu %Cumulative Infections Undergrad Off Campus
% 29: Cd %Cumulative Infections Dorms
% 30: Cg %Cumulative Infections Grad Students
% 31: Cf %Cumulative Infections Faculty
 
%Susceptible: Undergrads, Dorm Students, Grads and Faculty 
dy(1) = -maskeffect * beta * y(1) * (ConMat(1,1) * (aS * y(9) + y(10) + y(23))/nU + ConMat(1,2) * (aS * y(11) + y(12) + y(24))/nD + ConMat(1,3) * (aS * y(13) + y(14) + y(25))/nG + ConMat(1,4)* (aS * y(15) + y(16) + y(26))/nF + ConMat(1,5) * (aS * pA * pM + (1 - pA) * pM) );
dy(2) = -maskeffect * beta * y(2) * (ConMat(2,1) * (aS * y(9) + y(10) + y(23))/nU + ConMat(2,2) * (aS * y(11) + y(12) + y(24))/nD + ConMat(2,3) * (aS * y(13) + y(14) + y(25))/nG + ConMat(2,4)* (aS * y(15) + y(16) + y(26))/nF + ConMat(2,5) * (aS * pA * pM + (1 - pA) * pM) );
dy(3) = -maskeffect * beta * y(3) * (ConMat(3,1) * (aS * y(9) + y(10) + y(23))/nU + ConMat(3,2) * (aS * y(11) + y(12) + y(24))/nD + ConMat(3,3) * (aS * y(13) + y(14) + y(25))/nG + ConMat(3,4)* (aS * y(15) + y(16) + y(26))/nF + ConMat(3,5) * (aS * pA * pM + (1 - pA) * pM) ); 
dy(4) = -maskeffect * beta * y(4) * (ConMat(4,1) * (aS * y(9) + y(10) + y(23))/nU + ConMat(4,2) * (aS * y(11) + y(12) + y(24))/nD + ConMat(4,3) * (aS * y(13) + y(14) + y(25))/nG + ConMat(4,4)* (aS * y(15) + y(16) + y(26))/nF + ConMat(4,5) * (aS * pA * pM + (1 - pA) * pM) ); 
%Exposed : Undergrads, Dorm Students, Grads and Faculty
dy(5) =  maskeffect * beta * y(1) * (ConMat(1,1) * (aS * y(9) + y(10) + y(23))/nU + ConMat(1,2) * (aS * y(11) + y(12) + y(24))/nD + ConMat(1,3) * (aS * y(13) + y(14) + y(25))/nG + ConMat(1,4)* (aS * y(15) + y(16) + y(26))/nF + ConMat(1,5) * (aS * pA * pM + (1 - pA) * pM) ) - sigmau * y(5); 
dy(6) =  maskeffect * beta * y(2) * (ConMat(2,1) * (aS * y(9) + y(10) + y(23))/nU + ConMat(2,2) * (aS * y(11) + y(12) + y(24))/nD + ConMat(2,3) * (aS * y(13) + y(14) + y(25))/nG + ConMat(2,4)* (aS * y(15) + y(16) + y(26))/nF + ConMat(2,5) * (aS * pA * pM + (1 - pA) * pM) ) - sigmad * y(6); 
dy(7) =  maskeffect * beta * y(3) * (ConMat(3,1) * (aS * y(9) + y(10) + y(23))/nU + ConMat(3,2) * (aS * y(11) + y(12) + y(24))/nD + ConMat(3,3) * (aS * y(13) + y(14) + y(25))/nG + ConMat(3,4)* (aS * y(15) + y(16) + y(26))/nF + ConMat(3,5) * (aS * pA * pM + (1 - pA) * pM) ) - sigmag * y(7);
dy(8) =  maskeffect * beta * y(4) * (ConMat(4,1) * (aS * y(9) + y(10) + y(23))/nU + ConMat(4,2) * (aS * y(11) + y(12) + y(24))/nD + ConMat(4,3) * (aS * y(13) + y(14) + y(25))/nG + ConMat(4,4)* (aS * y(15) + y(16) + y(26))/nF + ConMat(4,5) * (aS * pA * pM + (1 - pA) * pM) ) - sigmaf * y(8);
%Infected Undergrads
dy(9) = phiu * sigmau * y(5) - gammaau * y(9); 
dy(10) = (1 - phiu)* sigmau * y(5) - gammasu * y(10); 
%Infected Dorm Students
dy(11) = phid * sigmad * y(6) - gammaad * y(11); 
dy(12) = (1 - phid)* sigmad * y(6) - gammasd * y(12); 
%Infected Grads
dy(13) = phig * sigmag * y(7) - gammaag * y(13); 
dy(14) = (1 - phig)* sigmag * y(7) - gammasg * y(14); 
%Infected Faculty 
dy(15) = phif * sigmaf * y(8) - gammaaf * y(15); 
dy(16) = (1 - phif)* sigmaf * y(8) - gammasf * y(16);
%Cumulative Infections
dy(17) = sigmau * y(5) + sigmad * y(6) + sigmag * y(7) + sigmaf * y(8);
%In health services
dy(18) = alphau * gammasu * y(10) - hu * y(18);
dy(19) = alphad * gammasd * y(12) - hd * y(19);
dy(20) = alphag * gammasg * y(14) - hg * y(20);
dy(21) = alphaf * gammasf * y(16) - hf * y(21);
%Cumulative at health services 
dy(22) = alphau * gammasu * y(10) + alphad * gammasd * y(12) + alphag * gammasg * y(14) + alphaf * gammasf * y(16); 
%Not quarantine / health services
dy(23) = (1 - alphau) * gammasu * y(10)  - deltau * y(23);
dy(24) = (1 - alphad) * gammasd * y(12)  - deltad * y(24);
dy(25) = (1 - alphag) * gammasg * y(14)  - deltag * y(25);
dy(26) = (1 - alphaf) * gammasf * y(16)  - deltaf * y(26);
%Recovered 
dy(27) = gammaau * y(9) + gammaad * y(11) + gammaag * y(13) + gammaaf * y(15) ...
    + hu * y(18) + hd * y(19) + hg * y(20) + hf * y(21)...
    + deltau * y(23) + deltad * y(24) + deltag * y(25) + deltaf * y(26); 
%Cumulative Undergrad
dy(28) = sigmau * y(5);
%Cumulative Dorms
dy(29) = sigmad * y(6);
%Cumulative Grad Students
dy(30) = sigmag * y(7);
%Cumulative Faculty/Staff;
dy(31) = sigmaf * y(8);

dy = dy';
end