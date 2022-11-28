function dy = fun_covid_model_ode(t,y,ConMat,SocCon,nU,nD,nG,nF,PARAMS)
% % FUNCTION: FUN_COVID_MODEL_ODE
% % AUTHOR: Shilpa Khatri and Fabian Santiago
% % EMAIL: fsantiago3@ucmerced.edu
% % DATE: 8/27/2021
% % DESCRIPTION: Model of covid-19 transmission on a structured 
% %              college campus.

% COMMUNITY INFECTIONS
pM = PARAMS(4); pA = PARAMS(5);

% SOCIAL CONTACT
weekend_multiplier = PARAMS(8);
maskeffect = PARAMS(10); 

% INFECTION PARAMETERS
beta = PARAMS(9)/2614.552209348986;

% PROBABILITY OF TRANSMISSION BY ASYMPTOMATICS
aS = PARAMS(11);
sigmau = PARAMS(12); sigmad = PARAMS(13);
sigmag = PARAMS(14); sigmaf = PARAMS(15);

% PHI
phiu = PARAMS(16); phid = PARAMS(17); 
phig = PARAMS(18); phif = PARAMS(19); 

% GAMMA
gammasu = PARAMS(20); gammasd = PARAMS(21); 
gammasg = PARAMS(22); gammasf = PARAMS(23);
gammaau = PARAMS(24); gammaad = PARAMS(25); 
gammaag = PARAMS(26); gammaaf = PARAMS(27);

% ALPHA
alphau = PARAMS(28); alphad = PARAMS(29); 
alphag = PARAMS(30); alphaf = PARAMS(31);

% IMMUNITY
tR = PARAMS(40); % Natural immunity
tV = PARAMS(41); % Vaccine immunity

% EXPECTED TIME IN SELF-ISOLATION
hu = PARAMS(42); hd = PARAMS(43); hg = PARAMS(44); hf = PARAMS(45);

% EXPECTED TIME IN SYMPTOMATIC STATE IF w/o SELF-ISOLATION
deltau = PARAMS(46); deltad = PARAMS(47); 
deltag = PARAMS(48); deltaf = PARAMS(49);

% WEEKEND SOCIAL CONTACTS
% Weekends are assumed to go from Friday at 5pm to Sunday at 5pm 
if (rem(t,168)<=168 && rem(t,168)>=120)
    SocCon = weekend_multiplier * SocCon;
end
ConMat = ConMat + SocCon/maskeffect;

% MODEL COMPARTMENTS
% Susceptible: Undergrads, Dorm Students, Grads and Faculty 
dy(1) = tV*y(31)+tR*y(27)-maskeffect * beta * y(1) * (ConMat(1,1) * (aS * y(9) + y(10) + y(23))/nU + ConMat(1,2) * (aS * y(11) + y(12) + y(24))/nD + ConMat(1,3) * (aS * y(13) + y(14) + y(25))/nG + ConMat(1,4)* (aS * y(15) + y(16) + y(26))/nF + ConMat(1,5) * (aS * pA * pM + (1 - pA) * pM) );
dy(2) = tV*y(32)+tR*y(28)-maskeffect * beta * y(2) * (ConMat(2,1) * (aS * y(9) + y(10) + y(23))/nU + ConMat(2,2) * (aS * y(11) + y(12) + y(24))/nD + ConMat(2,3) * (aS * y(13) + y(14) + y(25))/nG + ConMat(2,4)* (aS * y(15) + y(16) + y(26))/nF + ConMat(2,5) * (aS * pA * pM + (1 - pA) * pM) );
dy(3) = tV*y(33)+tR*y(29)-maskeffect * beta * y(3) * (ConMat(3,1) * (aS * y(9) + y(10) + y(23))/nU + ConMat(3,2) * (aS * y(11) + y(12) + y(24))/nD + ConMat(3,3) * (aS * y(13) + y(14) + y(25))/nG + ConMat(3,4)* (aS * y(15) + y(16) + y(26))/nF + ConMat(3,5) * (aS * pA * pM + (1 - pA) * pM) ); 
dy(4) = tV*y(34)+tR*y(30)-maskeffect * beta * y(4) * (ConMat(4,1) * (aS * y(9) + y(10) + y(23))/nU + ConMat(4,2) * (aS * y(11) + y(12) + y(24))/nD + ConMat(4,3) * (aS * y(13) + y(14) + y(25))/nG + ConMat(4,4)* (aS * y(15) + y(16) + y(26))/nF + ConMat(4,5) * (aS * pA * pM + (1 - pA) * pM) ); 

% Exposed : Undergrads, Dorm Students, Grads and Faculty
dy(5) =  maskeffect * beta * y(1) * (ConMat(1,1) * (aS * y(9) + y(10) + y(23))/nU + ConMat(1,2) * (aS * y(11) + y(12) + y(24))/nD + ConMat(1,3) * (aS * y(13) + y(14) + y(25))/nG + ConMat(1,4)* (aS * y(15) + y(16) + y(26))/nF + ConMat(1,5) * (aS * pA * pM + (1 - pA) * pM) ) - sigmau * y(5); 
dy(6) =  maskeffect * beta * y(2) * (ConMat(2,1) * (aS * y(9) + y(10) + y(23))/nU + ConMat(2,2) * (aS * y(11) + y(12) + y(24))/nD + ConMat(2,3) * (aS * y(13) + y(14) + y(25))/nG + ConMat(2,4)* (aS * y(15) + y(16) + y(26))/nF + ConMat(2,5) * (aS * pA * pM + (1 - pA) * pM) ) - sigmad * y(6); 
dy(7) =  maskeffect * beta * y(3) * (ConMat(3,1) * (aS * y(9) + y(10) + y(23))/nU + ConMat(3,2) * (aS * y(11) + y(12) + y(24))/nD + ConMat(3,3) * (aS * y(13) + y(14) + y(25))/nG + ConMat(3,4)* (aS * y(15) + y(16) + y(26))/nF + ConMat(3,5) * (aS * pA * pM + (1 - pA) * pM) ) - sigmag * y(7);
dy(8) =  maskeffect * beta * y(4) * (ConMat(4,1) * (aS * y(9) + y(10) + y(23))/nU + ConMat(4,2) * (aS * y(11) + y(12) + y(24))/nD + ConMat(4,3) * (aS * y(13) + y(14) + y(25))/nG + ConMat(4,4)* (aS * y(15) + y(16) + y(26))/nF + ConMat(4,5) * (aS * pA * pM + (1 - pA) * pM) ) - sigmaf * y(8);

% Infected Undergrads
dy(9) = phiu * sigmau * y(5) - gammaau * y(9); % Asymptomatic
dy(10) = (1 - phiu)* sigmau * y(5) - gammasu * y(10); % Symptomatic

% Infected Dorm Students
dy(11) = phid * sigmad * y(6) - gammaad * y(11); % Asymptomatic 
dy(12) = (1 - phid)* sigmad * y(6) - gammasd * y(12); % Symptomatic

% Infected Grads
dy(13) = phig * sigmag * y(7) - gammaag * y(13); % Asymptomatic
dy(14) = (1 - phig)* sigmag * y(7) - gammasg * y(14); % Symptomatic

% Infected Faculty 
dy(15) = phif * sigmaf * y(8) - gammaaf * y(15); % Asymptomatic
dy(16) = (1 - phif)* sigmaf * y(8) - gammasf * y(16); % Symptomatic

% Cumulative Infections
dy(17) = sigmau * y(5) + sigmad * y(6) + sigmag * y(7) + sigmaf * y(8);

% In health services
dy(18) = alphau * gammasu * y(10) - hu * y(18);
dy(19) = alphad * gammasd * y(12) - hd * y(19);
dy(20) = alphag * gammasg * y(14) - hg * y(20);
dy(21) = alphaf * gammasf * y(16) - hf * y(21);

% Cumulative at health services 
dy(22) = alphau * gammasu * y(10) + alphad * gammasd * y(12) + alphag * gammasg * y(14) + alphaf * gammasf * y(16); 

% Not quarantine / health services
dy(23) = (1 - alphau) * gammasu * y(10) - deltau * y(23);
dy(24) = (1 - alphad) * gammasd * y(12) - deltad * y(24);
dy(25) = (1 - alphag) * gammasg * y(14) - deltag * y(25);
dy(26) = (1 - alphaf) * gammasf * y(16) - deltaf * y(26);

% Recovered 
dy(27) = gammaau * y(9)  + hu * y(18) + deltau * y(23)-tR*y(27); % Ru
dy(28) = gammaad * y(11) + hd * y(19) + deltad * y(24)-tR*y(28); % Rd
dy(29) = gammaag * y(13) + hg * y(20) + deltag * y(25)-tR*y(29); % Rg
dy(30) = gammaaf * y(15) + hf * y(21) + deltaf * y(26)-tR*y(30); % Rf

% Vaccinated: U,D,G,F
% 31: Vu (new) % 32: Vd (new)
% 33: Vg (new) % 34: Vf (new)
dy(31) = -tV*y(31);
dy(32) = -tV*y(32);
dy(33) = -tV*y(33);
dy(34) = -tV*y(34);

dy = dy';
end