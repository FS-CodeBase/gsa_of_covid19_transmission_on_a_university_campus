function PRMS = fun_lhs_sampling_of_parameters(prms_info,N)
% % FUNCTION: FUN_LHS_SAMPLING_OF_PARAMETERS
% % AUTHOR: Fabian Santiago
% % EMAIL: fsantiago3@ucmerced.edu
% % DATE: 10/27/2022
% % 
% % DESCRIPTION: Selected model parameters are sampled N times using Latin
% %              hypercube sampling. Each parameter is sampled from a  
% %              uniform distribution within the range specified in the 
% %              input prms_info.
% % 
% % ABBREVIATIONS: 
% %         U: undergraduates living off campus
% %         D: undergraduates living in the doorms
% %         G: graduate students
% %         F: faculty and staff
% %

% Error detection: Determine number parameters not being varied but still 
% contain a range for the parameter.
lb_num = sum(prms_info(~prms_info(:,1),2) ~= prms_info(~prms_info(:,1),3));
if sum(lb_num)~=0
    warning(['flag1 is set to false for ',num2str(sum(lb_num)),...
              ' parameters, and the lower_bound value is ',... 
             'used as the model parameter value for each, even though',...
             'an upper bound has been provided. ',... 
             'See fun_model_parameter_ranges.'])
end
% Determine the number of parameters and the number of parameters that 
% are being varied for sensitivity. 
numPRMS = size(prms_info,1);
gsaPRMS = sum(prms_info(:,1));

% Latin hypercube sampling (LHS) of parameters
LHS_PRMS = lhsdesign(N,gsaPRMS);

% Pre-allocate space for all of the paramters that will be sampled. 
PRMS = zeros(N,numPRMS);

% Sample parameter values from a uniform distribution in the specified
% range.
prm_lhs = 1;
for param = 1:numPRMS   
    % Determine range (a,b)
    [a,b] = deal(prms_info(param,2),prms_info(param,3));
    
    % Determine if parameter is in days or not, and if we are assesing the
    % sensitivity of that parameter.
    if logical(prms_info(param,4))
        if logical(prms_info(param,1))
            PRMS(:,param) = 1./(24*((b-a)*LHS_PRMS(:,prm_lhs)+a));
            prm_lhs = prm_lhs + 1;
        else
            PRMS(:,param) = 1./(24*(a*ones(N,1)));
        end
    else
        if logical(prms_info(param,1))
            PRMS(:,param) = (b-a)*LHS_PRMS(:,prm_lhs)+a;
            prm_lhs = prm_lhs + 1;
        else
            PRMS(:,param) = a*ones(N,1);
        end
    end
end

% Make pA is the same pM (M = pA + pM) 
PRMS(:,5) = PRMS(:,4);

% Make SIGMA the same for all U, D, G, and F 
PRMS(:,13:15) = repmat(PRMS(:,12),[1 3]);

% Make PHI the same for all U, D, G, and F 
PRMS(:,17:19) = repmat(PRMS(:,16),[1 3]);

% Make ALPHA the same for all U, D, G, and F
PRMS(:,29:31) = repmat(PRMS(:,28),[1 3]);