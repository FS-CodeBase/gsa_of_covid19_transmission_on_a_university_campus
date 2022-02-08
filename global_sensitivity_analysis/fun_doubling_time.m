function doubling_time_res = fun_doubling_time(cumulative_cases)
% FUNCTION: fun_measure
% AUTHOR: Fabian Santiago
% EMAIL: fsantiago3@ucmerced.edu
% DATE: 8/27/2021
% DESCRIPTION: Computes the average doubling time based on the cumulative 
%              cases in the first num_days.
%
num_days = 30;
doubling_time_res = mean(...
                    log(2)./log(cumulative_cases(2:num_days,:)./...
                                cumulative_cases(1:(num_days-1),:))...
                            ,1); 
end

