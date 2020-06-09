function [answer] = log_sum_exp_positive (arr)
    % Calculates the sum of logs.
    % 
    % INPUTS
    % arr: list of log values
    % 
    % OUTPUTS
    % The sum of log values, avoiding the underflow issues.
    
    arr_max = max(arr(:)); 
    answer = arr_max + log(sum(exp(arr(:)-arr_max)));
    if (isnan(answer))
        answer = -Inf;        
    end