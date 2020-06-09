function [answer] = log_sum_exp (arr, signs)
    % Calculates the sum of logs with specified signs.
    % 
    % INPUTS
    % arr: list of log values
    % sign: list of signs corresponding to each log values
    % 
    % OUTPUTS
    % The sum of log values with specified signs to allow subtraction.
    % The calculation resolves the underflow issue.
    
    arr_max = max(arr(:));
    term2_array = signs.*exp(arr-arr_max);
    term2 = sum(term2_array(:));
    answer = [arr_max + log(abs(term2)), sign(term2)];
    if (isnan(answer(1)))
        answer = [-Inf, 0];        
    end