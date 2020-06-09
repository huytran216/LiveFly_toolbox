function diff = difference_wL(fluo_logs, fluo_signs, time, K, wL, ...
                            state, v_logs, v_signs, naive_count_logs)
    % Calculates the log of (X_t - V(S_t))^2 for all values of t.
    % 
    % INPUTS
    % fluo_logs: the logs of the fluorescence values at all times
    % fluo_sign: the signs of the fluorescence values at all times
    % time: the length of the fluorescence sequence
    % K: number of naive states
    % wL: MS2 configuration vector
    % state: compound state of the system    
    % v_logs: log of the absolute value of the emission parameters
    % v_signs: signs of the emission parameters
    % naive_count_logs: the list of the logs of naive counts for all states    
    % 
    % OUTPUTS
    % diff: logs of (X_t - V(S_t))^2 for all values of t, where V(S) is 
    %       the compound flurescence at state S. Note: underflow issues are 
    %       avoided in the calculation.
    
    v_multi_log = v_logs + naive_count_logs(state,:);
    v_multi_log_time = repmat(v_multi_log, [time, 1]);
    
    naive_counts = exp(naive_count_logs(state,:));
    state1_count = naive_counts(1);   
    w=numel(wL);
    ms2_coeff = wL;
    
    if state1_count > 0
        for t = 1:(w-1)
            count_reduction = sum(ms2_coeff((t+1):end));
            v_multi_log_time(t, 1) = v_multi_log_time(t, 1) - log(state1_count) + ...
                log(abs(state1_count - count_reduction));
        end
    end
    
    v_multi_sign_time = repmat(v_signs, [time, 1]);
    logs_combined = [fluo_logs', v_multi_log_time];
    signs_combined = [fluo_signs', -v_multi_sign_time];
    
    terms_max = max(logs_combined, [], 2); 
    terms_diff = logs_combined - repmat(terms_max, [1, (K+1)]);    
    terms_diff_with_signs = signs_combined.*exp(terms_diff);
    terms_diff_sum_with_signs = sum(terms_diff_with_signs, 2);
    diff = 2*(terms_max + log(abs(terms_diff_sum_with_signs)));
    
    % account for cases when all terms are -Inf
    diff(isnan(diff)) = -Inf;
    %error('Cleared')