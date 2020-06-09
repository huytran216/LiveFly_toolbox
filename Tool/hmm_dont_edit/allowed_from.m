function [l] = allowed_from(state, K, w) 
    % List of all compound states that the specified compound state can transition to.
    % 
    % INPUTS
    % state: compound state
    % K: number of naive states
    % w: memory
    % 
    % OUTPUTS
    % List of all compound states that the specified compound state can transition to.
    % The fact that should be an overlap of (w-1) naive states in the present and future
    % is used in the calculations.
    
    last_naive = mod(state-1, K);
    front = floor( (state-1-last_naive) / K);
    l = front + K^(w-1) * (0:(K-1)) + 1;