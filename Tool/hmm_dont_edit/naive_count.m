function [counts] = naive_count (state, K, wL)
    % Returns the number of each naive state at the given compound state.
    % 
    % INPUTS
    %   state: compound state
    %   K: number of naive states
    %   wL: ms2 configurations
    % 
    % OUTPUTS
    %   counts: The number of times each naive state is accessed in the given compound state.
    
    naive = compound_to_naive(state, K, numel(wL));
    counts = zeros(1,K);
    for i = 1:K
        counts(i) = sum((naive==i).*wL(end:-1:1));    % Fix to account for nonuniform wL
    end