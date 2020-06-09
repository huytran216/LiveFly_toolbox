function [R] = prob_to_rate (A, deltaT)
    % Calculates the transition rate matrix using the transition probabilities
    % 
    % INPUTS
    % A: transition probability matrix
    % deltaT: time interval between consecutive transitions
    % 
    % OUTPUTS
    % The transition rate matrix
    
    R = logm(A)/deltaT;