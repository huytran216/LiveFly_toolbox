function [A] = rate_to_prob (R, deltaT)
    % Calculates the transition probability matrix using the transition rates
    % 
    % INPUTS
    % R: transition rate matrix
    % deltaT: time interval between consecutive transitions
    % 
    % OUTPUTS
    % Transition probability matrix.
    
    A = expm(R*deltaT);