function [y] = hybrid_power_exponential(x0,lambda, beta0, ax)
    % equal 1: at x0
        % expo function: y = y1_coeff*exp((-x)/lambda)
        % expo function: y = y2_coeff*((x+beta1)^(-beta2))
    y1_coeff = exp(x0/lambda);
    y2_coeff = (x0+beta0(1))^beta0(2);
    
    
    ax1 = ax<x0;
    y = ax;
    y(ax1) = y1_coeff*exp(-ax(ax1)/lambda);
    y(~ax1) = y2_coeff*(ax(~ax1)+beta0(1)).^(-beta0(2));
    y=y/max(y)