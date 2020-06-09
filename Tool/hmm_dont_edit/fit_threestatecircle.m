function [Afit,beta_]=fit_threestatecircle(A,DT)
% Find the constraint matrix Afit that maximize sum(log(Afit).*A)
% Afit is generated from a threestatecircle model
% INPUT:
%       A: input transition matrix, representing the posterior distribution
%       of (Zt,Z(t-1)).
%       DT: sampling time
% OUPUT:
%       Afit: Fitted matrix, generated with the parameter beta_
%       beta_: fitted parameters, that maximize the likelihood function sum(log(Afit).*A)
            % beta(1): sumk
            % beta(2): pon
            % beta(3): k1k2
    options = optimset('TolX',1e-6,'TolF',1e-6);
    fun = @(beta) param_to_mat('threestatecircle',DT,beta);
    obj = @(x) -sum(sum(log(fun(x)).*A));
    beta_=fminsearchbnd(obj,[1e-2 1e-2 1e-2],[1e-7 1e-7 1e-7],[1e1 1e1 1e1],options);
    Afit=fun(beta_);