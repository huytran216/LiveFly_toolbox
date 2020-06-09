%This code investigates checks accuracy of EM inferences, using Hessians to
%estimate lower bound on variance of estimates.
%       1) Do true values fall withing precision bounds of inferred values?
%       
%       2) How sensitive is the accuracy to data characteristics? Is
%       it consistently a solid estimate? Or can we only employ for certain
%       trace counts and/or lengths? This will be explored analytically as
%       well.
%       3) How does accuracy vary for different parameter sets?
%% I. Check accuracy: 1 inferences on 10 loops length:100
V_values = [0, 3, 6];
Noise = 1/sqrt(exp(1));
Pi0 = [0.4, 0.3, 0.3];
KK = 3;
ww = 5;
AA = [0.3, 0.2, 0.6; 0.2, 0.55, 0.3 ; 0.5, 0.25, 0.1];
%   Define max number of steps and likelihood threshold
%   for each EM inference run
N_steps_max = 1000;
Eps = 10^(-4);
%   Set number of independent inference loops to run
INF_L = [1];
trace_length = [100];
N_traces = [50];
%   comp_inf_precision(seq_length, inf_loops, n_traces, n_steps_max, eps, ...
%   v_synth, noise_synth, pi0_synth, K, w, A)
tic
hessian_spread = comp_inf_precision(trace_length, INF_L,...
    N_traces, N_steps_max, Eps, V_values, Noise, Pi0, KK, ww, AA);
toc
%% I continued: Assess results 

%Transition Parameters
h_std_A = hessian_spread{3};
c_mean_A = reshape(hessian_spread{1},1,[]);
true_mean_A = reshape(AA,1,[]);

errorbar(1:KK^2, c_mean_A,h_std_A,'b O');
hold on
plot(1:KK^2,true_mean_A, 'black X');

grid on
title('Inference Accuracy: Inferred Parameter Values vs. Actual (Transition Probabilities)');
ylabel('Transition Parameter (linear indexing)');
xlabel('Transition Probability');

% legend('10 EM Runs','100 EM Runs','1000 EM Runs','Location','southeast')
