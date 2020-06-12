addpath '..\..\Tool\hmm_dont_edit';
addpath '..\..\Tool\hmm_dont_edit\utilities\';
%%
Dataset='Hb-M1-New';
cycleno=13;
filename = ['trace_' Dataset '_nc' num2str(cycleno) ];
%% Load traces
load(['trace_data/' filename '.mat']);

bg=-1e5;            % Background intensity (w no bg intensity, set bg to -inf)
%% Inference setting
% setting
model_range={};
model_range{1}='twostate';
%model_range{2}='threestate';
n_steps_max = 50;   % maximum number of forward-backward iterations
eps = 10^(-4);      % tolerance variable for convergence
n_local = 100;       % number of local em runs
% Params to infer
tobeinfer=[1 1 1];                  % Which parameters need infering (A,v,noise)
%% Begin infering based on trace_rec
local_em_outputs_best={};

for model_inference = 1:numel(model_range)
    logL_max = -Inf;A_inf=NaN;R_inf=NaN;v_inf=NaN;
    updated=0;
    % Find the default emission vector of model_inference
    [A_]=param_to_mat(model_range{model_inference},dt);
    K=size(A_,1);
    tr=numel(trace_all);
    parfor i_local = 1:n_local
        % random initialization of model parameters
        init_list = initialize_random_compound (model_range{model_inference},dt, K, trace_all(1:tr), tr);

        pi0_init = init_list{1};
        K_init=numel(pi0_init);
        
        A_init = init_list{2};
        v_init = init_list{3};
        noise_init = init_list{4};

        % localEM call
        local_em_outputs = local_em_compound_bg_wd (model_range{model_inference},trace_all(1:tr), time_all(1:tr), bg, tr, ...
            A_init, dt, v_init, noise_init, tobeinfer, ...
            pi0_init', K_init, wL, n_steps_max, eps);
        display([num2str(i_local) ' - Done-' num2str(max(local_em_outputs.logL(end)))]);
        local_em_outputs.logL_max = max(local_em_outputs.logL);
        %local_em_outputs = rmfield(local_em_outputs,{'p_s_log_cmp','p_ss_log_cmp'});
        em_output_rec(i_local)=local_em_outputs;
    end
    % save better inference results
    for i_local=1:n_local
        local_em_outputs = em_output_rec(i_local);
        if max(local_em_outputs.logL_max) > logL_max
            %if ~any(~isreal(prob_to_rate(exp(local_em_outputs.A_log),dt)))
                logL_max = max(local_em_outputs.logL);
                
                local_em_outputs.A = exp(local_em_outputs.A_log);
                local_em_outputs.R = prob_to_rate(local_em_outputs.A,dt);
                local_em_outputs.v = exp(local_em_outputs.v_logs).*local_em_outputs.v_signs;
                local_em_outputs.lambda = exp(local_em_outputs.lambda_log);
                local_em_outputs.sigma = 1/sqrt(local_em_outputs.lambda);
                
                local_em_outputs_best{model_inference} = local_em_outputs;
                updated=updated+1;
            %end
        end
    end
    clear local_em_outputs;
end

%% Calculate noise level
% Get data relative noise level
T = min(cell2mat(time_all));
tr=numel(time_all);
stmp = zeros(1,tr);
for i=1:tr
    stmp(i)=mean(trace_all{i}(1:T));
end
noise_data = sqrt(var(stmp))/mean(stmp);
% Get model prediction noise
for model_inference = 1:numel(model_range)
    output=local_em_outputs_best{model_inference};
    % Find noise at SS
        % solve promoter state probability
        K = size(output.R,1);
        R_=[output.R;ones(1,K)];
        B=zeros(K+1,1);
        B(end)=1;
        mX = linsolve(R_,B);
        % Find variance from gaussian noise
        S2_gauss = output.sigma^2/T;
        % Find mean and noise
        [CV,S2,S1]=Pattern_CV_SS(output.R,mX,output.v,T*dt);
        noise_model = (S2+S2_gauss*(T*dt)^2-S1^2)/S1^2;
        % Find noise out of SS
        X0 = (output.A^w)*exp(output.pi0_log)';
        [CV,S2,S1]=Pattern_CV_oSS(output.R,X0,output.v,T*dt);
        noise_model = (S2+S2_gauss-S1^2)/S1^2;        
end
%% Save results
save(['data/result_' Dataset '_nc' num2str(cycleno) '_bg' num2str(bg)],'local_em_outputs_best','model_range','noise_data','noise_model','wL','dt');
%% Display results
model_inference=1;
display(['Noise for model ' num2str(model_inference)]);
display([noise_data noise_model])
display('Rate matrix: ');
display(local_em_outputs_best{1}.R)
kon=local_em_outputs_best{1}.R(2,1);
koff=local_em_outputs_best{1}.R(1,2);
display(['Emission vector: ' num2str(local_em_outputs_best{1}.v)]);
%local_em_outputs_best{1}.sigma
%pon_oSS=(local_em_outputs_best{1}.A^w)*exp(local_em_outputs_best{1}.pi0_log)';
%pon_oSS=pon_oSS(2)/sum(pon_oSS);
%display(['pon oo SS: ' num2str(pon_oSS)]);
pon_SS=kon/(kon+koff);
display(['pon from model at SS: ' num2str(pon_SS)]);
maxI = sum(local_em_outputs_best{1}.v)*sum(wL);
tmp=[];
nON=0;
for i=1:total
    tmp =[tmp  mean(trace_all{i})/maxI];
    nON = nON + any(trace_all{i}>0);
end
display(['pon from data at SS: ' num2str(mean(tmp))]);
pSpot = ((1/koff)+w*dt)/(1/koff+1/kon);
display(['pSpot at SS: ' num2str(pSpot)]);
mIntensity = sum([1-pon_SS pon_SS].*local_em_outputs_best{1}.v)*mean(wL);
display(['mean intensity at SS: ' num2str(mIntensity)]);
T = min(cell2mat(time_all));
pON = kon/(kon+koff)+koff/(kon+koff)*expcdf(T*dt,1/kon);
display(['pON (any spot or not over whole trace) from model: ' num2str(pON)]);
display(['pON (any spot or not over whole trace) from data: ' num2str(nON/total)]);
% Period
period = 1/kon + 1/koff;
display(['period: ' num2str(period)]);
%% Veterbi to find fitted signal
idx=60;
viterbi_out = viterbi_wL(trace_all{idx}(:)', local_em_outputs_best{1}.v(:), local_em_outputs_best{1}.sigma, local_em_outputs_best{1}.pi0_log(:), ...
                                     local_em_outputs_best{1}.A_log, numel(local_em_outputs_best{1}.v), wL);
plot(trace_all{idx});hold on;
plot(viterbi_out.fluo_viterbi(w+1:end));hold on;
stem((viterbi_out.z_viterbi(w+1:end)-1)*sum(wL));