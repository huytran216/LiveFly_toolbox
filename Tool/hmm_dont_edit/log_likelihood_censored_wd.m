function [eta_log_list_cmp,difference_list_cmp]=log_likelihood_censored_wd(n_traces,n_unique,K,w,naive_count_map,fluo_logs_cmp,fluo_sign_cmp,time_cmp,v_logs,v_signs,naive_count_log_list,lambda_log,beta,censored_val,fluo_censored_cmp)
    % list of log (X_t - V_t)^2 terms that appear in the emission pdf's
            difference_list_cmp = cell([n_traces, 1]);
            eta_log_list_cmp = cell([n_traces, 1]);
            for i_tr = 1:n_traces
                difference_list = zeros(K^w, time_cmp{i_tr});
                for i = 1:n_unique
                    states = naive_count_map{i};
                    difference_list(states, :) = repmat(difference_new(fluo_logs_cmp{i_tr}, fluo_sign_cmp{i_tr}, time_cmp{i_tr}, K, ...
                       states(1), v_logs, v_signs, naive_count_log_list)', [length(states), 1]);
                end
                difference_list_cmp{i_tr} = difference_list;
                difference_list_cmp{i_tr}(:,1:w) = 0;
                % Calculate the log-likelihood value:
                    % For observed value
                    tmp1 = 0.5*(lambda_log - log(2*pi)) ...
                        -0.5*exp(lambda_log + difference_list_cmp{i_tr});
                    % For censored value
                    tmp2 = -(exp(difference_list_cmp{i_tr}/2)-beta-censored_val)*exp(lambda_log/2);
                    %tmp3 = log(0.5) + log(1+sign(tmp2).*sqrt(1-exp(-2*(tmp2.^2)/pi)));
                    tmp3=normcdf(tmp2);
                    % log-likelihood
                    tmp1(:,fluo_censored_cmp{i_tr})=0;
                    tmp1(:,1:w)=0;
                    tmp3(:,~fluo_censored_cmp{i_tr})=0;
                    eta_log_list_cmp{i_tr} = tmp1+tmp3;
                    if any(isnan(eta_log_list_cmp{i_tr}))
                        %'ops'
                    end
            end