function local_em_outputs = local_em_compound_bg_wd (model,fluo_values_cmp, time_cmp, beta, ...
    n_traces,  A, DT, v, noise, tobeinfer, pi0, K, wL, n_steps_max, eps)

    % The function returns the maximum likelihood estimates of the model
    % parameters, given multiple data sequences and initial conditions
    % 
    % INPUTS
    % model: model type
    % fluo_values_cmp: cell data type with multiple time series fluorescence
    %                  values that come from the same AP position
    % time_cmp: cell data type of the lengths of the fluorescence sequences
    % n_traces: number of different traces at the same AP position
    % A: initial transition matrix
    % DT: sampling time
    % v: initial emission values
    % noise: initial gaussian noise
    % pi0: initial naive state cmf at t=0
    % K: number of naive states
    % wL: memory configuration
    % n_steps_max: maximum number of backward-forward iterations
    % eps: relative change in model parameters for which the backward-forward
    % iterations stop and the algorithm is assumed to converge
    %
    % OUTPUTS
    % Inferred model parameters when:
    % - the relative change in model parameters is smaller than eps,
    % - or the number of backward-forward iterations exceeds n_steps_max
    % - account for censored signal due to background noise
    % - inference over time window (not requiring F(t=0)=0)
    % 
    % pi0_log: log of the inferred initial state cmf
    % A_log: log of the inferred transition matrix
    % v_logs: logs of the absolute values of the emission parameters
    % v_sings: sings of the inferred emission parameters
    % lambda_log: log of the inferred lambda parameter defined as 1/noise^2
    % baum_welsh: the number of forward-backward iterations done before
    %             convergence
    % log_likelihoods(1:baum_welsh): log-likelihood of observing the
    %                                fluoresence data evaluated using the
    %                                inferred model parameters
    % p_s_log_cmp: expected log-probabilities of each compound state at all
    %              time points for each trace
    % p_ss_log_cmp: expected log-probabilities of joint compound states at all
    %               time points for each trace
    
    w=numel(wL);
    % ------------------------ Log converstions ------------------------
    pi0_log = log(pi0);
    A_log = log(A);
    v_logs = log(abs(v));
    v_signs = sign(v);
    lambda_log = -2*log(noise);
    
    fluo_logs_cmp = cell([n_traces, 1]);
    fluo_sign_cmp = cell([n_traces, 1]);
    fluo_censored_cmp = cell([n_traces, 1]);
    
    censored_val=-100;
    
    censored_ever=false;

    % Check if signal is censored and adding unknown value at beginning
    for i_tr = 1:n_traces
        % patching signal by unknown initial values
        fluo_values_cmp{i_tr}=[zeros(1,w) fluo_values_cmp{i_tr}];
        time_cmp{i_tr}=time_cmp{i_tr}+w;        
        % Values
        fluo_censored_cmp{i_tr} = fluo_values_cmp{i_tr}==beta;
        if sum(fluo_censored_cmp{i_tr})
            censored_ever=true;
        end
        fluo_values_cmp{i_tr}(fluo_censored_cmp{i_tr})=censored_val;
        fluo_logs_cmp{i_tr} = log(abs(fluo_values_cmp{i_tr}));
        fluo_sign_cmp{i_tr} = sign(fluo_values_cmp{i_tr});
    end

    
    % ---------------------- Accounting variables ----------------------
    
    % C_{zm} = 1 iff the first digit of compound state m is equal to z
    C = zeros(K, K^w);
    for i = 1:K
        columns = ((i-1)*K^(w-1) : ((i*K^(w-1))-1)) + 1;
        C(i, columns(:)) = 1;
    end

    % log of the accounting variable C
    logC = log(C);

    % F_{zm} is the number of times the naive state z is encountered in the
    % naive state representation of the compound state m
    F = zeros(K, K^w);
    for i = 1:K^w
        F(:,i) = naive_count(i, K, wL);        %%%%% Fix here to account for non uniform wL
    end

    
    % ----------------------------- Lists ------------------------------
    
    % list of states that can transition to each of 1:K^w compounds states
    allowed_to_list = zeros(K^w, K, 'int16');
    
    % list of states allowed to go from each of 1:K^w compound states
    allowed_from_list = zeros(K^w, K, 'int16');
    
    % list of the first digits (naive states) of each compound state
    digit_first_list = zeros(1, K^w, 'int16');
    
    % list of the number of each naive state encountered in compound states
    naive_count_list = zeros(K^w, K, 'int16');
    
    for i = 1:K^w
        allowed_to_list(i,:) = allowed_to(i, K, w);
        allowed_from_list(i,:) = allowed_from(i, K, w);
        digit_first_list(i) = digit(i, 1, K, w);
        naive_count_list(i, :) = naive_count (i, K, wL);
    end          
    
    % log of the naive count list
    naive_count_log_list = log(double(naive_count_list));
    
    % list of compound states that have a unique sequence of the number of
    % naive states contained in the base-K representation
    unique_naive_list = unique(naive_count_list,'rows');
    [n_unique, ~] = size(unique_naive_list);    
    
    % a list that maps the unique naive count combinations with compound
    % states
    naive_count_map = cell(n_unique, 1);   
    for i = 1:n_unique
        naive = unique_naive_list(i,:);
        naive_count_map{i} = find(ismember(naive_count_list,naive,'rows'));        
    end
    
    % list of possible compound state at all times
    time_max = max(cell2mat(time_cmp));
    possible_states_list = cell(time_max, 1);    
    for t = 1:time_max
        possible_states_list{t} = possible_states (t, K, w);        
    end
    
    % calculation of p_ss_log indices used in A_log maximization
    A_log_maximization_ind_cmp = cell([n_traces, 1]);
    for i_tr = 1:n_traces
        A_log_maximization_ind_cmp{i_tr} = cell(K);
        
        ind_addition_mat = repmat((1:(time_cmp{i_tr}-1))*K^(2*w), [K^(w-1), 1]);
        ind_addition_list = ind_addition_mat(:);
        
        for m = 1:K
            for n = 1:K
                d_n_list = find(digit_first_list == n);
                d_m_list = find(digit_first_list == m);
                allowed_from_n = allowed_from_list(d_n_list, :);
                column_ind = ismember(allowed_from_n(1,:), d_m_list) == 1;
                ind_2d = sub2ind([K^w, K^w], allowed_from_n(:, column_ind), d_n_list');
                ind_3d = repmat(ind_2d, [1, time_cmp{i_tr}-1]);
                A_log_maximization_ind_cmp{i_tr}{m,n} = ind_addition_list + ind_3d(:);
            end
        end
    end
    
    % calculation of the log of the F variable used in v maximization
    log_F_terms_cmp = cell([n_traces, 1]);
    for i_tr = 1:n_traces
        log_F_terms_cmp{i_tr} = cell([K, 1]);
        for n = 2:K
            log_F_terms_cmp{i_tr}{n} = repmat(log(F(n,:))', 1, time_cmp{i_tr});
        end
    end
    
    % log's and signs of the fluorescence data used in v maximization
    x_term_log_cmp = cell([n_traces, 1]);
    x_term_sign_cmp = cell([n_traces, 1]);
    for i_tr = 1:n_traces
        x_term_log_cmp{i_tr} = repmat(log(abs(fluo_values_cmp{i_tr})), K^w, 1);
        x_term_sign_cmp{i_tr} = repmat(sign(fluo_values_cmp{i_tr}), K^w, 1);
    end
    
    
    % ------------------------ Pre-allocations -------------------------
    % log likelihoods of observing the fluorescence sequence at each EM
    % forward-backward iteration    
    log_likelihoods_cmp = cell([n_traces, 1]);
    for i_tr = 1:n_traces
        log_likelihoods_cmp{i_tr} = -Inf(1, n_steps_max);
    end
    log_likelihoods = -Inf(1, n_steps_max);
    
    % pre-allocation of alpha and beta coefficient matrices
    alpha_matrix_log_cmp = cell([n_traces, 1]);
    beta_matrix_log_cmp = cell([n_traces, 1]);
    for i_tr = 1:n_traces
        alpha_matrix_log_cmp{i_tr} = zeros(K^w,time_cmp{i_tr});                
        beta_matrix_log_cmp{i_tr} = zeros(K^w,time_cmp{i_tr});
    end

    % pre-allocation of the p_ss_log array used in the expectation step
    p_ss_log_cmp = cell([n_traces, 1]);
    for i_tr = 1:n_traces
        p_ss_log_cmp{i_tr} = -Inf(K^w, K^w, time_cmp{i_tr});
    end
    
    
    % ------------- Expectation-Maximizaiton iterations ----------------
    for baum_welsh = 1:n_steps_max
        
        %%%%%%%%%%%%%%%%%%%%%%%%% EXPECTATION %%%%%%%%%%%%%%%%%%%%%%%%%%
    
        % ------------- lists used in the expectation step -------------
        [eta_log_list_cmp,difference_list_cmp]=log_likelihood_censored_wd(n_traces,n_unique,K,w,naive_count_map,fluo_logs_cmp,fluo_sign_cmp,time_cmp,v_logs,v_signs,naive_count_log_list,lambda_log,beta,censored_val,fluo_censored_cmp);
        
        % row and column subscripts used for indexing A_log elements for
        % alpha matrix calculation
        A_log_alpha_rowSubs = digit_first_list(repmat(1:K^w, [1, K]));
        A_log_alpha_colSubs = digit_first_list(allowed_to_list(:));
        
        % A_log element indexing in a 2d slice used for alpha matrix
        % calculation
        A_log_alpha_subs_1d = A_log(sub2ind([K, K], A_log_alpha_rowSubs, A_log_alpha_colSubs));
        alpha_A_log_list = reshape(A_log_alpha_subs_1d, [K^w K]);        
        
        % row and column subscripts used for indexing A_log elements for
        % alpha matrix calculation
        A_log_beta_rowSubs = digit_first_list(allowed_from_list(:));
        A_log_beta_colSubs = digit_first_list(repmat(1:K^w, [1, K]));
        
        % A_log element indexing in a 2d slice used for beta matrix
        % calculation
        A_log_beta_subs_1d = A_log(sub2ind([K, K], A_log_beta_rowSubs, A_log_beta_colSubs));
        beta_A_log_list = reshape(A_log_beta_subs_1d, [K^w K]);                       
        
        
        % ------------ alpha coefficient matrix calculation ------------ 
                
        for i_tr = 1:n_traces
            
            % initializes the alpha matrix as all zeros (logs all infinities)            
            alpha_matrix_log_cmp{i_tr}(:,:) = -Inf;
            
            % calculates the alpha matrix elements for t = 1
            for i = possible_states_list{1}
                alpha_matrix_log_cmp{i_tr}(i, 1) = eta_log_list_cmp{i_tr}(i, 1) + ...
                    pi0_log(digit_first_list(i));
            end                   
        
            % calculates the alpha matrix elements for t > 1
            for t = 2:time_cmp{i_tr}
                % possible states at time t
                i_possible = possible_states_list{t};

                % list of terms that are added to find the alpha matrix
                % elements
                alpha_terms_list = alpha_A_log_list(i_possible, :) + ...
                    reshape(alpha_matrix_log_cmp{i_tr}(allowed_to_list(i_possible,:),t-1), [], K);             

                % local execution of the vectorized log_sum_exp_positive
                alpha_terms_max = max(alpha_terms_list, [], 2); 
                alpha_terms_diff = alpha_terms_list - repmat(alpha_terms_max, [1, K]);
                alpha_terms_diff_sum_exp_log = log(sum(exp(alpha_terms_diff), 2));           

                % assignment of alpha matrix element values
                alpha_matrix_log_cmp{i_tr}(i_possible,t) = alpha_terms_diff_sum_exp_log + ...
                    alpha_terms_max + eta_log_list_cmp{i_tr}(i_possible, t);                             
            end
        end

        
        % ------------ beta coefficient matrix calculation -------------
        
        for i_tr = 1:n_traces
            
            % initializes the alpha matrix as all zeros (logs all infinities)
            beta_matrix_log_cmp{i_tr}(:,:) = -Inf;       

            % assigns 1 to beta matrix elements at t=time (or 0's to log's)
            beta_matrix_log_cmp{i_tr}(:, time_cmp{i_tr}) = 0;

            % calculates the alpha matrix elements for t < time
            for t = (time_cmp{i_tr}-1: -1: 1)
                % possible states at time t
                i_possible = possible_states_list{t};

                % list of terms that are added to find the beta matrix elements
                beta_terms_list = beta_A_log_list(i_possible, :) + ...
                    reshape(eta_log_list_cmp{i_tr}(allowed_from_list(i_possible,:), t+1), [], K) + ...
                    reshape(beta_matrix_log_cmp{i_tr}(allowed_from_list(i_possible,:), t+1), [], K);

                % local execution of the vectorized log_sum_exp_positive
                beta_terms_max = max(beta_terms_list, [], 2);
                beta_terms_diff = beta_terms_list - repmat(beta_terms_max, [1, K]);
                beta_terms_diff_sum_exp_log = log(sum(exp(beta_terms_diff), 2));

                % assignment of beta matrix element values
                beta_matrix_log_cmp{i_tr}(i_possible,t) = beta_terms_max + ...
                    beta_terms_diff_sum_exp_log; 
            end
        end
        
        
        % --------------- log-likelihood calcuation --------------------
        log_likelihoods(baum_welsh) = 0;
        for i_tr = 1:n_traces
            log_likelihoods_cmp{i_tr}(baum_welsh) = ...
                log_sum_exp_positive(alpha_matrix_log_cmp{i_tr}(:,time_cmp{i_tr}));
            log_likelihoods(baum_welsh) = log_likelihoods(baum_welsh) + ...
                log_likelihoods_cmp{i_tr}(baum_welsh);
        end
        
        
        % --------------------- <S_t> calculation ----------------------
        p_s_log_cmp = cell([n_traces, 1]);
        for i_tr = 1:n_traces
            p_s_log_cmp{i_tr} = alpha_matrix_log_cmp{i_tr} + beta_matrix_log_cmp{i_tr} ...
                - log_likelihoods_cmp{i_tr}(baum_welsh);
        end
        

        % ------------------- <S_t, S_{t-1}> calculation ---------------

        for i_tr = 1:n_traces
            % list of additional 1d indices to account for multiple time
            % layers in the 3d representation of the p_ss_log matrix
            ind_addition_mat = repmat((1:(time_cmp{i_tr}-1))*K^(2*w), [K^(w+1), 1]);
            ind_addition_list = ind_addition_mat(:);

            % list of 1d indices of the positions of pairs in the 2d slice of
            % p_ss_log which correspond to valid col -> row transitions
            ind_positions_2d = sub2ind([K^w, K^w], allowed_from_list(:), repmat(1:K^w, [1, K])');

            % replication of the 1d pair indices along the time axis
            ind_positions_3d = repmat(ind_positions_2d, [1, (time_cmp{i_tr}-1)]);

            % combination of the pair indices with the additional 1d indices
            ind_3d = ind_addition_list + ind_positions_3d(:);

            % replication of the alpha matrix K times along one of the 2d axes
            alpha_matrix_log_minus = alpha_matrix_log_cmp{i_tr}(:,1:(time_cmp{i_tr}-1));
            alpha_matrix_log_minus_rep = repmat(alpha_matrix_log_minus, [K, 1]);                        

            % 2d indexing of eta_log and beta_log matrices
            % note: the 2d slice components are accounted for through the
            %       1d representation of the matrix of possible transitions.
            %       Also, note that the 1d indexing of the p_ss_log matrix
            %       is done by rows, which is used in calculating the terms
            eta_log_list_minus = eta_log_list_cmp{i_tr}(allowed_from_list(:), 2:time_cmp{i_tr});
            beta_matrix_log_minus = beta_matrix_log_cmp{i_tr}(allowed_from_list(:), 2:time_cmp{i_tr});        

            % row and column subscripts used for indexing A_log elements
            A_log_rowSubs = digit_first_list(allowed_from_list(:));
            A_log_colSubs = digit_first_list(repmat(1:K^w, [1, K]));

            % A_log element indexing in a 2d slice
            A_log_sub_single = A_log(sub2ind([K, K], A_log_rowSubs, A_log_colSubs));

            % A_log element indexing in a 3d slice along the time axis
            A_log_colSubs_multi = repmat(A_log_sub_single, [1, time_cmp{i_tr}-1]);                

            % calcuation of p_ss_log elements
            p_ss_log_cmp{i_tr}(ind_3d) = alpha_matrix_log_minus_rep(:) + ...
                eta_log_list_minus(:) + beta_matrix_log_minus(:) + ...
                A_log_colSubs_multi(:) - log_likelihoods_cmp{i_tr}(baum_welsh);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%% MAXIMIZATION %%%%%%%%%%%%%%%%%%%%%%%%%

        % -------------------------- pi0_log ---------------------------
        for m = 1:K
            pi0_terms = [];
            for i_tr = 1:n_traces
                pi0_terms = [pi0_terms, logC(m, :) + p_s_log_cmp{i_tr}(:,1)'];
            end
            pi0_log(m) = log_sum_exp_positive(pi0_terms) - log(n_traces);
        end


        % --------------------------- A_log ----------------------------
        if tobeinfer(1)
            % Collapse the compound transition matrix to weight matrix B (not normalized)
            A_old = exp(A_log);
            for m = 1:K
                for n = 1:K
                    arr = [];
                    for i_tr = 1:n_traces
                        p_ss_log_ith = p_ss_log_cmp{i_tr}(A_log_maximization_ind_cmp{i_tr}{m,n});
                        arr = [arr, p_ss_log_ith(:)'];
                    end
                    arr_max = max(arr(:));
                    if arr_max==-Inf
                        B_log_(m,n)=-Inf;
                    else
                        B_log_(m,n) = arr_max + log(sum(exp(arr(:)-arr_max)));                        
                    end
                    B_log(m,n) = log(sum(exp(arr(:))));
                end
            end
            for n = 1:K
                arr = B_log_(:,n);
                arr_max = max(arr(:));
                B_log_(:,n) = B_log_(:,n) - (arr_max + log(sum(exp(arr(:)-arr_max))));
            end
            % Normalized by max value
            arr_max=max(B_log(:));
            B_log = B_log - arr_max;
            
            B=exp(B_log);
            B_ = exp(B_log_);
            % Adjusting B:
            
            % Adding model constraint on A
            switch model
                case 'twostate'
                    A_log = B_log_;
                case 'threestate'
                    A_log = B_log_;
                case 'threestatecircle'
                    A_log = fit_threestatecircle(B,DT);
                case 'twocopycorr'
                    A_log = B_log_;
                case 'twocopyuncorr'
                    A_log = fit_twocopyuncorr(B,DT);
                case 'twocopyalpha'
                    A_log = fit_twocopyalpha(B,DT);
                otherwise
                    A_log = B_log_;                        
            end
            %A_log=log(A);
            A=exp(A_log);
            A_norm_change = abs(norm(A) - norm(A_old))/K;            
        else
            A_norm_change = 0;
        end


        % ----------------------------- v ------------------------------
        if tobeinfer(2)
            if ~censored_ever
                v_logs_old = v_logs;

                m_sign = ones(K-1, K-1);
                m_log = zeros(K-1, K-1);

                for m = 2:K
                    for n = 2:K
                        terms = [];
                        for i_tr = 1:n_traces
                            terms_ith = p_s_log_cmp{i_tr} + log_F_terms_cmp{i_tr}{n} + log_F_terms_cmp{i_tr}{m};
                            terms = [terms, terms_ith(:)'];
                        end
                        terms_max = max(terms(:));
                        m_log(m-1, n-1) = terms_max + log(sum(exp(terms(:)-terms_max)));
                    end
                end

                b_sign = ones(1, K-1);
                b_log = zeros(1, K-1);

                for m = 2:K
                    terms_b_log = [];
                    terms_b_sign = [];
                    for i_tr = 1:n_traces
                        terms_b_log_ith = x_term_log_cmp{i_tr} + p_s_log_cmp{i_tr} + ...
                            + log_F_terms_cmp{i_tr}{m};
                        terms_b_log = [terms_b_log, terms_b_log_ith(:)'];
                        terms_b_sign = [terms_b_sign, x_term_sign_cmp{i_tr}(:)'];
                    end
                    terms_b_sum = log_sum_exp(terms_b_log, terms_b_sign);
                    b_log(m-1) = terms_b_sum(1);
                    b_sign(m-1) = terms_b_sum(2);
                end
            
                v_updated = v_log_solve(m_log, m_sign, b_log, b_sign);
                v_logs = v_updated(1,:);
                v_signs = v_updated(2,:);
            else
                v_logs_old = v_logs;
                
                % Optimization here:
                v=fminsearchbnd(@obj_v,exp(v_logs).*v_signs,zeros(size(v_logs)),1e5*ones(size(v_logs)));
                v_logs = log(v);
                v_signs = sign(v);
            end

            v_norm_change = abs(norm(exp(v_logs_old)) - norm(exp(v_logs)))/K/norm(exp(v_logs_old));
        else
            v_logs=log(v+1e-10);
            v_signs = ones(size(v));
            v_norm_change = 0;
        end

        
                % ------------------------- lambda_log -------------------------
        if tobeinfer(3)
            if ~censored_ever
                noise_old = exp(-lambda_log/2);

                arr = [];
                for i_tr = 1:n_traces
                    term_ith = p_s_log_cmp{i_tr} + difference_list_cmp{i_tr};
                    arr = [arr, term_ith(:)'];
                end
                arr_max = max(arr(:));
                lambda_log = log(sum(cell2mat(time_cmp))) - (arr_max + log(sum(exp(arr(:)-arr_max))));

                noise = exp(-lambda_log/2);
                noise_change = abs(noise - noise_old);
            else
                noise_old = exp(-lambda_log/2);
                % Optimization here
                lambda_log=fminsearchbnd(@obj_lambda,lambda_log,-1e5,10);
                noise = exp(-lambda_log/2);
                noise_change = abs(noise - noise_old)/noise_old;
            end
        else
            noise_change = 0;
        end
        
        %[num2str(baum_welsh) '. '   num2str(noise_change) ' ' num2str(A_norm_change) '  ' num2str(log_likelihoods(baum_welsh))]
        
        % ------------------- convergence criterion --------------------
        if (max([noise_change, v_norm_change, A_norm_change]) < eps)
            break
        end
    end

    % -------------- collection of outputs into a cell -----------------

    local_em_outputs = struct('pi0_log', pi0_log, 'A_log', A_log, ...
        'v_logs', v_logs, 'v_signs', v_signs, 'lambda_log', lambda_log, ...
        'logL', log_likelihoods(1:baum_welsh), 'max_bw', baum_welsh);
    % Only keep the last P_s_log_cmp and P_ss_log_cmp
    local_em_outputs.p_s_log_cmp = p_s_log_cmp{end};
    local_em_outputs.p_ss_log_cmp = p_ss_log_cmp{end};

    % ------------- auxiliary function for inference of censored signal
    function obj=obj_v(v)
        eta_log_list_cmp_tmp=log_likelihood_censored_wd(n_traces,n_unique,K,w,naive_count_map,fluo_logs_cmp,fluo_sign_cmp,time_cmp,log(v),sign(v),naive_count_log_list,lambda_log,beta,censored_val,fluo_censored_cmp);
        obj=0;
        for jtmp=1:n_traces
            term_ith_tmp = exp(p_s_log_cmp{jtmp}).*eta_log_list_cmp_tmp{jtmp};
            obj=obj-sum(term_ith_tmp(:));
        end
    end

    function obj=obj_lambda(lambda_log)
        obj=0;
        eta_log_list_cmp_tmp = cell([n_traces, 1]);
        for jtmp=1:n_traces
            % Calculate the log-likelihood value:
            % For observed value
            tmp1 = 0.5*(lambda_log - log(2*pi)) ...
                -0.5*exp(lambda_log + difference_list_cmp{jtmp});
            % For censored value
            tmp2 = -(exp(difference_list_cmp{jtmp}/2)-beta-censored_val)*exp(lambda_log/2);
            tmp3=normcdf(tmp2);
            % Merging the censored and observed likelihood
            tmp1(:,fluo_censored_cmp{jtmp})=0;
            tmp3(:,~fluo_censored_cmp{jtmp})=0;
            eta_log_list_cmp_tmp{jtmp} = tmp1+tmp3;
            % Calculate the log-likihood
            term_ith_tmp = exp(p_s_log_cmp{jtmp}).*eta_log_list_cmp_tmp{jtmp};
            obj=obj-sum(term_ith_tmp(:));
        end
    end
end