function init_list = initialize_random_compound (model, DT, K, fluo_values_cmp, n_traces)
    % Returns random initialization values for the HMM model
    % 
    % INPUTS
    % model: model type, see generate_input_matrix.m
    % K: number of naive states - ignore this if model ~= Nstate
    % fluo_values_cmp: fluorescence sequences for all traces in the same
    %                  AP position
    % n_traces: number of traces in the same ap position
    % 
    % OUTPUTS
    % Randomly initialized values for model parameters:
    % pi0: initial probability distribution of naive states
    % A: transition matrix
    % v: emission values where v[1]=0, corresponding to the OFF state
    % noise: S.D. of the gaussian noise

    
    
    % random A generation
    A=param_to_mat(model,DT);
    K=size(A,1);
    
    % random pi0 generation
    pi0 = rand(K, 1);
    pi0 = pi0/sum(pi0);
    
    % v and noise auxiliary
    fluo_values = [];
    for i_tr = 1:n_traces
        fluo_values = [fluo_values fluo_values_cmp{i_tr}(:)'];
    end
    
    % random v generation
    v = rand() * max(fluo_values(:)) * (0:(K-1))/(K-1);
    
    % random noise genertion
    noise = rand() * mean(fluo_values);
    
    % combine the initialized values into a list
    init_list = cell(4, 1);
    init_list{1} = pi0;
    init_list{2} = A;
    init_list{3} = v;
    init_list{4} = noise;