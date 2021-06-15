idx_list = [1 2 3 4 10];

isBcd1X = [0 1];

fea_range = [1];
trimmed_trace = [0];
nc_range = [12];

beta_best_rec = [];
xborder_rec = [];
for i=1:numel(idx_list)
    compare_list = [idx_list(i) idx_list(i)];
    LiveFly_COMPARE_1x2x(compare_list,isBcd1X,fea_range,trimmed_trace,nc_range);
end
