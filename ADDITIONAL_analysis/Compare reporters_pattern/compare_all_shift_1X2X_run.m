idx_list = [1 3 1 3 1 3 1 3 1 3 1 3 2 2 2 4 4 4 10 10 10 12 12 12];
isBcdE1 = [0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
nc_range = [11 11 12 12 13 13 11 11 12 12 13 13 11 12 13 11 12 13 11 12 13 11 12 13];

fea_range = [1];
time_range = [0 500];   % Set to 0 if whole trace


beta_best_rec = [];
xborder_rec = [];
tic
for i=1:numel(idx_list)
    if isBcdE1(i)
        isBcd1X = [0 1];
    else
        isBcd1X = [0 2];
    end
    compare_list = [idx_list(i) idx_list(i)];
    %try
        LiveFly_COMPARE_1x2x_Kymo(compare_list,isBcd1X,fea_range,time_range,nc_range(i));
    %catch
    %end
end
toc