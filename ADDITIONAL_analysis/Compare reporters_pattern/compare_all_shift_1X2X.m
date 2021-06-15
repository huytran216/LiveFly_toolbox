idx_list = [1 2 3 4 10];
isBcd1X = [0 1];
nc_range = 12;
beta_best_rec = [];
xborder_rec = [];
for i=1:numel(idx_list)
    filename = [num2str([idx_list(i) idx_list(i)]) '_' num2str(isBcd1X) '_nc' num2str(nc_range) '.mat'];
    dattmp=load(['shift_rec/' filename]);
    beta_best_rec(i) = dattmp.beta_best_*log(2);
    xborder_rec(i,:) = dattmp.xborder(1,dattmp.fea,dattmp.nc,1:3);
end

% Make label:
DatasetLabel = {dattmp.dtset(idx_list).label};
figure;
subplot(121);
bar(beta_best_rec);
set(gca,'XTickLabel',DatasetLabel);
subplot(122);
plot(xborder_rec(:,1),beta_best_rec,'o');