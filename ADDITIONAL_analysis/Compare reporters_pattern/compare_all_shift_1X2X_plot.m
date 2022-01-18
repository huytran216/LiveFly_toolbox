idx_list = [2 3 12 1 2 3 12 1];
isBcdE1 = [0 0 0 0 1 1 1 1];
nc_list = [13 13 13 13 13 13 13 13];

%idx_list = [1 2 3 4 10 12];
%isBcdE1 = [1 1 1 1 1 1];
%nc_list = [13 13 13 13 13 13];

%idx_list = [2 12 2 12];
%isBcdE1 = [0 0 1 1];
%nc_list = [13 13 13 13];

time_range =0; % Set to 0 if using whole trace
if time_range==0
    suffix = '';
else
    suffix = num2str(time_range);
end
beta_best_rec = [];
sbeta_best_rec = [];
xborder_rec = [];
for i=1:numel(idx_list)
    if isBcdE1(i)
        isBcd1X = [0 1];
    else
        isBcd1X = [0 2];
    end
    nc_range = nc_list(i);
    filename = [num2str([idx_list(i) idx_list(i)]) '_' num2str(isBcd1X) '_nc' num2str(nc_range) '_' suffix '.mat'];
    dattmp=load(['shift_rec/' filename]);
    beta_best_rec(i) = dattmp.beta_best_*log(2);
    sbeta_best_rec(i) = dattmp.sbeta_best_*log(2);
    xborder_rec(i,:) = dattmp.xborder(1,dattmp.fea,dattmp.nc,1:3);
end
beta_best_rec = round(beta_best_rec/0.5)*0.5;
sbeta_best_rec = round(sbeta_best_rec/0.5)*0.5;

%close all;
% Make label:
DatasetLabel = {dattmp.dtset(idx_list).label};
figure;
subplot(121);
bar(beta_best_rec,'FaceColor',[1 1 1]);hold on;
errorbar(1:numel(idx_list),beta_best_rec,sbeta_best_rec,'LineStyle','none');
set(gca,'XTickLabel',DatasetLabel);
ylim([3 14]);
subplot(122);
errorbar(xborder_rec(:,1),beta_best_rec,sbeta_best_rec,sbeta_best_rec,(xborder_rec(:,3)-xborder_rec(:,2))/2,(xborder_rec(:,3)-xborder_rec(:,2))/2,'Marker','.','LineStyle','none');
ylim([3 14]);
%% Make specific figure compare between Ebcd1 and dBcd
x={'dbcd';'bcdE1';};

y = reshape(beta_best_rec,[numel(beta_best_rec)/2 2])';
errorplus = reshape(sbeta_best_rec,[numel(beta_best_rec)/2 2]);
errorminus=errorplus;

figure;
bar(categorical(x),y);
hBar = bar(y, 0.8);
ctr = [];
ydt = [];
for k1 = 1:size(y,2)
    ctr(k1,:) = bsxfun(@plus, hBar(k1).XData, hBar(k1).XOffset');     
    ydt(k1,:) = hBar(k1).YData;                    
end
hold on
errorbar(ctr, ydt, errorplus, 'LineStyle','none','LineWidth',1,'color','k');
hold off
ylim([0 15]);
set(gca,'XTickLabel',x)