DatasetList={'B6_Huy','H6B6_Huy','B6C6H6_Huy','WT_Huy'};
load('feature_label.mat');
% Feature to plot
fea_range=[1 5 8 9 10];
nc_range=[11 12 13];

cnt=0;

for nc=nc_range
    for fea=fea_range
        cnt=cnt+1;
        subplot(numel(nc_range),numel(fea_range),cnt);
        for i=1:numel(DatasetList)
            load(['tmp_/' DatasetList{i}]);
            figure(11);
            errorbar(pos_range,mf_rec{fea,nc},sf_rec{fea,nc}./sqrt(nf_rec{fea,nc}),'Display',DatasetList{i});
            figure(12);
            plot(pos_range,mf_rec{fea,nc},'Display',DatasetList{i},'LineWidth',2);
            hold on;
        end
        figure(11);
        xlabel('AP axis (%EL)');
        ylabel(feature_label{fea});
        legend('Show');
        figure(12);
        xlabel('AP axis (%EL)');
        ylabel(feature_label{fea});
        legend('Show');
    end
end

    