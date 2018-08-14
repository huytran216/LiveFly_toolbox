DatasetList_={'B6_Huy','H6B6_Huy','B6C6H6_Huy','WT_Huy','WT_Zelda'};
DatasetList_={'B6_trim','H6B6_trim','B6C6H6_trim','WT_trim','WT_Zelda_trim'};
load('feature_label.mat');
% Feature to plot
fea_range=[1 5 8 9 10];
nc_range=[11 12 13];

cnt=0;

for nc=nc_range
    for fea=fea_range
        cnt=cnt+1;
        
        for i=1:numel(DatasetList_)
            load(['tmp_/' DatasetList_{i}]);
            figure(11);
            subplot(numel(nc_range),numel(fea_range),cnt);
            errorbar(pos_range,mf_rec{fea,nc},sf_rec{fea,nc}./sqrt(nf_rec{fea,nc}),'Display',DatasetList_{i});
            figure(12);
            subplot(numel(nc_range),numel(fea_range),cnt);
            plot(pos_range,mf_rec{fea,nc},'Display',DatasetList_{i},'LineWidth',2);
            hold on;
        end
        figure(11);
        xlabel('AP axis (%EL)');
        ylabel(feature_label{fea});
        if cnt==1
            legend('Show');
        end
        figure(12);
        xlabel('AP axis (%EL)');
        ylabel(feature_label{fea});
        if cnt==1
            legend('Show');
        end
    end
end

    