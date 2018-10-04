DatasetList_={'B6_all','H6B6_all','B6C6H6_all','WT_CP'}; posborder = [-17.2 -15 -15 -5];
offset=[240 220 200 270];
pspotmax=[37 80 70 70]/100;
%pspotmax=[1 1 1 1];
%DatasetList_={'B6_all_trim','H6B6_all_trim','B6C6H6_all_trim','WT_trim'};
%DatasetList_={'B6_CP','B6_AB'};
%DatasetList_={'H6B6_CP','H6B6_AB'};
%DatasetList_={'B6C6H6_CP','B6C6H6_AB'};
load('feature_label.mat');
% Feature to plot
fea_range=[1 5 8 9 10];
nc_range=[12];

cnt=0;
wd=12;
for nc=nc_range
    for fea=fea_range
        cnt=cnt+1;
        
        for i=1:numel(DatasetList_)
            load(['tmp/' DatasetList_{i}]);
            % Plot mean curve with error
            figure(11);
            subplot(numel(nc_range),numel(fea_range),cnt);
            errorbar(pos_range,mf_rec{fea,nc},sf_rec{fea,nc}./sqrt(nf_rec{fea,nc}),'Display',DatasetList_{i});
            hold on;
            % Plot mean curve without error
            figure(12);
            subplot(numel(nc_range),numel(fea_range),cnt);
            plot(pos_range,mf_rec{fea,nc},'Display',DatasetList_{i},'LineWidth',2);
            hold on;
            % Plot heatmap (relative time)
            if fea==fea_range(1)
                figure(13);
                % Find position id in range
                posidx=(heatmapI(nc-8).pos_range-posborder(find(nc_range==nc),i)-wd/2).*(heatmapI(nc-8).pos_range-posborder(find(nc_range==nc),i)+wd/2)<0;
                tmp=mean(heatmapI(nc-8).Rel_map(:,posidx),2);
                tmp_=heatmapI(nc-8).Rel_time-offset(find(nc_range==nc),i);
                % Draw
                plot(tmp_,tmp/pspotmax(i),'Display',DatasetList_{i}); hold on;
            end
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

figure(13);
xlabel('time (s)');
ylabel('P_{Spot}');
legend('Show');