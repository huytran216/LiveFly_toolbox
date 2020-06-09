cnt=0;
cycleno=13;
figure;
isIntensity = 1;
for i=1:size(heatmapI(5).Rel_map_indi,1)
    if numel(heatmapI(5).Rel_map_indi{i,1,2})
        cnt=cnt+1;
        subplot(3,3,cnt);
        HeatMap_(heatmapI(cycleno-8).Rel_map_indi{i,1,isIntensity+1},1:size(heatmapI(cycleno-8).Rel_map_indi{i,1,isIntensity+1},2),1:size(heatmapI(cycleno-8).Rel_map_indi{i,1,isIntensity+1},1),[0 15*isIntensity+1])
    end
end