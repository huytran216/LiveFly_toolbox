function []=local_show_feature(cycleno,ts_spec,DatasetList,DatasetFeature,AP,feature_label,fea_plot,showoption_plot,fea_ref,defcolor)

%% Show data:
 %% Extract all the feature values per embryo 
if numel(ts_spec)
    figure('Name',['nc' num2str(cycleno)]);
    cnt=0;
    yrange={};
    timefea=[2 3 4 5 6];
    for tsidx=ts_spec
        for feaidx=fea_plot(:)'
            if numel(DatasetFeature.fearec_all{feaidx,tsidx})
                cnt=cnt+1;
                % Check for border position if needed
                if fea_ref & DatasetList.nc_ref & showoption_plot(1)
                    xborder=DatasetList.xborder_rec(fea_ref,tsidx)
                end
                % Begin plotting                
                subplot(numel(ts_spec),numel(fea_plot),cnt)
                if showoption_plot(6)||(~any(timefea==feaidx))
                    plot(DatasetFeature.xaxis_all{feaidx,tsidx},DatasetFeature.fearec_all{feaidx,tsidx},'.b');hold on;
                else
                    plot(DatasetFeature.xaxis_all{feaidx,tsidx},DatasetFeature.fearec_all{feaidx,tsidx}...
                        .*DatasetFeature.fearec_all{10,tsidx},'.b');hold on;
                end
                tmpx=[AP(1):AP(2)];
                % Show fitted Hill curve if needed
                if showoption_plot(2)
                    plot(tmpx,2*DatasetFeature.vborder_rec(feaidx,tsidx)*sigmf(tmpx,[ -DatasetFeature.hborder_rec(feaidx,tsidx)*0.04 DatasetFeature.xborder_rec(feaidx,tsidx)]),'--k');
                end
                % Set axis limit
                xlim(AP);
                yrange{feaidx}=[];
                if tsidx==ts_spec(1)
                    yrange{feaidx}=get(gca,'YLim');
                    title(feature_label{feaidx});
                else
                    if numel(yrange{feaidx})
                        ylim(yrange{feaidx});
                    else
                        yrange{feaidx}=get(gca,'YLim');
                        title(feature_label{feaidx});
                    end
                end
                if feaidx==fea_plot(1)
                    ylabel(['embryo ' num2str(tsidx)]);
                end
            end
        end    
    end
end