function [idrec,tsrec,xrec,yrec,fearec,xborder_rec,hborder_rec,wborder_rec,vborder_rec,xaxis_all,yaxis_all,fearec_all]=local_extract_feature(datamat,cycleno,ts_spec,normalize_feature,feature_label,AP,fitoption)
%% Define features to be extracted: Check feature_label.mat for reference
    fea_range=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15];
%% Initiating the feature storage
    idrec=[];
    fearec=cell(1,max(fea_range));
    xrec=[];
    yrec=[];
    tsrec=[];
%% Extract the feature values for border detection
    isplot=0;
    if isplot
        figure('Name',['nc' num2str(cycleno)]);
    end
    xborder_rec=[];         % Position of the border
    wborder_rec=[];         % Width of the border
    vborder_rec=[];         % Value of the feature at the border
    hborder_rec=[];         % Hill function for the border
    fearec_all={};
    xaxis_all={};
    yaxis_all={};
    cnt=0;
    for feaidx=fea_range
        for tsidx=ts_spec
            % Take cell id
            idselect=find(([datamat(:).cycle]==cycleno)&(([datamat(:).tscnt]==tsidx)));
            % Take cell position and feature value
            xaxis=[datamat(idselect).x]*100-50; % No embryo alignment yet
            yaxis=[datamat(idselect).y]*100; % No embryo alignment yet
            tmp=arrayfun(@(x) subindex(datamat(x).Feature,feaidx),idselect);
            % Record all the features
            fearec_all{feaidx,tsidx}=tmp(tmp>=0);
            xaxis_all{feaidx,tsidx}=xaxis(tmp>=0);
            yaxis_all{feaidx,tsidx}=yaxis(tmp>=0);
        end
        % 
        cnt=cnt+1;
        if isplot
            subplot(3,3,cnt);
        end
        % Find the border region by fitting with Hill function
        limit=0;
        if numel(xaxis_all)
            [xborder,hborder,wborder,vborder]=detect_border_all(xaxis_all(feaidx,ts_spec),fearec_all(feaidx,ts_spec),limit,isplot,fitoption);
            % Record the border information
            xborder_rec(feaidx,ts_spec)=xborder';
            hborder_rec(feaidx,ts_spec)=hborder';
            wborder_rec(feaidx,ts_spec)=wborder';
            vborder_rec(feaidx,ts_spec)=vborder';
        end
    end
 %% Save the data with normalized feature values:
    for tsidx=ts_spec
        % Take cell id
        idselect=find(([datamat(:).cycle]==cycleno)&(([datamat(:).tscnt]==tsidx)));
        % Check if the cell is valid:
        tmp=arrayfun(@(x) subindex(datamat(x).Feature,1),idselect);        
        xaxis=[datamat(idselect).x]*100-50; % No embryo alignment yet
        % Record id, time series and nuclei absolute position
        idrec=[idrec idselect(tmp>=0)];
        tsrec=[tsrec tsidx+idselect(tmp>=0)*0];
        xrec=[xrec xaxis(tmp>=0)];
        % Take cell feature
        cnt=0;
        for feaidx=fea_range
            cnt=cnt+1;
            tmp=arrayfun(@(x) subindex(datamat(x).Feature,feaidx),idselect);
            % Record all the features (for border detection):
            if normalize_feature&&(feaidx>=6)
                fearec{feaidx}=[fearec{feaidx} tmp(tmp>=0)/vborder_rec(feaidx,tsidx)/2];
            else
                fearec{feaidx}=[fearec{feaidx} tmp(tmp>=0)];
            end
        end
    end