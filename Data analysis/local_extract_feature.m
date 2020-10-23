function [idrec,tsrec,xrec,yrec,fearec,xborder_rec,hborder_rec,wborder_rec,vborder_rec,xborder_CI,hborder_CI,wborder_CI,vborder_CI,xaxis_all_full,yaxis_all_full,fearec_all_full]=local_extract_feature(datamat,mf_indi,cycleno,ts_spec,normalize_feature,feature_label,AP,fitoption,skipfit)
%% Skipfit: skip fitting or not
    if ~exist('skipfit','var')
        skipfit = false;
    end
%% Define features to be extracted: Check feature_label.mat for reference
    fea_range=1:numel(feature_label);
%% Initiating the feature storage
    idrec=[];
    fearec=cell(1,max(fea_range));
    xrec=[];
    yrec=[];
    tsrec=[];
%% Extract the feature values for border detection
    Nsample_min = 3;    % number of nuclei required per bin to calculate mean properly.
    
    isplot=0;
    if isplot
        figure('Name',['nc' num2str(cycleno)]);
    end
    xborder_rec=[];         % Position of the border
    wborder_rec=[];         % Width of the border
    vborder_rec=[];         % Value of the feature at the border
    hborder_rec=[];         % Hill function for the border
    xborder_CI=[]; % CI for xborder
    hborder_CI=[]; % CI for hborder
    wborder_CI=[]; % CI for wborder
    vborder_CI=[]; % CI for vborder
    
    fearec_all={};
    xaxis_all={};
    yaxis_all={};
    
    fearec_all_full={};
    xaxis_all_full={};
    yaxis_all_full={};
    
    cnt=0;
    mfearec_all=cell(max(fea_range),max(ts_spec));
    mxaxis_all=cell(max(fea_range),max(ts_spec));
    
    for feaidx=fea_range
        for tsidx=ts_spec
            % Take mean curves
            if fitoption(4)
                mxaxis_all{feaidx,tsidx}=AP(1):AP(2);
                mfearec_all{feaidx,tsidx}=mf_indi{feaidx,cycleno,tsidx};
                tmp=isnan(mf_indi{feaidx,cycleno,tsidx});
                mxaxis_all{feaidx,tsidx}=mxaxis_all{feaidx,tsidx}(~tmp);
                mfearec_all{feaidx,tsidx}=mfearec_all{feaidx,tsidx}(~tmp);
            end
            % Take individual time point
                % For specific AP window (for inference only):
                    % Take cell id
                    idselect=find(([datamat(:).cycle]==cycleno)&(([datamat(:).tscnt]==tsidx))...
                        &([datamat(:).x]*100-50)>AP(1)&([datamat(:).x]*100-50)<AP(2));
                    % Take cell position and feature value
                    xaxis=[datamat(idselect).x]*100-50; % No embryo alignment yet
                    yaxis=[datamat(idselect).y]*100; % No embryo alignment yet
                    tmp=arrayfun(@(x) subindex(datamat(x).Feature,feaidx),idselect);
                    % Record all the features
                    fearec_all{feaidx,tsidx}=tmp(tmp>=0);
                    xaxis_all{feaidx,tsidx}=xaxis(tmp>=0);
                    yaxis_all{feaidx,tsidx}=yaxis(tmp>=0);
                % For all AP window (for showing only, not for inference)
                    % Take cell id
                    idselect=find(([datamat(:).cycle]==cycleno)&(([datamat(:).tscnt]==tsidx)));
                    % Take cell position and feature value
                    xaxis=[datamat(idselect).x]*100-50; % No embryo alignment yet
                    yaxis=[datamat(idselect).y]*100; % No embryo alignment yet
                    tmp=arrayfun(@(x) subindex(datamat(x).Feature,feaidx),idselect);
                    % Record all the features
                    fearec_all_full{feaidx,tsidx}=tmp(tmp>=0);
                    xaxis_all_full{feaidx,tsidx}=xaxis(tmp>=0);
                    yaxis_all_full{feaidx,tsidx}=yaxis(tmp>=0);
        end
        % 
        cnt=cnt+1;
        if isplot
            subplot(3,3,cnt);
        end
        % Find the border region by fitting with Hill function
        limit=0;
        if numel(xaxis_all)
            fitoption_=fitoption;
            %if feaidx==1
            %    fitoption_(1)=1;
            %end
            % Fit only features with sigmoid like patterns
            % if (~any([1 4 5 7 8 9 13 14 16 20]==feaidx))||(skipfit)
            if (skipfit)
                % Find interphase duration
                for ts=ts_spec
                    % Default value
                    xborder_rec(feaidx,ts)=1e10;
                    hborder_rec(feaidx,ts)=1e10;
                    wborder_rec(feaidx,ts)=0;
                    vborder_rec(feaidx,ts)=mode(fearec_all{feaidx,ts})/2;
                    % Default CI
                    xborder_CI(feaidx,1:2)=[1e10 1e10];
                    hborder_CI(feaidx,1:2)=[1e10 1e10];
                    wborder_CI(feaidx,1:2)=[0 0];
                    vborder_CI(feaidx,1:2)=ones(1,2)*mode(fearec_all{feaidx,ts})/2;
                end
            else
                erroroption_=[1 1 1 1];
                if ~fitoption(4)
                    [xborder,hborder,wborder,vborder,CIxborder,CIhborder,CIwborder,CIvborder]=detect_border_all(xaxis_all(feaidx,ts_spec),fearec_all(feaidx,ts_spec),limit,isplot,fitoption_,erroroption_);
                else
                    [xborder,hborder,wborder,vborder,CIxborder,CIhborder,CIwborder,CIvborder]=detect_border_all(mxaxis_all(feaidx,ts_spec),mfearec_all(feaidx,ts_spec),limit,isplot,fitoption_,erroroption_);
                end
                % Record the border information
                xborder_rec(feaidx,ts_spec)=xborder';
                hborder_rec(feaidx,ts_spec)=hborder';
                wborder_rec(feaidx,ts_spec)=wborder';
                vborder_rec(feaidx,ts_spec)=vborder';
                % Record the CI
                xborder_CI(feaidx,1:2)=CIxborder;
                hborder_CI(feaidx,1:2)=CIhborder;
                wborder_CI(feaidx,1:2)=CIwborder;
                vborder_CI(feaidx,1:2)=CIvborder;
            end            
        end
    end
 %% Save the data with normalized feature values - for all embryo, no specific AP window
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