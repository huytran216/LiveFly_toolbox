%% Re-extract all features in all strains before running the script
% The data should be in Data_analysis/
%% Params
addpath '..\..\Tool\hmm_dont_edit';
addpath '..\..\Tool\hmm_dont_edit\utilities\';
addpath '../Compare reporters_pattern';
fld='../../Data analysis/final_dataset'; % Location of the dataset

%% Select data
%% Set up data list

dtset = struct('filename','','label','','pos_SS',[],'time_oSS',[]);
dtset(1).filename = 'hb-vk33';  dtset(1).label = 'hb-P2'; dtset(1).pos_SS=[-32 -27];dtset(1).time_oSS=[0 800]; dtset(1).pos_boundary = -8.5;

dtset(2).filename = 'B6-near';  dtset(2).label = 'B6'; dtset(2).pos_SS=[-32 -27];dtset(2).time_oSS=[0 850]; dtset(2).pos_boundary = -14;
dtset(3).filename = 'B9-near';  dtset(3).label = 'B9'; dtset(3).pos_SS=[-32 -27];dtset(3).time_oSS=[0 800]; dtset(3).pos_boundary = -7;
dtset(4).filename = 'B12-near'; dtset(4).label = 'B12'; dtset(4).pos_SS=[-32 -27];dtset(4).time_oSS=[0 850]; dtset(4).pos_boundary = -5;
dtset(5).filename = 'B6-far';   dtset(5).label = 'B6-far'; dtset(5).pos_SS=[-35 -28];dtset(5).time_oSS=[0 800];
dtset(6).filename = 'B9-far';   dtset(6).label = 'B9-far'; dtset(6).pos_SS=[];dtset(6).time_oSS=[];
dtset(7).filename = 'B12-far';  dtset(7).label = 'B12-far'; dtset(7).pos_SS=[-32 -20];dtset(7).time_oSS=[0 800];

dtset(8).filename = 'hb-II';  dtset(8).label = 'rand II';  dtset(8).pos_SS=[-32 -25];dtset(8).time_oSS=[0 800];dtset(8).pos_boundary = -8.4;
dtset(9).filename = 'hb-III-Lucas2018';  dtset(9).label = 'rand. III';  dtset(9).pos_SS=[-32 -27];dtset(9).time_oSS=[0 800];dtset(9).pos_boundary = -10;

dtset(10).filename = 'H6B6-near';   dtset(10).label = 'H6B6';  dtset(10).pos_SS=[-32 -25];dtset(10).time_oSS=[0 800];dtset(10).pos_boundary = -15;

dtset(11).filename = 'Z6';  dtset(11).label = 'Z6'; dtset(11).pos_SS=[-20 20];dtset(11).time_oSS=[0 800];
dtset(12).filename = 'Z2B6-near';  dtset(12).label = 'Z2B6'; dtset(12).pos_SS=[-32 -25];dtset(12).time_oSS=[0 800];dtset(12).pos_boundary = -15;
dtset(13).filename = 'Z7B6-near';  dtset(13).label = 'Z7B6'; dtset(13).pos_SS=[-32 -22];dtset(13).time_oSS=[0 800];

compare_list =  [10];                 % For B6-B9-B12 comparison

isBcd1X =    zeros(size(compare_list));  % 1 if load Bcd1x , 0 if not

nc_range = [13];                % Interphase duration
avr = [600 750 1100];                 % Mean nc13 duration

check_boundary = 0;                   % Scan at the anterior at the boundary
    dw = 5; % Set boundary width for analysis of time to reach boundary.
plot_intensity =1;                   % 0 for pspot, 1 for loci intensity, 2 for spot intensity

only_ON = 1;                          % Apply only to nuclei with ON signals
%% Cook label_list
DatasetLabel = {dtset(compare_list).label};
DatasetFile = {dtset(compare_list).filename};
for i=1:numel(compare_list)
    if isBcd1X(i)
        DatasetLabel{i}=[DatasetLabel{i} '-Bcd1X'];
        DatasetFile{i}=[DatasetFile{i} '-Bcd1X'];
        compare_1x2x = true;
    end
    if check_boundary
        dtset(compare_list(i)).pos_SS = dtset(compare_list(i)).pos_boundary + [-dw +dw];
    end
end
%% MS2 configuration and memory length
dt=10;
scale_time = 1;                     % Scaling interphase by mean interphase duration
    switch plot_intensity
        case 0
            ylb = 'PSpot';
        case 1
            ylb = 'Loci Intensity';
        case 2
            ylb = 'Spot Intensity';
    end
    switch only_ON
        case 0 
            on_label = '';
        case 1 
            on_label = ', ON';
    end
%% Begin analyzing
Irec_indi = {};
h21=[];
maxmI = [];
for i = 1:numel(compare_list)
    Dataset = DatasetFile{i};
    % Load the dataset
    load(fullfile(fld,Dataset),'Nmov','datamat','DatasetFeature','DatasetList','heatmapI');
    % Define aux function
    subindex = @(x,y) x(y);
    %% Refine interphase
    tinterphase=[];
    trace_length=[];
    tax={};
    mIax={};
    sIax={};
    nIax={};
    nI=[];
    time_first=[];
    %figure;
    for cycleno=nc_range
        ts_spec=find(arrayfun(@(x) getfield(DatasetList,{x},['nc' num2str(cycleno)]),1:Nmov));
        % Get interphase:
        tinterphase(cycleno) = max(heatmapI(cycleno-8).Rel_time);
        time_all={};
        trace_all={};
        time_all_aligned={};
        trace_all_aligned={};
        
        if numel(ts_spec)
            % Get interphase duration
            tphase=DatasetFeature(cycleno-8).vborder_rec(10,:)*2;
            if scale_time
                tphase_set = ones(size(tphase))*mean(tphase(tphase>0));
            else
                tphase_set = tphase;
            end
            total=0;
            tmax=0;
            cntempty=0; 
            minpositive=100000;
            for tsidx = ts_spec
                % Get nuclei id that match cycleno and time series and position
                idselect=find(([datamat(:).cycle]==cycleno)&([datamat(:).tscnt]==tsidx)...
                    &(100*([datamat(:).x]-0.5)>=dtset(compare_list(i)).pos_SS(1))&...
                    (100*([datamat(:).x]-0.5)<=dtset(compare_list(i)).pos_SS(2)));
                ison=arrayfun(@(x) subindex(datamat(x).Feature,1),idselect); % Is cell valid
                idselect = idselect(ison>=0);
                time_ax = 0:dt:tphase_set(1);
                cnt=0;
                for id=idselect
                    total=total+1;
                    cnt=cnt+1;
                    % Get traces
                    tr = interp1([datamat(id).time-datamat(id).time(1) 1e5],[datamat(id).Intensity -1e10],time_ax*tphase(tsidx)/tphase_set(tsidx));
                    tr = tr(tr>=0);
                    
                    % Record first spot appearance - error in time alignment
                    tfirst = find(tr>0,1,'first')*dt;
                    if numel(tfirst)&(tfirst<500)&(tfirst>150)
                        time_first(cycleno-8,total)=tfirst;
                        trace_all_aligned{total} = tr(tfirst:min(450,numel(tr)));
                        time_all_aligned{total}= min(450,numel(tr))-tr;
                    else
                        time_first(cycleno-8,total)=NaN;
                        if tfirst<=150
                            tr=NaN*tr;
                        end
                    end
                    
                    % Extract feature: pspot, loci or spot intensty
                    if plot_intensity==0
                        tr = double(tr>0);
                    end
                    if plot_intensity==2
                        tr(tr==0) = -1;
                    end
                    % select only ON nuclei?
                    if only_ON
                        if ~any(tr>0)
                            tr = tr*0-1;
                        end
                    end
                    % Record time and trace
                    time_all{total}=sum(~isnan(tr));
                    tmax=max(time_all{total},tmax);
                    trace_all{total}=tr;
                    
                    
                    
                    if sum(tr>0)==0
                        cntempty=cntempty+1;
                    end
                    if minpositive>min(tr(tr>0))
                        if min(tr(tr>0))>0.01
                            minpositive = min(tr(tr>0));
                        end
                    end
                    %plot(time_ax(tr>=0),tr(tr>=0));hold on;
                end
            end
            % Extract mean curve and standard deviation
            Irec=cell(1,tmax);
            for k=1:total
                for j=1:time_all{k}
                    Irec{j}=[Irec{j} trace_all{k}(j)];
                end
            end
            mIrec=[];
            sIrec=[];
            nIrec=[];
            nIrec_real=[];
            for j=1:numel(Irec)
                if numel(Irec{j})>3
                    mIrec(j)=nanmean(Irec{j}(Irec{j}>=0));
                    sIrec(j)=sqrt(nanvar(Irec{j}(Irec{j}>=0)));
                    nIrec(j)=sum(~isnan(Irec{j}(Irec{j}>=0)));
                    nIrec_real(j)=sum(~isnan(Irec{j}));
                else
                    mIrec(j)=0;
                    sIrec(j)=0;
                    nIrec(j)=0;
                    nIrec_real(j)=0;
                end
            end
        end
        tax{cycleno-8}=time_ax(1:numel(mIrec));
        mIax{cycleno-8}=mIrec;
        sIax{cycleno-8}=sIrec;
        nIax{cycleno-8}=nIrec;
        nI(cycleno-8)=total;
        % Trim traces till mitosis
        last_time = find(nIrec_real<=max(nIrec_real)/3,1,'first');
        if last_time
            trace_length(cycleno-8) = last_time;
            tax{cycleno-8}=tax{cycleno-8}(1:trace_length(cycleno-8));
            mIax{cycleno-8}=mIax{cycleno-8}(1:trace_length(cycleno-8));
            sIax{cycleno-8}=sIax{cycleno-8}(1:trace_length(cycleno-8));
            nIax{cycleno-8}=nIax{cycleno-8}(1:trace_length(cycleno-8));
        else
            trace_length(cycleno-8)=numel(mIrec);
        end
        
        
        display([num2str(total) ' traces extracted']);
        display([num2str(cntempty),' empty traces']);
        
        % Plot distribution of first spot appearance:
        figure(10);
        if numel(nc_range)>1
            subplot(1,numel(nc_range),find(nc_range==cycleno));
        end
        [tmp,ax] = ecdf(time_first(cycleno-8,:));
        
        plot(ax,tmp,'Display',DatasetLabel{i},'color',corder(compare_list(i)),'LineWidth',2);
        ylabel('CDF(T_0)');
        xlabel('Time (s)');
        hold on;
        
        % Calculate mean and std:
        pd = fitdist(time_first(cycleno-8,:)','Normal');
        pd_ci = paramci(pd);
        mT0(cycleno-10,i)=pd.mu;
        sT0(cycleno-10,i)=pd.sigma;
        mT0_ub(cycleno-10,i)=pd_ci(2,1);
        mT0_lb(cycleno-10,i)=pd_ci(1,1);
        sT0_ub(cycleno-10,i)=pd_ci(2,2);
        sT0_lb(cycleno-10,i)=pd_ci(1,2);
    end
end
%% Plot data:
figure;
cnt=0;
for i=randsample(total,total)'
    if (cnt<24)&&(any(trace_all{i}>0))
        cnt=cnt+1;
        subplot(3,8,cnt);
        plot((1:numel(trace_all{i}))*dt,trace_all{i});
        ylim([0 60]);
        xlim([0 avr(nc_range-10)]);
    end
end
set(gcf,'Position',[565   518   1208   360]);
tightfig;
%% Plot specific traces with specific ID: