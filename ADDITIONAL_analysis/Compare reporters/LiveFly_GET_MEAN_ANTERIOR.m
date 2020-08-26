%% Re-extract all features in all strains before running the script
% The data should be in Data_analysis/
%% Params
addpath '..\..\Tool\hmm_dont_edit';
addpath '..\..\Tool\hmm_dont_edit\utilities\';
fld='../../Data analysis/final_dataset'; % Location of the dataset

%% Select data
%% Set up data list

dtset = struct('filename','','label','','pos_SS',[],'time_oSS',[]);
dtset(1).filename = 'hb-vk33';  dtset(1).label = 'hb-P2'; dtset(1).pos_SS=[-32 -27];dtset(1).time_oSS=[0 800]; dtset(1).pos_boundary = -15;

dtset(2).filename = 'B6-near';  dtset(2).label = 'B6'; dtset(2).pos_SS=[-32 -27];dtset(2).time_oSS=[0 850]; dtset(2).pos_boundary = -15;
dtset(3).filename = 'B9-near';  dtset(3).label = 'B9'; dtset(3).pos_SS=[-32 -27];dtset(3).time_oSS=[0 800]; dtset(3).pos_boundary = -15;
dtset(4).filename = 'B12-near'; dtset(4).label = 'B12'; dtset(4).pos_SS=[-32 -27];dtset(4).time_oSS=[0 850]; dtset(4).pos_boundary = -15;
dtset(5).filename = 'B6-far';   dtset(5).label = 'B6-far'; dtset(5).pos_SS=[-35 -28];dtset(5).time_oSS=[0 800];
dtset(6).filename = 'B9-far';   dtset(6).label = 'B9-far'; dtset(6).pos_SS=[];dtset(6).time_oSS=[];
dtset(7).filename = 'B12-far';  dtset(7).label = 'B12-far'; dtset(7).pos_SS=[-32 -20];dtset(7).time_oSS=[0 800];

dtset(8).filename = 'hb-II';  dtset(8).label = 'hb-ii';  dtset(8).pos_SS=[-32 -25];dtset(8).time_oSS=[0 800];dtset(8).pos_boundary = -8.4;
dtset(9).filename = 'hb-III-Lucas2018';  dtset(9).label = 'hb-iii';  dtset(9).pos_SS=[-32 -27];dtset(9).time_oSS=[0 800];dtset(9).pos_boundary = -10;

dtset(10).filename = 'H6B6-near';   dtset(10).label = 'H6B6';  dtset(10).pos_SS=[-32 -25];dtset(10).time_oSS=[0 800];dtset(10).pos_boundary = -15;

dtset(11).filename = 'Z6';  dtset(11).label = 'Z6'; dtset(11).pos_SS=[-32 32];dtset(11).time_oSS=[0 800];
dtset(12).filename = 'Z2B6-near';  dtset(12).label = 'Z2B6'; dtset(12).pos_SS=[-32 -22];dtset(12).time_oSS=[0 800];dtset(12).pos_boundary = -15;
dtset(13).filename = 'Z7B6-near';  dtset(13).label = 'Z7B6'; dtset(13).pos_SS=[-32 -22];dtset(13).time_oSS=[0 800];

compare_list =  [1 2 3];                 % For B6-B9-B12 comparison

isBcd1X =    zeros(size(compare_list));  % 1 if load Bcd1x , 0 if not

nc_range = [11 12 13];                % Interphase duration
avr = [600 750 1100];                 % Mean nc13 duration

check_boundary = 0;                   % Scan at the anterior at the boundary
    dw = 5; % Set boundary width for analysis of time to reach boundary.
plot_intensity = 1;                   % 1 for intensity, 0 for pspot.
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
            ylb = 'Intensity per nucleus';
        case 2
            ylb = 'Intensity per spot';
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
                    if plot_intensity==0
                        tr = tr>0;
                    end
                    if plot_intensity==2
                        tr(tr==0) = nan;
                    end
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
                    
                    trace_all{total}=tr;
                    time_all{total}=sum(~isnan(tr));
                    tmax=max(time_all{total},tmax);
                    
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
            for j=1:numel(Irec)
                if numel(Irec{j})>3
                    mIrec(j)=nanmean(Irec{j});
                    sIrec(j)=sqrt(nanvar(Irec{j}));
                    nIrec(j)=sum(~isnan(Irec{j}));
                else
                    mIrec(j)=0;
                    sIrec(j)=0;
                    nIrec(j)=0;
                end
            end
        end
        tax{cycleno-8}=time_ax(1:numel(mIrec));
        mIax{cycleno-8}=mIrec;
        sIax{cycleno-8}=sIrec;
        nIax{cycleno-8}=nIrec;
        nI(cycleno-8)=total;
        % Trim traces till mitosis
        last_time = find(nIrec<=max(nIrec)/3,1,'first');
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
    %% Save the data
    mkdir('data');
    save(['data/mean_' Dataset],'tax','mIax','sIax','nIax','trace_length','dt','Dataset','nI','time_first');
    %% Plot summary on mean and CV between Dataset in all cell cycle
    tlinear = [];
    tnormalized = [];
    tlinear_prev=0;
    tnormalized_prev = 0;
    for cycleno=nc_range
        tnormalized = [tnormalized tnormalized_prev+tax{cycleno-8}/tax{cycleno-8}(end)*avr(cycleno-10)];
        tnormalized_prev = tnormalized(end);
        tlinear = [tlinear tlinear_prev+tax{cycleno-8}];
        tlinear_prev = tlinear(end);
    end
    maxt = numel(cell2mat(tax));
    
%     figure;
%     shadedErrorBar(tlinear,cell2mat(mIax),cell2mat(sIax),{'LineWidth',2,'color','k'},0.5);
%     xlabel('Time (s)');
%     ylabel(ylb);
%     title(Dataset);
    figure(19);
    
    h19{i}=shadedErrorBar(tnormalized,cell2mat(mIax),cell2mat(sIax)./sqrt(nI(cycleno-8)),{'Display',DatasetLabel{i},'color',corder(compare_list(i)),'LineWidth',2},0.5);hold on;
    xlabel('Time (s)');
    ylabel(ylb);
    % Record for comparison:
    mI_rec_tmp{i}=cell2mat(mIax);
    tI_rec_tmp{i}=tnormalized;
    
    figure(20);
    subplot(211);
    plot(tlinear,cell2mat(sIax)./cell2mat(mIax),'Marker','.','Display',DatasetLabel{i});hold on;
    xlabel('Time (s)');
    ylabel('CV');
    subplot(212);
    plot(tnormalized,cell2mat(sIax)./cell2mat(mIax),'Marker','.','Display',DatasetLabel{i});hold on;
    xlabel('time %');
    ylabel('CV');
    % Compare time to reach steady stades between Dataset in each nc
    figure(30+i);
    for cycleno=nc_range
        sItmp = sIax{cycleno-8}./sqrt(nIax{cycleno-8});
        mItmp = mIax{cycleno-8};
        h30i{cycleno}=shadedErrorBar(tax{cycleno-8}/max(tax{cycleno-8})*avr(cycleno-10),mItmp,sItmp,{'Display',['nc' num2str(cycleno)],'color',corder(cycleno)},0.5);hold on;
    end
    xlabel('normalized time');
    ylabel(['Normalized ' ylb]);
    title(DatasetLabel{i});
    legend boxoff
    figure(21);
    for cycleno=nc_range
        if numel(nc_range)>1
            subplot(numel(nc_range),1,find(nc_range==cycleno));
        end
        sItmp = sIax{cycleno-8}./sqrt(nIax{cycleno-8});
        mItmp = mIax{cycleno-8};
        h21{i}=shadedErrorBar(tax{cycleno-8}/max(tax{cycleno-8})*avr(cycleno-10),mItmp,sItmp,{'Display',DatasetLabel{i},'color',corder(compare_list(i))},0.5);hold on;
        xlabel('Time (s)');
        ylabel([ylb]);
        %plot(tax{cycleno-8}/max(tax{cycleno-8}),mIax{cycleno-8}/max(mIax{cycleno-8}),'Display',Dataset);hold on;
    end
    
    figure(22);
    for cycleno=nc_range
        if numel(nc_range)>1
            subplot(numel(nc_range),1,find(nc_range==cycleno));
        end
        peak_range = 1:50;
        maxmI(i,cycleno-10) = max(mIax{cycleno-8}(peak_range));
        sItmp = sIax{cycleno-8}/maxmI(i,cycleno-10)/sqrt(nI(cycleno-8));
        mItmp = mIax{cycleno-8}/maxmI(i,cycleno-10);
        h22{i}=shadedErrorBar(tax{cycleno-8}/max(tax{cycleno-8})*avr(cycleno-10),mItmp,sItmp,{'Display',DatasetLabel{i},'color',corder(compare_list(i))},0.5);hold on;
        %plot(tax{cycleno-8}/max(tax{cycleno-8}),mIax{cycleno-8}/max(mIax{cycleno-8}),'Display',Dataset);hold on;
        xlabel('Time (s)');
        ylabel(['normalized ' ylb]);
    end
    
    % Create legend for comparing between nuclear cycle
    figure(30+i);
    h30i_=[];
    for cycleno=nc_range
        h30i_=[h30i_ h30i{cycleno}.mainLine];
        leg{cycleno}=['nc' num2str(cycleno)];
    end
    legend(h30i_,leg(nc_range));
    legend boxoff;
    %% Get mean traces to get kon, koff
    cycleno=13;
    time_SS=dtset(compare_list(i)).time_oSS;
    time_ax = tax{cycleno-8}/max(tax{cycleno-8})*avr(cycleno-10);
    Irec = mIax{cycleno-8};
    srec = sIax{cycleno-8};
    newt = (time_ax>time_SS(1))&(time_ax<time_SS(2));
    time_ax = time_ax(newt);
    time_ax=time_ax-time_ax(1);
    Irec = Irec(newt);
    srec = srec(newt);
    % Record for showing plot:
    Irec_full{i}=mIax{cycleno-8};
    err_rec_full{i}=sIax{cycleno-8}/sqrt(nI(cycleno-8)-1);
    time_rec_full{i}=tax{cycleno-8};
    % Record for fitting - mean and error
    Irec_{i}=Irec;
    err_rec_{i}=srec/sqrt(nI(cycleno-8)-1);
    time_rec_{i}=time_ax;
    % Record for individual traces:
    Irec_indi{i}=[];
    for k=1:total
        if max(find(newt))<=numel(trace_all{k})
            Irec_indi{i} = [Irec_indi{i};trace_all{k}(find(newt))];
        end
    end
    toffset(i) = time_SS(1);
end
%% Plot legend:
figure(21);
h21_ = [];
for i = 1:numel(compare_list)
    h21_ = [h21_ h21{i}.mainLine];
end
legend(h21_,DatasetLabel);
legend boxoff
figure(19);
h19_ = [];
for i = 1:numel(compare_list)
    h19_ = [h19_ h19{i}.mainLine];
end
legend(h19_,DatasetLabel,'Location','NorthWest');
legend boxoff
set(gcf,'Position',[362   741   818   237]);
figure(22);
h22_ = [];
for i = 1:numel(compare_list)
    h22_ = [h22_ h22{i}.mainLine];
end
%legend(h22_,DatasetLabel);
%legend boxoff
%% Plot peak intensity:
figure;bar(maxmI');
set(gca,'XTickLabel',leg(nc_range));
legend(DatasetLabel);
legend boxoff;
%% Expression ratios:
if numel(compare_list)==2
    tmp = min(numel(mI_rec_tmp{1}),numel(mI_rec_tmp{2}));
    figure;plot(tI_rec_tmp{1}(1:tmp),mI_rec_tmp{2}(1:tmp)./mI_rec_tmp{1}(1:tmp)); ylim([0 5]);
    hold on; plot(tI_rec_tmp{1},(tI_rec_tmp{1}*0+1)*1.5,'--')
    xlabel('time (s)');
    ylabel('Mean expression ratio');
end
%% Plot extracted params for T0 first spot appearance
figure(10);
% Plot prediction:
x = [100:10:500];
plot(x,normcdf(x,310,50),'LineStyle','--','LineWidth',1,'color','k','Display','N(305 s, 50 s)');
legend show;
legend('Location','NorthWest');
legend boxoff;
xlim([100 500]);
figure(11);
cnt = 0;
for cycleno=nc_range
    cnt = cnt + 1;
    subplot(numel(nc_range),2,cnt*2-1);
    barplot_bound(mT0(cycleno-10,:),mT0_lb(cycleno-10,:),mT0_ub(cycleno-10,:),DatasetLabel);
    ylabel('Mean T0');
    subplot(numel(nc_range),2,cnt*2);
    barplot_bound(sT0(cycleno-10,:),sT0_lb(cycleno-10,:),sT0_ub(cycleno-10,:),DatasetLabel);
    ylabel('Std T0');
end

%% setup the elongation time
load('theRightL');
kelongation=40;
L=ms(1:kelongation*dt:end);
L=L/sum(L);
sigmaT0=50/dt;
gaussfilt = normpdf(-sigmaT0*2:sigmaT0*2,0,sigmaT0);
gaussfilt = gaussfilt/sum(gaussfilt);
L = conv2(gaussfilt,L,'full');

%% Guessing the peak value:
x0=[];
nsample =0;
for i=1:numel(compare_list)
    x0 = [x0 Irec_{i}(end)];
    nsample = nsample + numel(Irec_indi{i});
end
x0 = [x0 240*ones(1,numel(compare_list))];
%% Fit all to one: - fit kon, koff
fun_1 = @(x,t) x(3)*(x(1)/abs(x(1)+x(2)) + x(2)/abs(x(1)+x(2))*exp(-(t-x(4))*(x(1)+x(2)))).*(t>=x(4)); % kon, koff, Imax, toffset
   
% All independent params 
    % Begin fitting
    xcrit1 = 1:numel(compare_list)*4;
    xord = 1:numel(xcrit1);
    fmin1 = 1e20;
    
    for cnt=1:20
        [x1_,fmin1_] = fminsearchbnd(@(x) fitall_(x,xcrit1,Irec_indi(:),time_rec_(:),fun_1,L),...
            [rand(1,2*numel(compare_list))*0.005 x0],...
            [repmat([1 1],1,numel(compare_list))*0.0001 x0*0],...
            [repmat([1 1],1,numel(compare_list)) x0*100]);
        if fmin1>fmin1_
            fmin1=fmin1_;
            x1=x1_;
        end
    end
    % Find confidence interval:
    ll1 = get_loglikelihood(fmin1,nsample);
    
% Similar kon:
    xcrit2 = 1:numel(compare_list)*4;
    for i=1:numel(compare_list)
        xcrit2(i*2-1)=1;
        xcrit2(i*2)=i+1;
        xcrit2(end-2*numel(compare_list)+i) = 1+numel(compare_list)+i;
        xcrit2(end-numel(compare_list)+i) = 1+2*numel(compare_list)+i;
    end
    fmin2 = 1e20;
    for cnt=1:20
        [x2_,fmin2_] = fminsearchbnd(@(x) fitall_(x,xcrit2,Irec_indi(:),time_rec_(:),fun_1,L),...
            [rand(1,1+numel(compare_list)).*0.01 x0],...
            [ones(1,1+numel(compare_list))*0.0001 x0*0],...
            [ones(1,1+numel(compare_list))*1 x0*100]);
        if fmin2>fmin2_
            fmin2=fmin2_;
            x2=x2_;
        end
    end
    % Find confidence interval:
    ll2 = get_loglikelihood(fmin2,nsample);
    
% Similar koff:
    xcrit3 = 1:numel(compare_list)*4;
    for i=1:numel(compare_list)
        xcrit3(i*2-1)=i;
        xcrit3(i*2)=numel(compare_list)+1;
        xcrit3(end-2*numel(compare_list)+i) = 1+numel(compare_list)+i;
        xcrit3(end-numel(compare_list)+i) = 1+2*numel(compare_list)+i;
    end
    fmin3 = 1e20;
    for cnt=1:20
        [x3_,fmin3_] = fminsearchbnd(@(x) fitall_(x,xcrit3,Irec_indi(:),time_rec_(:),fun_1,L),...
            [rand(1,1+numel(compare_list)).*0.01 x0],...
            [ones(1,1+numel(compare_list))*0.0001 x0*0],...
            [ones(1,1+numel(compare_list))*1 x0*100]);
        if fmin3>fmin3_
            fmin3=fmin3_;
            x3=x3_;
        end
    end
    % Find confidence interval:
    ll3 = get_loglikelihood(fmin3,nsample);
    

display(['Independent fitting']);
    show_info(x1(xcrit1));
    display(['Score: ' num2str(log(fmin1))]);
 display(['Keeping kon']);
    show_info(x2(xcrit2));
    display(['Score: ' num2str(log(fmin2))]);
 display(['Keeping koff']);
    show_info(x3(xcrit3));
    display(['Score: ' num2str(log(fmin3))]);
%% All independent params - fit ON, tau
    
% fun_2 = @(x,t) x(3)*(x(1) + (1-x(1))*exp(-t/x(2))); % PON, tau, Imax
%     % Begin fitting
%     xcrit1 = 1:numel(compare_list)*3;
%     fmin1_2 = 1e20;
% 
%     for cnt=1:10
%         [x1_,fmin1_] = fminsearchbnd(@(x) fitall_(x,xcrit1,Irec_indi(compare_list),time_rec_(compare_list),fun_2),...
%             [rand(1,2*numel(compare_list)).*repmat([1 150],1,numel(compare_list)) x0],...
%             [repmat([0.1 5],1,numel(compare_list)) x0*0],...
%             [repmat([1 200],1,numel(compare_list)) x0*3]);
%         if fmin1_2>fmin1_
%             fmin1_2=fmin1_;
%             x1_2=x1_;
%         end
%     end
%     % Find confidence interval:
%     ll1 = get_loglikelihood(fmin1_2,nsample);
%% Plot individual fits:
figure;
cnt= 0;
[~,Iout]=fitall_(x3,xcrit3,Irec_indi(:),time_rec_(:),fun_1,L);
for Datasetidx = 1:numel(compare_list)
    cnt = cnt+1;
    subplot(numel(compare_list),1,cnt);
    htmp = shadedErrorBar(time_rec_full{Datasetidx}/max(time_rec_full{Datasetidx})*avr(cycleno-10),Irec_full{Datasetidx},err_rec_full{Datasetidx},{'Display','data','LineStyle','-','color','k'},0.7);
    hold on;
    ax = [dtset(compare_list(Datasetidx)).time_oSS(1):dt:dtset(compare_list(Datasetidx)).time_oSS(2)];
    hline = plot(dtset(compare_list(Datasetidx)).time_oSS(1)+time_rec_{Datasetidx},Iout{cnt},'Display','model','LineStyle','--','color','k','LineWidth',2);
    ylabel(ylb);
    xlabel('Time (s)');
    yl = get(gca,'ylim');
    ylim([0 yl(2)]);
    xlim([0 1200]);
    legend([htmp.mainLine hline],{'data','fitted'},'box','off','Location','NorthWest');
end