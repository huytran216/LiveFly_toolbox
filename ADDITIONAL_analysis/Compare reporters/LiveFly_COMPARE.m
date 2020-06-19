%% Re-extract all features in all strains before running the script
% THe data should be in Data_analysis/
warning off;
%% Params
fld='../../Data analysis/';
load([fld 'feature_label.mat']);

Nsample=1;      % Minimum total nuclei per position
Nsample_indi=1; % Minimum nuclei per embryo per position

plot_embryo_error=1;    % plot error based on embryo diversity (1) or nuclei error (merged, 0)
shaded_error_bar = 1;   % Plot shaded errorbar or normal errorbar

trimmed_trace = 1;      % trimmed trace (1) or not (0)
smooth_curve = 1;
compare_1x2x = false;
%% Set up data list

dtset = struct('filename','','label','');
dtset(1).filename = 'hb-vk33';  dtset(1).label = 'hb-P2';

dtset(2).filename = 'B6-near';  dtset(2).label = 'B6';
dtset(3).filename = 'B9-near';  dtset(3).label = 'B9';
dtset(4).filename = 'B12-near'; dtset(4).label = 'B12';
dtset(5).filename = 'B6-far';   dtset(5).label = 'B6-far';
dtset(6).filename = 'B9-far';   dtset(6).label = 'B9-far';
dtset(7).filename = 'B12-far';  dtset(7).label = 'B12-far';

dtset(8).filename = 'hb-II';  dtset(8).label = 'hb-P2 ii';
dtset(9).filename = 'hb-III-Lucas2018';  dtset(9).label = 'hb-P2 iii';

dtset(10).filename = 'H6B6-near';   dtset(10).label = 'H6B6';

dtset(11).filename = 'Z6';  dtset(11).label = 'Z6';
dtset(12).filename = 'Z2B6-near';  dtset(12).label = 'Z2B6';
dtset(13).filename = 'Z7B6-near';  dtset(13).label = 'Z7B6';

%compare_list = [1 2 5 3 7];isBcd1X = [0 0 0 0 0]; % For B6-B9-B12 comparison
%compare_list = [1 2 10]; isBcd1X = [0 0 0]; % For hb-B6-H6B6 comparison
compare_list = [1 2 3 7 1 2 3 7];isBcd1X=[0 0 0 0 1 1 1 1]; % For hb-B6-H6B6 comparison, 1x2x
%compare_list = [1 8 9]; isBcd1X =[0 0 0 ];% For vk33 vs random insertion
%compare_list = [7 7];isBcd1X=[0 1];
compare_list = [1 2 3 4]; isBcd1X = compare_list*0;

% Set folder containing mean data (contain dash)
folder={};
folder{1}='tmp/';
folder{2}='tmp_trimmed/';

avr = [600 750 1100];                % Mean nc13 duration
%% Cook label_list
DatasetLabel = {dtset(compare_list).label};
DatasetFile = {dtset(compare_list).filename};
for i=1:numel(compare_list)
    if isBcd1X(i)
        DatasetLabel{i}=[DatasetLabel{i} '-Bcd1X'];
        DatasetFile{i}=[DatasetFile{i} '-Bcd1X'];
        compare_1x2x = true;
    end
end
%% Feature to plot, plot settings
fea_range=[1];
nc_range=[13];

%AP_limit = [-35 20]; % for B6-B9-B12
AP_limit = [-32 20]; % for zld
%% Plot stuffs
ymax=zeros(1,16);
ymin=zeros(1,16)+1e5;
rec_fea=[];

h=[];
hplot = [];
xborder = []; hborder = []; vborder = []; noembryo = [];
cnt=0;
cnt2=0;
for nc=nc_range
    for fea=fea_range
        cnt=cnt+1;
        for i=1:numel(compare_list)
            cnt2=cnt2+1;
            load(fullfile(fld,folder{trimmed_trace+1},DatasetFile{i}),'pos_range','mf_rec','nf_rec','nf_indi','mf_indi','sf_indi','FitRes');
            % Get confidence interval of inferreable params
                tsfirst = find(FitRes(nc-8).xborder_rec(fea,:),1,'first');
                noembryo(i,nc) = sum(FitRes(nc-8).vborder_rec(fea,:)~=0);
                xborder(i,fea,nc,1) = FitRes(nc-8).xborder_rec(fea,tsfirst);
                xborder(i,fea,nc,2) = FitRes(nc-8).xborder_CI(fea,1);
                xborder(i,fea,nc,3) = FitRes(nc-8).xborder_CI(fea,2);
                hborder(i,fea,nc,1) = FitRes(nc-8).hborder_rec(fea,tsfirst);
                hborder(i,fea,nc,2) = FitRes(nc-8).hborder_CI(fea,1);
                hborder(i,fea,nc,3) = FitRes(nc-8).hborder_CI(fea,2);
                vborder(i,fea,nc,1) = FitRes(nc-8).vborder_rec(fea,tsfirst)*2;
                vborder(i,fea,nc,2) = FitRes(nc-8).vborder_CI(fea,1)*2;
                vborder(i,fea,nc,3) = FitRes(nc-8).vborder_CI(fea,2)*2;
            % Get mean and standard deviation curve
                ne_rec = zeros(size(pos_range));
                mf_indi_ = zeros(size(pos_range));
                sf_indi_ = zeros(size(pos_range));
            % Calculate number of embryo per position
                for j=1:size(nf_indi,3)
                    if numel(nf_indi{1,nc,j})
                        ne_rec=ne_rec+(nf_indi{1,nc,j}>Nsample_indi);
                    end
                end
            % Calculate mean of indi curve:
                for j=1:size(nf_indi,3)
                    if numel(nf_indi{1,nc,j})
                        tmp=mf_indi{fea,nc,j};
                        tmp(isnan(tmp))=0;
                        nf_indi{1,nc,j}(isnan(nf_indi{1,nc,j}))=0;
                        mf_indi_=mf_indi_ + tmp.*(nf_indi{1,nc,j}>Nsample_indi);
                        sf_indi_=sf_indi_ + (tmp.*(nf_indi{1,nc,j}>Nsample_indi)).^2;
                    end
                end
                mf_indi_=mf_indi_./ne_rec;
                sf_indi_=sqrt(sf_indi_./ne_rec - mf_indi_.^2);
            % Plot mean curve with error
            
                figure(40);
                subplot(numel(nc_range),numel(fea_range),cnt);
                color_order = get(gca,'colororder');
                set(gcf,'Position',[300   10   400*numel(fea_range)   250*numel(nc_range)]);
                flttmp = nf_rec{fea,nc}>Nsample;
                if plot_embryo_error
                    mtmp = mf_indi_(flttmp);
                    stmp = sf_indi_(flttmp)./sqrt(ne_rec(flttmp)-1);
                else
                    mtmp = mf_rec{fea,nc}(flttmp);
                    stmp = sf_rec{fea,nc}(flttmp)./sqrt(nf_rec{fea,nc}(flttmp)-1);
                end
                if smooth_curve
                    stmp_tmp = isnan(stmp);
                    mtmp = smooth(mtmp(end:-1:1));
                    mtmp = mtmp(end:-1:1);
                    %stmp = smooth(stmp(end:-1:1));
                    %stmp = stmp(end:-1:1);
                    %stmp(stmp_tmp)=nan;
                    mtmp(stmp_tmp)=nan;
                end
                color = corder(compare_list(i));
                if shaded_error_bar
                    if compare_1x2x & numel(compare_list)==2
                        if cnt2==1
                            color=corder(2);
                        else
                            color=corder(4);
                        end
                    end
                    h40(cnt2)=shadedErrorBar(pos_range(flttmp),mtmp,stmp,{'color',color,'Display',DatasetLabel{i}},0.8,1);
                else
                    h40(cnt2)=errorbar(pos_range(flttmp),mtmp,stmp,'Display',DatasetLabel{i});
                end
                hold on;
            % Plot mean curve with error, normalized by maximum value
                figure(41);
                subplot(numel(nc_range),numel(fea_range),cnt);
                set(gcf,'Position',[300   10   400*numel(fea_range)   250*numel(nc_range)]);
                if shaded_error_bar
                    if compare_1x2x & numel(compare_list)==2
                        if cnt2==1
                            color=corder(2);
                        else
                            color=corder(4);
                        end
                    end
                    h41(cnt2)=shadedErrorBar(pos_range(flttmp),mtmp/vborder(i,fea,nc,2),stmp/vborder(i,fea,nc,2),{'color',color,'Display',dtset(i).label},0.8,1);
                else
                    h41(cnt2)=errorbar(pos_range(flttmp),mtmp/vborder(i,fea,nc,2),stmp/vborder(i,fea,nc,2),'Display',dtset(i).label);
                end
                hold on;
            % Plot mean curve without error
                figure(42);
                subplot(numel(nc_range),numel(fea_range),cnt);
                set(gcf,'Position',[300   10   400*numel(fea_range)   250*numel(nc_range)]);
                plot(pos_range(flttmp),mtmp,'Display',DatasetLabel{i},'LineWidth',2);
                hold on;
        end
        % Get legends and correct axis
            figure(40);
            subplot(numel(nc_range),numel(fea_range),cnt);
            if nc==nc_range(end)
                xlabel('AP axis (%EL)');
            end
            ylabel(feature_label{fea});
            % get similar ylim for all nc
            yax = get(gca,'Ylim');
            ymax(fea)=max(ymax(fea),yax(2));
            ymin(fea)=min(ymin(fea),yax(1));
            rec_fea(cnt) = fea;
            
            figure(41);
            subplot(numel(nc_range),numel(fea_range),cnt);
            if nc==nc_range(end)
                xlabel('AP axis (%EL)');
            end
            ylabel(['Normalized' feature_label{fea}]);
            
            figure(42);
            subplot(numel(nc_range),numel(fea_range),cnt);
            if nc==nc_range(end)
                xlabel('AP axis (%EL)');
            end
            ylabel(feature_label{fea});
            if cnt==1
                legend('Show');
            end
            xlim(AP_limit);ylim([0 1.3]);
        % Plot fitted params (Boundary position,Hill,maximum expression) with confidence error bar:
            if compare_1x2x
                barplot_fun = @barplot_bound_1x2x;
            else
                barplot_fun = @barplot_bound;
            end
            
            figure(43);            
            subplot(numel(nc_range),numel(fea_range),cnt);
            set(gcf,'Position',[900   10   400*numel(fea_range)   250*numel(nc_range)]);
            btmp=barplot_fun(xborder(:,fea,nc,1),xborder(:,fea,nc,2),xborder(:,fea,nc,3),DatasetLabel,AP_limit(1),'h');
            box on;
            xlim(AP_limit);
            xlabel('AP axis (%EL)');
            xtickangle(0);
            
            figure(44);
            subplot(numel(nc_range),numel(fea_range),cnt);
            set(gcf,'Position',[900   10   400*numel(fea_range)   250*numel(nc_range)]);
            btmp=barplot_fun(hborder(:,fea,nc,1),hborder(:,fea,nc,2),hborder(:,fea,nc,3),DatasetLabel,0);
            box on;
            xtickangle(0);
            ylabel('Hill coefficient');
            
            figure(45);
            v1=-log(1/19)./hborder(:,fea,nc,1)*20*2;
            v2=-log(1/19)./hborder(:,fea,nc,2)*20*2;
            v3=-log(1/19)./hborder(:,fea,nc,3)*20*2;
            subplot(numel(nc_range),numel(fea_range),cnt);
            set(gcf,'Position',[900   10   400*numel(fea_range)   250*numel(nc_range)]);
            btmp=barplot_fun(v1,v2,v3,DatasetLabel,0);
            box on;
            xtickangle(0);
            ylabel('Boundary width (%EL)');
            
            figure(46);
            subplot(numel(nc_range),numel(fea_range),cnt);
            set(gcf,'Position',[900   10   400*numel(fea_range)   250*numel(nc_range)]);
            btmp=barplot_fun(vborder(:,fea,nc,1),vborder(:,fea,nc,2),vborder(:,fea,nc,3),DatasetLabel,0);
            box on;
            xtickangle(0);
            ylabel('Maximum expression');
            % Full boundary position & Width
            
            figure(47);
            subplot(numel(nc_range),numel(fea_range),cnt);
            set(gcf,'Position',[900   10   400*numel(fea_range)   250*numel(nc_range)]);
            if compare_1x2x
                btmp=barplot_bound_width_1x2x(xborder(:,fea,nc,1),xborder(:,fea,nc,2),xborder(:,fea,nc,3),v1,v2,v3,DatasetLabel,-50,'h');
            else
                btmp=barplot_bound_width(xborder(:,fea,nc,1),xborder(:,fea,nc,2),xborder(:,fea,nc,3),v1,v2,v3,DatasetLabel,-50,'h');
            end
            box on;
            xtickangle(45);
            xlabel('Boundary position (%EL)');
        % Make a table out of it: constructs are put vertically
            figure(30);
            header={};
            header{1,1}='FEATURE:';
            header{1,2}=feature_label{fea};
            header{2,1}='Construct';
            header{2,2}='N';
            header{2,3}='xborder';
            header{2,4}='Hill';
            header{2,5}='Width';
            header{2,6}='HalfMax';
            cnti=2; % Line number
            for i=1:numel(compare_list)
                cnti=cnti+1;
                header{cnti,1}=DatasetLabel{i};
                header{cnti,2}=noembryo(i,nc);
                header{cnti,3}=[num2str(xborder(i,fea,nc,1),'%.1f') ' ' char(177) ' ' ...
                    num2str((xborder(i,fea,nc,3) - xborder(i,fea,nc,2))/2,'%.1f')];
                header{cnti,4}=[num2str(hborder(i,fea,nc,1),'%.1f') ' ' char(177) ' ' ...
                    num2str((hborder(i,fea,nc,3) - hborder(i,fea,nc,2))/2,'%.1f')];
                v=-log(1/19)./hborder(i,fea,nc,:)*20*2;
                header{cnti,5}=[num2str(v(1),'%.1f') ' ' char(177) ' ' ...
                    num2str((v(3) - v(2))/2,'%.1f')];
                switch fea
                    case 9
                        header{cnti,6}=[num2str(vborder(i,fea,nc,1)/100,'%.1f') ' ' char(177) ' ' ...
                            num2str((vborder(i,fea,nc,3) - vborder(i,fea,nc,2))/100,'%.1f')];
                    case 5
                        header{cnti,6}=[num2str(vborder(i,fea,nc,1),'%.2f') ' ' char(177) ' ' ...
                            num2str((vborder(i,fea,nc,3) - vborder(i,fea,nc,2)),'%.2f')];
                    case 4
                        header{cnti,6}=[num2str(vborder(i,fea,nc,1),'%.2f') ' ' char(177) ' ' ...
                            num2str((vborder(i,fea,nc,3) - vborder(i,fea,nc,2)),'%.2f')];
                    case 1
                        header{cnti,6}=[num2str(vborder(i,fea,nc,1),'%.2f') ' ' char(177) ' ' ...
                            num2str((vborder(i,fea,nc,3) - vborder(i,fea,nc,2)),'%.2f')];    
                end
            end
            % Make the table:
            htmp=figure(30);set(htmp,'Name',['nc ' num2str(nc),', Feature: ' feature_label{fea}],'Position',[100 100 800 500]);
            htmptable =   uitable('Parent',htmp,'Position',[0 0 800 500]);
            htmptable.Data = header;
            set(htmptable,'ColumnEditable',true(1,10))
    end
end

% Apply similar axis to all subplots in Fig 11 and Fig 12
for i=1:cnt
    figure(40);
    subplot(numel(nc_range),numel(fea_range),i);
    xlim(AP_limit);
    if rec_fea(i)==1
        ylim([0 1.1])
    else
        ylim([ymin(rec_fea(i)) ymax(rec_fea(i))]);
    end
    h=legend(h40(1:numel({dtset(compare_list).label})),DatasetLabel);
    legend boxoff
    set(gca,'XTick',round((AP_limit(1):10:AP_limit(2))/5)*5);    

    figure(42);
    subplot(numel(nc_range),numel(fea_range),i);
    xlim(AP_limit);
    if rec_fea(i)==1
        ylim([0 1.1])
    else
        ylim([ymin(rec_fea(i)) ymax(rec_fea(i))]);
    end
    set(gca,'XTick',round((AP_limit(1):5:AP_limit(2))/5)*5);
end
%% Extract interphase sort by nuclear cycles
tphase_all=cell(1,15);
ymax_ = zeros(1,15);
ymin_ = ones(1,15)*1e5;

time_limit{11} = [0 600];
time_limit{12} = [0 800];
time_limit{13} = [0 1100];
hall=cell(1,numel(compare_list));
tall=cell(1,numel(compare_list));
avr_step = [0 cumsum(avr(nc_range-10))];
for i=1:numel(compare_list)
    load(fullfile(fld,folder{trimmed_trace+1},DatasetFile{i}),'heatmapI','pos_range');
    hall{i}=[];
    tall{i}=[];
    for nc=nc_range
        ymax_(compare_list(i)) = max(ymax_(i),max(heatmapI(nc-8).Rel_map{1,2}(:)));
        tphase_all{nc} = [tphase_all{nc} heatmapI(nc-8).tphase];
        hall{i} = [hall{i};heatmapI(nc-8).Rel_map{1,2}];
        tall{i} = [tall{i} avr_step(nc - nc_range(1) + 1)+heatmapI(nc-8).Rel_time/heatmapI(nc-8).Rel_time(end)*avr(nc-10)];
    end
end

for i=1:numel(compare_list)
    load(fullfile(fld,folder{trimmed_trace+1},DatasetFile{i}),'heatmapI','pos_range');
    h=figure(i);set(h,'Name',['Kymograph ' DatasetLabel{i}]);
    for nc=nc_range        
%       Plot into subplot
%         subplot(numel(nc_range),1,find(nc_range==nc));
%         HeatMap_(heatmapI(nc-8).Rel_map{1,2}*0,pos_range,1.5*heatmapI(nc-8).Rel_time/heatmapI(nc-8).Rel_time(end)*avr(nc-10),[0 ymax_(i)]);
%         hold on;
%         HeatMap_(heatmapI(nc-8).Rel_map{1,2},pos_range,heatmapI(nc-8).Rel_time/heatmapI(nc-8).Rel_time(end)*avr(nc-10),[0 ymax_(i)]);
%         tphase_all{nc} = [tphase_all{nc} heatmapI(nc-8).tphase];
%         xlim(AP_limit);
%         ylim(time_limit{nc});
%         caxis([0 ymax_(compare_list(i))]);
%         xlabel('AP axis (%EL)');
%         ylabel('Time (s)');
    end
    %       Plot into same plot
        HeatMap_(hall{i}*0,pos_range,tall{i}*1.5,[0 ymax_(i)]);
        hold on;
        HeatMap_(hall{i},pos_range,tall{i},[0 ymax_(i)]);
        xlim(AP_limit);
        ylim([0 sum(avr(nc_range-10))]);
        caxis([0 ymax_(compare_list(i))]);
        xlabel('AP axis (%EL)');
        ylabel('Time (s)');
        for nc=nc_range(1:end-1)
            plot3([-50 50],[avr_step(nc-9) avr_step(nc-9)],[1 1],'LineStyle','--','color',[1 1 1],'LineWidth',2);
        end
        set(gcf,'Position',[600   10   400   250*numel(nc_range)]);
end
%% Extract interphase sort by reporters
ymax_ = zeros(1,15);
ymin_ = ones(1,15)*1e5;

time_limit{11} = [0 600];
time_limit{12} = [0 800];
time_limit{13} = [0 1100];

for i=1:numel(compare_list)
    load(fullfile(fld,folder{trimmed_trace+1},DatasetFile{i}),'heatmapI','pos_range');
    for nc=nc_range
        ymax_(nc) = max(ymax_(nc),max(heatmapI(nc-8).Rel_map{1,2}(:)));
    end
end

for i=1:numel(compare_list)
    load(fullfile(fld,folder{trimmed_trace+1},DatasetFile{i}),'heatmapI','pos_range');
    for nc=nc_range
        h=figure(nc+10);set(h,'Name',['Kymograph nc' num2str(nc)]);
        subplot(numel(compare_list),1,i);
        HeatMap_(heatmapI(nc-8).Rel_map{1,2}*0,pos_range,1.5*heatmapI(nc-8).Rel_time/heatmapI(nc-8).Rel_time(end)*avr(nc-10),[0 ymax_(i)]);
        hold on;
        HeatMap_(heatmapI(nc-8).Rel_map{1,2},pos_range,heatmapI(nc-8).Rel_time/heatmapI(nc-8).Rel_time(end)*avr(nc-10),[0 ymax_(i)]);
        xlim(AP_limit);
        ylim(time_limit{nc});
        caxis([0 ymax_(nc)]);
        xlabel('AP axis (%EL)');
        ylabel('Time (s)');
        set(gcf,'Position',[1200   10   400   250*numel(nc_range)]);
    end
end