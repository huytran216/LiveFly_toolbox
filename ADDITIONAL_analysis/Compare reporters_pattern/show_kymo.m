%% Prepare for fit of sigmoid function:
if ~fixed_plateau
    ft = fittype( 'a/(1+exp(b*(x-c)))', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = [0.957166948242946 0.485375648722841 0];
    if kymo_intensity
        opts.Upper  = [100 10 20];
    else
        opts.Upper  = [1 10 20];
    end
    opts.Lower = [0 0 -30];
    opts.Upper = [1.0 1e10 30];
else
    % Find the maximum level at the plateau
    if kymo_intensity
        fea = 16;   % Mean intensity
    else
        fea = 21;   % PSpot
    end
end


%% Extract interphase sort by nuclear cycles
tphase_all=cell(1,15);
ymax_ = zeros(1,15);            % Maximum intensity For individual reporters
ymin_ = ones(1,15)*1e5;
ymax__ = zeros(1,15);           % Maximum intensity For each nuclear cycle
ymin__ = ones(1,15)*1e5;

time_limit{11} = [0 600];
time_limit{12} = [0 800];
time_limit{13} = [0 1100];
hall=cell(1,numel(compare_list));
tall=cell(1,numel(compare_list));
xall=cell(1,numel(compare_list));
xall_ub=cell(1,numel(compare_list));
xall_lb=cell(1,numel(compare_list));
avr_step = [0 cumsum(avr(nc_range-10))];
hborder = [];
for i=1:numel(compare_list)
    load(fullfile(fld,folder{2},DatasetFile{i}),'heatmapI','pos_range','FitRes');
    if fixed_plateau
        % Extract peaked intensity
        nc = 13;
        tsfirst = find(FitRes(nc-8).xborder_rec(fea,:),1,'first');
        xborder(i,fea,nc,1) = FitRes(nc-8).xborder_rec(fea,tsfirst);
        xborder(i,fea,nc,2) = FitRes(nc-8).xborder_CI(fea,1);
        xborder(i,fea,nc,3) = FitRes(nc-8).xborder_CI(fea,2);
        hborder(i,fea,nc,1) = FitRes(nc-8).hborder_rec(fea,tsfirst);
        hborder(i,fea,nc,2) = FitRes(nc-8).hborder_CI(fea,1);
        hborder(i,fea,nc,3) = FitRes(nc-8).hborder_CI(fea,2);
        vborder(i,fea,nc,1) = FitRes(nc-8).vborder_rec(fea,tsfirst)*2;
        vborder(i,fea,nc,2) = FitRes(nc-8).vborder_CI(fea,1)*2;
        vborder(i,fea,nc,3) = FitRes(nc-8).vborder_CI(fea,2)*2;
    end
    % holder for merged kymograph
    hall{i}=[];
    tall{i}=[];
    % holder for fitted params
    xall{i}=[];
    xall_ub{i}=[];
    xall_lb{i}=[];
    aall{i}=[];
    aall_ub{i}=[];
    aall_lb{i}=[];
    for nc=nc_range
        % Merge data from different nuclear cycle into 1 graph
        tphase_all{nc} = [tphase_all{nc} heatmapI(nc-8).tphase];
        hall{i} = [hall{i};heatmapI(nc-8).Rel_map{1,1+kymo_intensity}];
        tall{i} = [tall{i} avr_step(nc - nc_range(1) + 1)+heatmapI(nc-8).Rel_time/heatmapI(nc-8).Rel_time(end)*avr(nc-10)];
        if kymo_intensity
            ymax_(compare_list(i)) = max(ymax_(i),max(heatmapI(nc-8).Rel_map{1,1+kymo_intensity}(:)));
            ymax__(nc) = max(ymax__(nc),max(heatmapI(nc-8).Rel_map{1,2}(:)));
        else
            ymax_(compare_list(i)) = 1;
            ymax__(nc) = 1;
        end        
    end
    % Find the boundary position at specific frame
    for tcnt = 1:numel(tall{i})
        htmp = hall{i}(tcnt,:);
        idselect = (pos_range(:)>=-35)&(pos_range(:)<=20);
        xall{i}(tcnt) = -50;
        xall_ub{i}(tcnt) = -50;
        xall_lb{i}(tcnt) = -50;
        aall{i}(tcnt) = 0;
        aall_ub{i}(tcnt) = 0;
        aall_lb{i}(tcnt) = 0;
        if fixed_plateau
            maxtmp = vborder(i,fea,nc,1);
            if maxtmp>0
                % find anterior most change points:
                tmp =  find(((htmp(end:-1:2)-maxtmp)<0)&((htmp(end-1:-1:1)-maxtmp)>0),1,'first');
                if ~numel(tmp)
                    xall{i}(tcnt) = -50;
                else
                    xall{i}(tcnt) = pos_range(numel(htmp)-tmp);
                end
            end
        else            
            [fitresult, gof] = fit( pos_range(idselect)', htmp(idselect)', ft, opts );
            if (fitresult.a>1e-2)&&(fitresult.c>-35)
                xall{i}(tcnt) = fitresult.c;
                ci = confint(fitresult);
                xall_ub{i}(tcnt) = ci(2,3);
                xall_lb{i}(tcnt) = ci(1,3);
                aall{i}(tcnt) = fitresult.a;
                ci = confint(fitresult);
                aall_ub{i}(tcnt) = ci(2,1);
                aall_lb{i}(tcnt) = ci(1,1);
            else
                xall{i}(tcnt) = -50;
                xall_ub{i}(tcnt) = -50;
                xall_lb{i}(tcnt) = -50;
                aall{i}(tcnt) = 0;
                aall_ub{i}(tcnt) = 0;
                aall_lb{i}(tcnt) = 0;
            end
            [tcnt fitresult.a xall{i}(tcnt)];
        end
    end
    figure(190);
    subplot(2,1,1);
    color = corder(compare_list(i));
    shadedErrorBar(tall{i}',xall{i}',[xall{i}-xall_lb{i};xall_ub{i}-xall{i}]',...
        {'color',color,'Display',DatasetLabel{i},'LineWidth',2},0.8,1); hold on;
    ylim(AP_limit);
    subplot(2,1,2);
    color = corder(compare_list(i));
    shadedErrorBar(tall{i}',aall{i}',[aall{i}-aall_lb{i};xall_ub{i}-aall{i}]',...
        {'color',color,'Display',DatasetLabel{i},'LineWidth',2},0.8,1); hold on;
    ylim([0 1]);
end
legend show;
%% Show kymograph
h= figure(1);
[ha]=tight_subplot(1, numel(compare_list), [.01 .03],[.1 .03],[.1 .1]); % Set small distance between plots
for i=1:numel(compare_list)
    load(fullfile(fld,folder{1},DatasetFile{i}),'heatmapI','pos_range');
    %h=figure(i);set(h,'Name',['Kymograph ' DatasetLabel{i}]);
    figure(1);
    axes(ha(i));
    
    
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
        HeatMap_(hall{i}*0,pos_range,tall{i}*1.5,[0 max(ymax__(:))]);
        hold on;
        HeatMap_(hall{i},pos_range,tall{i},[0 max(ymax__(:))]);
        if kymo_intensity
            caxis([0 15])
        else
            caxis([0 1]);
        end
        if i~=numel(compare_list)
            colorbar off
        end
        if i==1
            ylabel('Time (s)');
        else
            set(gca,'YTick',[]);
        end
        xlim(AP_limit);
        ylim([0 sum(avr(nc_range-10))]);
        xlabel('AP axis (%EL)');        
        for nc=nc_range(1:end-1)
            plot3([-50 50],[avr_step(nc-9) avr_step(nc-9)],[1 1],'LineStyle','--','color',[1 1 1],'LineWidth',2);
        end
        set(gcf,'Position',[500   100   250*numel(compare_list)   250*numel(nc_range)]);
end
%% Extract interphase sort by reporters
for i=1:numel(compare_list)
    load(fullfile(fld,folder{1},DatasetFile{i}),'heatmapI','pos_range');
    for nc=nc_range
        h=figure(nc+10);set(h,'Name',['Kymograph nc' num2str(nc)]);
        subplot(numel(compare_list),1,i);
        HeatMap_(heatmapI(nc-8).Rel_map{1,1+kymo_intensity}*0,pos_range,1.5*heatmapI(nc-8).Rel_time/heatmapI(nc-8).Rel_time(end)*avr(nc-10),[0 ymax__(i)]);
        hold on;
        HeatMap_(heatmapI(nc-8).Rel_map{1,1+kymo_intensity},pos_range,heatmapI(nc-8).Rel_time/heatmapI(nc-8).Rel_time(end)*avr(nc-10),[0 ymax__(i)]);
        xlim(AP_limit);
        ylim(time_limit{nc});
        caxis([0 ymax__(nc)]);
        xlabel('AP axis (%EL)');
        ylabel('Time (s)');
        set(gcf,'Position',[500   400   400   250*numel(nc_range)]);
    end
end