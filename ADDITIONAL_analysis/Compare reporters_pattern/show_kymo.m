%% Prepare for fit of sigmoid function:
if fit_boundary
    if fit_sigmoid
        ft = fittype( '1/(1+exp(b*(x-c)))', 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.StartPoint = [0.485375648722841 0];
        opts.Lower = [0 -40];
        opts.Upper = [1e10 30];
    else
        % Find the maximum level at the plateau
        switch kymo_intensity
            case 0
                fea = 16;   % Mean intensity
            case 1
                fea = 21;   % PSpot
            case 2
                fea=1;
            end
    end
end


%% Extract interphase sort by nuclear cycles
h195= figure(195);
[ha195]=tight_subplot(1, numel(compare_list), 0.02); % Set small distance between plots

tphase_all=cell(1,15);
ymax_ = zeros(1,15);            % Maximum intensity For individual reporters
ymin_ = ones(1,15)*1e5;
ymax__ = zeros(1,15);           % Maximum intensity For each nuclear cycle
ymin__ = ones(1,15)*1e5;

time_limit{11} = [0 600];
time_limit{12} = [0 800];
time_limit{13} = [0 1100];
hall=cell(1,numel(compare_list));
sall=cell(1,numel(compare_list));
tall=cell(1,numel(compare_list));
ncall=cell(1,numel(compare_list));
xall=cell(1,numel(compare_list));
xall_ub=cell(1,numel(compare_list));
xall_lb=cell(1,numel(compare_list));
avr_step = [0 cumsum(avr(nc_range-10))];
hborder = [];
xborder = [];
vborder = [];

hrec1 = [];
hrec2 = [];
for i=1:numel(compare_list)
    if kymo_intensity<2 % At steady state
        load(fullfile(fld,folder{2},DatasetFile{i}),'heatmapI','pos_range','FitRes');
    else % taking whole trace
        load(fullfile(fld,folder{1},DatasetFile{i}),'heatmapI','pos_range','FitRes');
    end
    %if fit_sigmoid
        % Extract peaked intensity
    for nc = nc_range
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
    sall{i}=[];
    tall{i}=[];
    ncall{i} =[];
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
        tmptime = heatmapI(nc-8).Rel_time/heatmapI(nc-8).Rel_time(end)*avr(nc-10);
        switch kymo_intensity
            case 0
                tmpmap = heatmapI(nc-8).Rel_map{1,1};   % Intensity
                stmpmap = heatmapI(nc-8).Rel_map_std{1,1};  % Error
            case 1
                tmpmap = heatmapI(nc-8).Rel_map{1,2};   % PSpot
                stmpmap = heatmapI(nc-8).Rel_map_std{1,2};  % Error
            case 2
                tmpmap = heatmapI(nc-8).Rel_map{2,1};   % ON
                stmpmap = heatmapI(nc-8).Rel_map_std{2,1};  % Error
                % Trim only to expression window for this feature
                tmpmap(tmptime>avr_cut(nc-10),:) = 0;
                stmpmap(tmptime>avr_cut(nc-10),:) = 0;
        end
        hall{i} = [hall{i};tmpmap];
        sall{i} = [sall{i};stmpmap];
        tall{i} = [tall{i} avr_step(nc - nc_range(1) + 1)+tmptime];
        ncall{i} = [ncall{i} tmptime*0+nc];
        if kymo_intensity==1
            ymax_(compare_list(i)) = max(ymax_(i),max(heatmapI(nc-8).Rel_map{1,1+kymo_intensity}(:)));
            ymax__(nc) = max(ymax__(nc),max(heatmapI(nc-8).Rel_map{1,2}(:)));
        else
            ymax_(compare_list(i)) = 1;
            ymax__(nc) = 1;
        end
        % Plot curve at specific position
        if plot_boundary
            bd_pos = xborder(i,fea,13)+[-2 -1 0 1 2];
        else
            bd_pos = [-35:-27];
        end
        [~,~,bd_idx] = intersect(round(bd_pos),pos_range);
        figure(130);
        subplot(2,ceil(numel(compare_list)/2),i);
        expression_at_positionsum=nanmean(tmpmap(:,bd_idx),2);
        expression_at_positionsum(tmptime>avr_cut(nc-10))=0;
        plot(tmptime,expression_at_positionsum);
        hold on;
        set(gcf,'Position',[680   568   305   246]);
        xlabel('Time (s)');
        ylim([0 1]);
    end
    % Find the boundary position at specific frame
    for tcnt = 1:numel(tall{i})
        htmp = hall{i}(tcnt,:);
        idselect = (pos_range(:)>=-30)&(pos_range(:)<=20);
        xall{i}(tcnt) = -50;
        xall_ub{i}(tcnt) = -50;
        xall_lb{i}(tcnt) = -50;
        aall{i}(tcnt) = 0;
        aall_ub{i}(tcnt) = 0;
        aall_lb{i}(tcnt) = 0;
        
        maxtmp = vborder(i,fea,nc_range(end),1)/2;
        if fit_boundary
            if ~fit_sigmoid
                if maxtmp>0
                    % find anterior most inflection points:
                    tmp =  find(((htmp(end:-1:2)-maxtmp)<0)&((htmp(end-1:-1:1)-maxtmp)>0),1,'first');
                    if ~numel(tmp)
                        xall{i}(tcnt) = -50;
                    else
                        xall{i}(tcnt) = pos_range(numel(htmp)-tmp);
                    end
                end
            else
                [fitresult, gof] = fit( pos_range(idselect)', htmp(idselect)'/maxtmp/2, ft, opts );
                [fitresult]
                if (fitresult.c>-30)
                    xall{i}(tcnt) = fitresult.c;
                    ci = confint(fitresult);
                    xall_ub{i}(tcnt) = ci(2,2);
                    xall_lb{i}(tcnt) = ci(1,2);
                    aall{i}(tcnt) = maxtmp*2;
                    ci = confint(fitresult);
                    aall_ub{i}(tcnt) = maxtmp*2;
                    aall_lb{i}(tcnt) = maxtmp*2;
                else
                    xall{i}(tcnt) = -50;
                    xall_ub{i}(tcnt) = -50;
                    xall_lb{i}(tcnt) = -50;
                    aall{i}(tcnt) = 0;
                    aall_ub{i}(tcnt) = 0;
                    aall_lb{i}(tcnt) = 0;
                end
                [tcnt xall{i}(tcnt)]
            end
        end
    end
    % Plot boundary position over time
    if fit_boundary
        if fit_sigmoid
            fborder = figure(190);
        else
            fborder = figure(191);
        end
        %subplot(2,1,1);
        color = corder(compare_list(i));
        hrec1(i) = shadedErrorBar(tall{i}',xall{i}',[xall{i}-xall_lb{i};xall_ub{i}-xall{i}]',...
            {'color',color,'LineWidth',2},0.8,1); hold on;
        % Get the time to reach final decision (+-1%)
            treach = find(xall{i}>max(xall{i}-2),1,'first');
            plot([tall{i}(treach) tall{i}(treach)],[-30 xall{i}(treach)],...
                'LineStyle','--','color',color,'LineWidth',1);
        xlim([200 800]);
        ylim([-30 10]);
    %     subplot(2,1,2);
    %     color = corder(compare_list(i));
    %     shadedErrorBar(tall{i}',aall{i}',[aall{i}-aall_lb{i};xall_ub{i}-aall{i}]',...
    %         {'color',color,'Display',DatasetLabel{i},'LineWidth',2},0.8,1); hold on;
    %     ylim([0 1]);
        % Check over lap
        figure(195);
        axes(ha195(i));
        for nc=nc_range
            ttmp = ncall{i}==nc;
            tax = tall{i}(ttmp)'-min(tall{i}(ttmp)');
            hrec2(i,nc)=shadedErrorBar(tax,xall{i}(ttmp)',[xall{i}(ttmp)-xall_lb{i}(ttmp);xall_ub{i}(ttmp)-xall{i}(ttmp)]',...
            {'LineWidth',2},0.8,1); hold on;            
        end
        xlim([200 800]);
        ylim([-30 10]);
    end
end
if fit_boundary
    figure(fborder);
        legend(hrec1,DatasetLabel);
        ylabel('Boundary position');
        xlabel('Time (s)');
   figure(195);
        ncleg = {};
        for nc=nc_range
            ncleg{end+1} = num2str(nc);
        end
        for i=1:numel(compare_list)
            axes(ha195(i));
            set(gca,'XTick',[],'YTick',[]);
            %legend(hrec2(i,nc_range),ncleg); 
            %ylabel('Boundary position');
            %xlabel('Time (s)');
        end
end
%% Show kymograph of reporters horizontally
h= figure(1);
[ha]=tight_subplot(1, numel(compare_list), [.01 .03],[.1 .03],[.1 .1]); % Set small distance between plots
for i=1:numel(compare_list)
    load(fullfile(fld,folder{1},DatasetFile{i}),'heatmapI','pos_range');
    %h=figure(i);set(h,'Name',['Kymograph ' DatasetLabel{i}]);
    figure(1);
    axes(ha(i));
    
    % Set cut and mitosis timing
    nc_draw = cumsum(avr(nc_range-10));
    nc_cut = [0 nc_draw(1:end-1)] + avr_cut(nc_range-10);
    
    % Plot into same plot
        HeatMap_(hall{i}*0,pos_range,tall{i}*1.5,[0 max(ymax__(:))]);
        hold on;
        HeatMap_(double(hall{i}),pos_range,tall{i},[0 max(ymax__(:))]);
        if kymo_intensity==1
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
        for j=1:numel(nc_cut)
            plot3([-50 50],[nc_draw(j) nc_draw(j)],[1 1],'LineStyle','--','color',[1 1 1],'LineWidth',2);
        end
        set(gcf,'Position',[500   100   250*numel(compare_list)   250*numel(nc_range)]);
end
%% Plot reporters vertically by nc (to compare e.g. boundary position)
if plot_vertically
    for i=1:numel(compare_list)
        load(fullfile(fld,folder{1},DatasetFile{i}),'heatmapI','pos_range');
        for nc=nc_range
            h=figure(nc+10);set(h,'Name',['Kymograph nc' num2str(nc)]);
            subplot(numel(compare_list),1,i);
            HeatMap_(heatmapI(nc-8).Rel_map{1,1}*0,pos_range,1.5*heatmapI(nc-8).Rel_time/heatmapI(nc-8).Rel_time(end)*avr(nc-10),[0 ymax__(i)]);
            hold on;
            switch kymo_intensity
                case 0
                    HeatMap_(heatmapI(nc-8).Rel_map{1,1+kymo_intensity},pos_range,heatmapI(nc-8).Rel_time/heatmapI(nc-8).Rel_time(end)*avr(nc-10),[0 ymax__(i)]);
                case 1
                    HeatMap_(heatmapI(nc-8).Rel_map{1,1+kymo_intensity},pos_range,heatmapI(nc-8).Rel_time/heatmapI(nc-8).Rel_time(end)*avr(nc-10),[0 ymax__(i)]);
                case 2
                    HeatMap_(heatmapI(nc-8).Rel_map{2,1},pos_range,heatmapI(nc-8).Rel_time/heatmapI(nc-8).Rel_time(end)*avr(nc-10),[0 ymax__(i)]);
            end
            xlim(AP_limit);
            ylim(time_limit{nc});
            caxis([0 ymax__(nc)]);
            xlabel('AP axis (%EL)');
            ylabel('Time (s)');
            set(gcf,'Position',[500   400   400   250*numel(nc_range)]);
        end
    end
end