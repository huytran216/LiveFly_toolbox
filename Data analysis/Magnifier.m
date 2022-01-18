function [tlower,tupper,cycle_range,posborder]=Magnifier(h,heatmapI,wd,nc_range,tcut1,tcut2)
% Output: 
% tlower: lowerbound for interphase: size: Nembryo x 5 (cycle)
% tupper: upperbound for interphase: size: Nembryo x 5 (cycle)
debug=1;
%% Author
% Huy Tran
%% ideas
% Check pattern within specific time range
%% Output: range of time window to take into account
figX1 = 100;
figX2 = 1200;
figY1 = 100;
figY2 = 600;

set(h,'Position',[figX1,figY1,figX2,figY2]);
set ( gcf, 'Color', [0.7 0.7 0.7] )

%% Panel for movie info
% Interface - trim kymograph with minimum sample number

Nsample=15;
    hmessages_Nsample = uicontrol('Style','text','String','Min sample:',...
        'Position',[30,40,100,20],'FontSize',10, 'HorizontalAlignment','left',...
        'ForegroundColor','white','BackgroundColor',[0.7 0.7 0.7]);
    hNSample = uicontrol('Style','edit','String',Nsample,...
        'Position',[110,40,120,20],'FontSize',10, 'HorizontalAlignment','left',...
        'Callback',@hset_NSample);

% Interface - plot absolute or relative time
    hmessages_Time = uicontrol('Style','text','String','Time axis:',...
        'Position',[30,10,100,20],'FontSize',10, 'HorizontalAlignment','left',...
        'ForegroundColor','white','BackgroundColor',[0.7 0.7 0.7]);
    hTime = uicontrol('Style','popup','String',{'Abs PSpot','Rel PSpot',...
        'Abs time - Accum PSpot','Rel time - Accum PSpot','Abs time - Intensity','Rel time - Intensity'},...
        'Position',[110,10,150,30],'FontSize',10, 'HorizontalAlignment','left',...
        'Callback',@hset_Time);

% Plot axis: 1-5 <=> 10-14
Nembryo = 0;
% Set value
auto_trim = true;
if ~exist('tcut1','var')
    auto_trim = false;
    tcut1 = [0 0 0 0 0];
    tcut2 = 1e4+[0 0 0 0 0];
end
for axidx=1:5
    % The axes for plotting kymograph
        axmap(axidx) = axes('Position',[(axidx-1)*0.2+0.03 0.40 0.15 0.55],'XTickLabel',[],'YTickLabel',[]);
        % Set beginning time
    % Set label for trimming:
        hmessage_fromTime(axidx)=uicontrol('Style','text','String','From (s):',...
            'Position',[(axidx-1)*240+36,170,100,20],'FontSize',10, 'HorizontalAlignment','left',...
            'ForegroundColor','white','BackgroundColor',[0.7 0.7 0.7]);
        hmessage_toTime(axidx)=uicontrol('Style','text','String','To (s):',...
            'Position',[(axidx-1)*240+36,140,100,20],'FontSize',10, 'HorizontalAlignment','left',...
        'ForegroundColor','white','BackgroundColor',[0.7 0.7 0.7]);
    % Find mean interphase duration
        tmp1=10000;
        tmp2=0;
        if numel(heatmapI)>=axidx+1
            if numel(heatmapI(axidx+1).mtphase)
                tmp1=heatmapI(axidx+1).mtphase;
                if isfield(heatmapI,'posborder')
                    tmp2=heatmapI(axidx+1).posborder;
                end
                if Nembryo < numel(heatmapI(axidx+1).tphase)
                    Nembryo = numel(heatmapI(axidx+1).tphase);
                end
            end
        end
    
    % Set begin and end time value for trimming
        h_fromTime(axidx) = uicontrol('Style','edit','String',num2str(max([0 tcut1(axidx)])),...
            'Position',[(axidx-1)*240+100,170,80,20],'FontSize',10, 'HorizontalAlignment','left');        
        h_toTime(axidx) = uicontrol('Style','edit','String',num2str(min([tmp1 tcut2(axidx)])),...
            'Position',[(axidx-1)*240+100,140,80,20],'FontSize',10, 'HorizontalAlignment','left');
    % Set border position
        hmessage_posBorder(axidx)=uicontrol('Style','text','String','Border (%EL):',...
            'Position',[(axidx-1)*240+16,110,120,20],'FontSize',10, 'HorizontalAlignment','left',...
            'ForegroundColor','white','BackgroundColor',[0.7 0.7 0.7]);
        h_posBorder(axidx) = uicontrol('Style','edit','String',num2str(tmp2),...
            'Position',[(axidx-1)*240+100,110,80,20],'FontSize',10, 'HorizontalAlignment','left');
end

% Trim button:
hQuit = uicontrol('Style','pushbutton','String', 'Quit',...
        'Position',[1000,20,90,40],...
        'Callback',@hQuit_CallBack);
hTrim = uicontrol('Style','pushbutton','String', 'Trim interphase',...
        'Position',[880,20,100,40],...
        'Callback',@hTrim_CallBack);
hDraw = uicontrol('Style','pushbutton','String', 'Draw mean intensity @border',...
        'Position',[680,20,180,40],...
        'Callback',@hDraw_CallBack);
% Plot the plot
set(hTime,'Value',2);
normalize_time=2;
Draw_kymo();

% Default outputs
tlower = [];
tupper = [];
cycle_range=[];

%% Functions:
    % Set mimum sample
    function hset_NSample(~,~)
        Nsample=str2double(get(hNSample,'String'));
        Draw_kymo;
    end
    % Set time axis:
    function hset_Time(~,~)
        normalize_time=get(hTime,'Value');
        Draw_kymo;
    end

    % Plot heatmap
    function Draw_kymo()
        % Get axis bigger intensity
        axmax=0;
        for i=2:6
            if numel(heatmapI)>=i
                if numel(heatmapI(i).Abs_map)
                    for hblendcnt=[1 2]
                        axmax(i,hblendcnt,1)=1;
                        axmax(i,hblendcnt,2)=max(heatmapI(i).Abs_map{hblendcnt,2}(:));
                    end
                end
            end
        end
        axmax=max(axmax);
        % Show heatmap
        for i=2:6
            if numel(heatmapI)>=i
                if numel(heatmapI(i).Abs_map)
                    axes(axmap(i-1));
                    plot3(heatmapI(i).pos_range,heatmapI(i).pos_range*0+heatmapI(i).mtphase,heatmapI(i).pos_range*0+1,'LineWidth',2,'LineStyle','--');
                    hold on;
                    switch normalize_time
                        case 1
                            hmtmp_1=heatmapI(i).Abs_time;
                            hmtmp_2=heatmapI(i).Abs_map{1,1};
                            ax_=[0 axmax(1,1,1)];
                        case 2
                            hmtmp_1=heatmapI(i).Rel_time;
                            hmtmp_2=heatmapI(i).Rel_map{1,1};
                            ax_=[0 axmax(1,1,1)];    
                        case 3    
                            hmtmp_1=heatmapI(i).Abs_time;
                            hmtmp_2=heatmapI(i).Abs_map{2,1};
                            ax_=[0 axmax(1,2,1)];
                        case 4
                            hmtmp_1=heatmapI(i).Rel_time;
                            hmtmp_2=heatmapI(i).Rel_map{2,1};
                            ax_=[0 axmax(1,2,1)];
                        case 5    
                            hmtmp_1=heatmapI(i).Abs_time;
                            hmtmp_2=heatmapI(i).Abs_map{1,2};
                            ax_=[0 axmax(1,1,2)];
                        case 6
                            hmtmp_1=heatmapI(i).Rel_time;
                            hmtmp_2=heatmapI(i).Rel_map{1,2};
                            ax_=[0 axmax(1,1,2)];
                    end
                    iscensored = heatmapI(i).Rel_count(1,:)<=Nsample;
                    hmtmp_2(:,iscensored)=0;
                    first_censored = find(~iscensored,1,'first');
                    last_censored = find(~iscensored,1,'last');
                    if ~numel(first_censored)
                        first_censored=1;
                        last_censored=numel(iscensored);
                    end
                    surf(heatmapI(i).pos_range,hmtmp_1,hmtmp_2,'EdgeColor','none');
                    caxis(ax_);
                    xlabel('AP position');
                    ylabel(['Time (s)']);
                    colormap(hot);
                    hold on;
                    plot3([heatmapI(i).pos_range(first_censored) heatmapI(i).pos_range(first_censored)],[0 1400],[1e3 1e3],'LineWidth',2,'color',[0.5 .5 .5],'LineStyle','--');
                    plot3([heatmapI(i).pos_range(last_censored) heatmapI(i).pos_range(last_censored)],[0 1400],[1e3 1e3],'LineWidth',2,'color',[0.5 .5 .5],'LineStyle','--');
                    xlim([min(heatmapI(i).pos_range) max(heatmapI(i).pos_range)]);
                    ylim([0 1400]);
                    view([0,0,1]);% Adopt horizontal view
                    title(['nc' num2str(i+8)]);
                    hold off;
                end
            end
        end
    end

    function hQuit_CallBack(~,~)
        close(h);
    end

    function hTrim_CallBack(~,~)
        cycle_range=[];
        for i=2:6
            if numel(heatmapI)>=i
                if numel(heatmapI(i).Abs_map)
                    cycle_range=[cycle_range i+8];
                    tlower(1:Nembryo,numel(cycle_range)) = str2double(get(h_fromTime(i-1),'String'));
                    tupper(1:Nembryo,numel(cycle_range)) = str2double(get(h_toTime(i-1),'String'));
                    posborder(1:Nembryo,numel(cycle_range))= str2double(get(h_posBorder(i-1),'String'));
                    if any([2 4 6]==normalize_time)
                        % Scale tlower and tupper by their relative tinterphase
                        for j=1:Nembryo
                            if numel(heatmapI(i).tphase)>=j
                                if i<6 % Only do this with nc<14
                                    tlower(j,numel(cycle_range))=tlower(j,numel(cycle_range))*heatmapI(i).tphase(j)/heatmapI(i).mtphase;
                                    tupper(j,numel(cycle_range))=tupper(j,numel(cycle_range))*heatmapI(i).tphase(j)/heatmapI(i).mtphase;
                                end
                            end
                        end
                    end
                end
            end
        end
        close(h);
    end

    function hDraw_CallBack(~,~)
        cycle_range=[];
        % Draw the intensity curve at the border over time
        for i=2:6
            if numel(heatmapI)>=i
                if numel(heatmapI(i).Abs_map)
                    cycle_range=[cycle_range i+8];
                    posborder(1:Nembryo,numel(cycle_range))= str2double(get(h_posBorder(i-1),'String'));
                    tlower_(1,numel(cycle_range)) = str2double(get(h_fromTime(i-1),'String'));
                    tupper_(1,numel(cycle_range)) = str2double(get(h_toTime(i-1),'String'));
                    % Begin drawing
                    posidx=(heatmapI(i).pos_range-posborder(1,end)-wd/2).*(heatmapI(i).pos_range-posborder(1,end)+wd/2)<0;
                    if normalize_time
                        tmp=mean(heatmapI(i).Rel_map(:,posidx)');
                        tmp_=heatmapI(i).Rel_time;
                    else
                        tmp=mean(heatmapI(i).Abs_map(:,posidx)');
                        tmp_=heatmapI(i).Abs_time;
                    end
                    if numel(cycle_range)==1
                        figure;
                    end
                    plot(tmp_-tlower_(1,end),tmp,'Display',['nc' num2str(cycle_range(end))]); hold on;
                end
            end
        end
        if numel(cycle_range)>0
            legend show;
            xlabel('Time (s)');
            ylabel('PSpot');
        end        
    end
%% Wait to close or automatic trim    
    if auto_trim
        hTrim_CallBack();
    else
        uiwait(h) 
    end
end