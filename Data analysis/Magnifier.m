function [tlower,tupper,cycle_range]=Magnifier(h,heatmapI)
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

% Interface - plot absolute or relative time
    hmessages_Time = uicontrol('Style','text','String','Time axis:',...
        'Position',[30,10,100,20],'FontSize',10, 'HorizontalAlignment','left',...
        'ForegroundColor','white','BackgroundColor',[0.7 0.7 0.7]);
    hTime = uicontrol('Style','popup','String',{'Absolute','Relative'},...
        'Position',[100,10,120,20],'FontSize',10, 'HorizontalAlignment','left',...
        'Callback',@hset_Time);

% Plot axis: 1-5 <=> 10-14
Nembryo = 0;
for axidx=1:5
    axmap(axidx) = axes('Position',[(axidx-1)*0.2+0.03 0.30 0.15 0.65],'XTickLabel',[],'YTickLabel',[]);
    hmessage_fromTime(axidx)=uicontrol('Style','text','String','From (s):',...
        'Position',[(axidx-1)*240+36,100,100,20],'FontSize',10, 'HorizontalAlignment','left',...
        'ForegroundColor','white','BackgroundColor',[0.7 0.7 0.7]);
    h_fromTime(axidx) = uicontrol('Style','edit','String',num2str(0),...
        'Position',[(axidx-1)*240+100,100,80,20],'FontSize',10, 'HorizontalAlignment','left');
    hmessage_toTime(axidx)=uicontrol('Style','text','String','To (s):',...
        'Position',[(axidx-1)*240+36,70,100,20],'FontSize',10, 'HorizontalAlignment','left',...
        'ForegroundColor','white','BackgroundColor',[0.7 0.7 0.7]);
    tmp=10000;
    if numel(heatmapI)>=axidx+1
        if numel(heatmapI(axidx+1).mtphase)
            tmp=num2str(heatmapI(axidx+1).mtphase);
            if Nembryo < numel(heatmapI(axidx+1).tphase)
                Nembryo = numel(heatmapI(axidx+1).tphase);
            end
        end
    end
    h_toTime(axidx) = uicontrol('Style','edit','String',tmp,...
        'Position',[(axidx-1)*240+100,70,80,20],'FontSize',10, 'HorizontalAlignment','left');
end

% Trim button:
hQuit = uicontrol('Style','pushbutton','String', 'Quit',...
        'Position',[1000,20,90,40],...
        'Callback',@hQuit_CallBack);
hTrim = uicontrol('Style','pushbutton','String', 'Trim interphase',...
        'Position',[880,20,100,40],...
        'Callback',@hTrim_CallBack);
% Plot the plot
normalize_time=0;
Draw_kymo();

% Default outputs
tlower = [];
tupper = [];
cycle_range=[];
%% Functions:
    % Set time axis:
    function hset_Time(~,~)
        normalize_time=get(hTime,'Value')-1;
        Draw_kymo;
    end
    % Plot heatmap
    function Draw_kymo()
        for i=2:6
            if numel(heatmapI)>=i
                if numel(heatmapI(i).Abs_map)
                    axes(axmap(i-1));
                    plot3(heatmapI(i).pos_range,heatmapI(i).pos_range*0+heatmapI(i).mtphase,heatmapI(i).pos_range*0+1,'LineWidth',2,'LineStyle','--');
                    hold on;
                    if normalize_time
                        surf(heatmapI(i).pos_range,heatmapI(i).Rel_time,heatmapI(i).Rel_map,'EdgeColor','none');
                    else
                        surf(heatmapI(i).pos_range,heatmapI(i).Abs_time,heatmapI(i).Abs_map,'EdgeColor','none');
                    end
                    caxis([0 1]);
                    xlabel('AP position');
                    ylabel(['Time (s)']);
                    colormap(hot);
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
        for i=2:6
            if numel(heatmapI)>=i
                if numel(heatmapI(i).Abs_map)
                    cycle_range=[cycle_range i+8];
                    tlower(1:Nembryo,numel(cycle_range)) = str2double(get(h_fromTime(i-1),'String'));
                    tupper(1:Nembryo,numel(cycle_range)) = str2double(get(h_toTime(i-1),'String'));
                    if normalize_time
                        % Scale tlower and tupper by their relative tinterphase
                        for j=1:Nembryo
                            if numel(heatmapI(i).tphase)>=j
                                tlower(j,numel(cycle_range))=tlower(j,numel(cycle_range))*heatmapI(i).tphase(j)/heatmapI(i).mtphase;
                                tupper(j,numel(cycle_range))=tupper(j,numel(cycle_range))*heatmapI(i).tphase(j)/heatmapI(i).mtphase;
                            end
                        end
                    end
                end
            end
        end
        close(h);
    end
%% Wait to close
    uiwait(h) 
end