function [tlower,tupper,cycle_range,posborder]=Magnifier(h,heatmapI,wd)
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
    % The axes for plotting kymograph
    axmap(axidx) = axes('Position',[(axidx-1)*0.2+0.03 0.40 0.15 0.55],'XTickLabel',[],'YTickLabel',[]);
    % Set beginning time
    hmessage_fromTime(axidx)=uicontrol('Style','text','String','From (s):',...
        'Position',[(axidx-1)*240+36,170,100,20],'FontSize',10, 'HorizontalAlignment','left',...
        'ForegroundColor','white','BackgroundColor',[0.7 0.7 0.7]);
    h_fromTime(axidx) = uicontrol('Style','edit','String',num2str(0),...
        'Position',[(axidx-1)*240+100,170,80,20],'FontSize',10, 'HorizontalAlignment','left');
    % Set ending time
    hmessage_toTime(axidx)=uicontrol('Style','text','String','To (s):',...
        'Position',[(axidx-1)*240+36,140,100,20],'FontSize',10, 'HorizontalAlignment','left',...
        'ForegroundColor','white','BackgroundColor',[0.7 0.7 0.7]);
    tmp1=10000;
    tmp2=0;
    if numel(heatmapI)>=axidx+1
        if numel(heatmapI(axidx+1).mtphase)
            tmp1=num2str(heatmapI(axidx+1).mtphase);
            if isfield(heatmapI,'posborder')
                tmp2=num2str(heatmapI(axidx+1).posborder);
            end
            if Nembryo < numel(heatmapI(axidx+1).tphase)
                Nembryo = numel(heatmapI(axidx+1).tphase);
            end
        end
    end
    h_toTime(axidx) = uicontrol('Style','edit','String',num2str(tmp1),...
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
        cycle_range=[];
        for i=2:6
            if numel(heatmapI)>=i
                if numel(heatmapI(i).Abs_map)
                    cycle_range=[cycle_range i+8];
                    tlower(1:Nembryo,numel(cycle_range)) = str2double(get(h_fromTime(i-1),'String'));
                    tupper(1:Nembryo,numel(cycle_range)) = str2double(get(h_toTime(i-1),'String'));
                    posborder(1:Nembryo,numel(cycle_range))= str2double(get(h_posBorder(i-1),'String'));
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
%% Wait to close
    uiwait(h) 
end