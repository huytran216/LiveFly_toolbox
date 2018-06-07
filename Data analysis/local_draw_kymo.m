function [heatmapI,pos_range,draw_time,mtphase,Cellcount]=local_draw_kymo(datamat,ts_spec,xborder,tphase,cycleno,binwidth,AP,time_normalize,time_align)
    % Draw the kymograph with data in datamat
    % Align the embryo by xborder value
    % Draw in nuclear cycle cycleno
    % Binwidth
    % time_normalize: Normalize interphase duration or not
    % time_align: align traces based on 1st spot appearance (left 0 for now)
mkdir('KymoFig');
blend_time=1;
show_movie=0;
barplot=0;
%% Setup up parameters for drawling
maxT=1400;
if time_normalize
    tmin=0;
    tmax=100;
    newdt=2.5;                   % Interval between snapshot (%)
    drawdt=2.5;                  % Drawing intervals (%)
else
    tmin=0;
    tmax=maxT;
    newdt=25;                    % Interval between drawing snapshot (s)
    drawdt=25;                   % Drawing intervals (s)
end
timeunit='s';
color='kgbk';                           % Color for phase in cell cycle
color_embryo='rgbykcrmb';               % Color for embryo
draw_time=[tmin:drawdt:tmax];    % When to show snapshot (for movie=1)
%% Initiating storage variables
snapshot_cellsize=[];
mean_snapshot={};
Irec_={};
pos_range=AP(1):AP(2);
Cellcount=zeros(numel(draw_time),numel(pos_range));
%% Get the data and do interpolation
xaxis_all={};yaxis_all={};
sizecell_all={};idselect_all={};
stay=[];
for tsidx=ts_spec
    % Select the cells of specific movies and cell cycle
    idselect=find(([datamat(:).cycle]==cycleno)&(([datamat(:).tscnt]==tsidx)));
    % Sampling interval
    if numel(idselect)
        dt_all(tsidx)=datamat(idselect(1)).dt;
        stay=[stay tsidx];
        % Take cell mean position and feature value
        xaxis=[datamat(idselect).x]*100;
        yaxis=[datamat(idselect).y]*100;
        yaxis_=arrayfun(@(x) mean(datamat(x).yrec),idselect); % x position, No embryo alignment yet
        scale=median(yaxis_./yaxis);
        sizecell=arrayfun(@(x) mean(datamat(x).sizerec),idselect)/scale;  % Size of cell
        ison=arrayfun(@(x) subindex(datamat(x).Feature,1),idselect); % Is cell valid
        % Record all the features
        xaxis_all{tsidx}=xaxis(ison>=0);
        yaxis_all{tsidx}=yaxis(ison>=0);
        sizecell_all{tsidx}=sizecell(ison>=0);
        idselect_all{tsidx}=idselect(ison>=0);
    end    
end
ts_spec=stay;

%% Process time series separately
% Process each time series separately
for tsidx=ts_spec
    tmp1=idselect_all{tsidx};
    inittime(tsidx)=100000;     % Time of beginning point of the cycle (or 1st spot appearance)
    endtime(tsidx)=0;
    if ~time_align
        % Option 1. Find mitosis of previous cell cycle
        for i=1:numel(tmp1)
            tmp3=datamat(tmp1(i)).time;    % Time record
            if inittime(tsidx)>tmp3(1)
                inittime(tsidx)=tmp3(1);
            end
        end
    else
        % Option 2. Find 1st spot appearance
        for i=1:numel(tmp1)
            tmp3=datamat(tmp1(i)).Intensity;
            tmp4=find(tmp3,1,'first');
            tmp5=datamat(tmp1(i)).time;
            if inittime(tsidx)>tmp5(tmp4)
                inittime(tsidx)=tmp5(tmp4);
            end
        end
    end
    % Extract snapshot intensity from each embryos
    Irec=[];
    Srec=[];
    for i=1:numel(tmp1)
        % Find the time after the 1st spot appearance / mitosis
        alignedidx=find(datamat(tmp1(i)).time>=inittime(tsidx));
        if alignedidx(1)>1
            alignedidx=[alignedidx(1)-1 alignedidx];
        end
        % Load data
        tmp_t=datamat(tmp1(i)).time(alignedidx)-datamat(tmp1(i)).time(alignedidx(1));   % Load Time record
        tmp_I=datamat(tmp1(i)).Intensity(alignedidx);          % Load Intensity record
        tmp_s=datamat(tmp1(i)).sizerec(alignedidx);            % Load spot size record
        % Apply the blending time
        hblend=ones(1,1+round(blend_time/dt_all(tsidx)))/(1+round(blend_time/dt_all(tsidx)));
        tmp_I=conv(tmp_I(end:-1:1),hblend);
        tmp_I=tmp_I(1:end-numel(hblend)+1);
        tmp_I=tmp_I(end:-1:1);
        % Normalize the time
        if time_normalize
            tmp_t=tmp_t/tphase(tsidx)*100;
            if any(tmp_t>=100)
                'alert'
            end
        end
        % Padding limiter at mitosis to
        tmp_I=[-1000000 tmp_I -1000000];
        tmp_s=[tmp_s(1) tmp_s tmp_s(end)];
        tmp_t=[-1000 tmp_t 10000];
        % Extract intensity at snapshot
        Irec(i,:)=interp1(tmp_t,tmp_I,draw_time,'linear');
        Srec(i,:)=interp1(tmp_t,tmp_s,draw_time,'linear');
    end
    Srec((Irec<0)&(Irec>-1))=NaN;    % Non-existing cells
    Srec(Irec<-1)=NaN;               % Non-existing cells

    Irec(Irec>0)=1;                 % Real intensity
    Irec((Irec<0)&(Irec>-1))=-1;    % Non-existing cells
    Irec(Irec<-1)=2;                % Non-existing cells

    Irec_{tsidx}=Irec;              % Record Spot existence
    Srec_{tsidx}=Srec;              % Record Spot size
end
%% Draw the kymograph
    cnt=0;
    draw_stop=0;
    snapshot_interpolated={};
    
    if show_movie
        figure;
    end
    for j=1:numel(draw_time)
        cnt=cnt+1;
        for tsidx=ts_spec
            Irec=Irec_{tsidx};
            Srec=Srec_{tsidx};
            if (j<=size(Irec,2))
                % Plot the mean Pon over time
                pon=[];
                Irectmp=Irec;
                Irectmp(Irec==2)=0;
                Irectmp(Irec<0)=0;
                cnt_=0;
                for pos=pos_range
                    cnt_=cnt_+1;
                    % Find cell with position position
                    tmp=(xaxis_all{tsidx}(:)-50-pos+binwidth/2).*(xaxis_all{tsidx}(:)-50-pos-binwidth/2)<=0;
                    % Calcylate mean Pon
                    pon=[pon mean(Irectmp(tmp,j))];
                    Cellcount(j,cnt_)=Cellcount(j,cnt_)+sum(tmp);
                end
                tmp1=[-50 pos_range-xborder(tsidx) 50];     % Position
                tmp2=[-1000000 pon -1000000];
                snapshot_interpolated{tsidx}=interp1(tmp1,tmp2,pos_range-xborder(tsidx));
                snapshot_interpolated{tsidx}(snapshot_interpolated{tsidx}<0)=NaN;
                if show_movie
                    if ~barplot
                        plot(pos_range,snapshot_interpolated{tsidx},'LineWidth',1,'color',color_embryo(tsidx),'DisplayName',['embryo' num2str(tsidx)]);hold on;
                    end
                end
                % Extract the cell size:
                snapshot_cellsize(tsidx,j)=nanmean(Srec(:,j));
            else
                snapshot_interpolated{tsidx}=pos_range*0;
            end
            % Draw legend and plot beautify
            if (tsidx==ts_spec(end))
                if blend_time==1
                    % Plot mean snap_shot between embryos:
                    tmp=cell2mat(snapshot_interpolated(:));
                    mean_snapshot{j}=nanmean(tmp,1);
                    var_snapshot{j}=nanvar(tmp,[],1);
                    ticktmp=1:1:numel(var_snapshot{j});
                    vartmp=var_snapshot{j}*0;vartmp(ticktmp)=var_snapshot{j}(ticktmp);
                    if show_movie
                        if barplot
                            h=errorbar(pos_range(ticktmp),mean_snapshot{j}(ticktmp),vartmp(ticktmp)/numel(ts_spec),'b');
                            set(h,'LineWidth',1,'LineStyle','-.','DisplayName',['hb-ms2']);hold on;
                        else
                            plot(pos_range,mean_snapshot{j},'LineWidth',2,'color',[0.5 0.5 0.5],'LineStyle','-.','DisplayName',['Mean pattern']);
                        end
                    end
                end
                if show_movie
                    hold on;
                    switch timeunit
                        case '%'
                            title([num2str(draw_time(j)) ' ' timeunit ', nc' num2str(cycleno)]);
                        case 's'
                            title([num2str(draw_time(j)) ' ' timeunit ', nc' num2str(cycleno)]);
                    end
                    xlabel('AP axis (%EL)');
                    ylabel('P_{Spot}');
                    axis([AP 0 1])
                    XTick=get(gca,'XTick');
                    XTickLabel=get(gca,'XTickLabel');
                    hold off;
                end
            end
        end
        if show_movie
            legend('show');
            if ~(draw_stop)
                F(cnt)=getframe(gcf);
                draw_time(j)
            end
        end
    end
    if show_movie
        filename=['KymoFig/mov_nc' num2str(cycleno) '_relativetime' num2str(time_normalize)];
        if exist('filename','var')
            imwrite(F(1).cdata,[filename '.tif'],'writemode', 'overwrite');
            for i=2:numel(F)
                imwrite(F(i).cdata,[filename '.tif'],'writemode', 'append');
            end
        end
    end
%% Plot the kymograph of snapshot vs time
nrowFISH=0;     % Set the width of the FISH stripev (not used)
numdraw=numel(nanmean(snapshot_cellsize,1));
heatmapI=cell2mat(mean_snapshot(:));
heatmapI(isnan(heatmapI))=0;
figure;
mtphase=mean(tphase(ts_spec));
if time_normalize
    draw_time=draw_time*mtphase/100;
end
% Heat map
subplot(222);
surf(pos_range,draw_time,heatmapI,'EdgeColor','none');
hold on;
plot3(pos_range,pos_range*0+mtphase,pos_range*0+1,'LineWidth',2,'LineStyle','--');
colorbar;
caxis([0 1]);
xlabel('AP position');
ylabel(['Time (s)']);
colormap(hot);
xlim(AP);
ylim([0 maxT]);
view([0,0,1]);% Adopt horizontal view
time_axis=get(gca,'ylim');
% Nuclei count:
subplot(224);
surf(pos_range,draw_time,Cellcount,'EdgeColor','none');
hold on;
plot3(pos_range,pos_range*0+mtphase,pos_range*0+max(Cellcount(:)),'LineWidth',2,'LineStyle','--');
colorbar;
%caxis([0 1]);
xlabel('AP position');
ylabel(['Time (s)']);
ylim(time_axis);
colormap(hot);
xlim(AP);
view([0,0,1]);% Adopt horizontal view


% Cell size
subplot(221);
snapshot_cellsize(snapshot_cellsize==0)=NaN;
plot(nanmean(snapshot_cellsize,1),draw_time(1:numdraw));
ylim(time_axis);
xlabel('Cell size (a.u)');
filename=['KymoFig/summary_nc' num2str(cycleno)];
title(['nc' num2str(cycleno)]);
saveas(gcf,[filename '_relativetime' num2str(time_normalize) '.fig']);
