function [heatmapI_output,heatmapI_std,heatmapI_indi,pos_range,draw_time,mtphase,Cellcount]=local_draw_kymo(datamat,ts_spec,xborder,tphase,cycleno,binwidth,AP,time_normalize,time_align,DatasetName,Ymax)
    % Draw the kymograph with data in datamat
    % Align the embryo by xborder value
    % Draw in nuclear cycle cycleno
    % Binwidth
    % time_normalize: Normalize interphase duration or not
    % time_align: align traces based on 1st spot appearance (left 0 for now)
% Input:
    % datamat: data matrix
    % ts_spec: time series selected for analysis
    % xborder: border position per movie
    % cycleno: cycle of interest
    % binwidth: bin size for calculating PSpot
    % AP: AP axis
    % time_normalize: normalize time axis or not
    % time_align: align time or not
% Output:
    % heatmapI_output: mean
        % {1} heatmapI_PSpot: heatmap PSpot
        % {2} heatmapI_PSpot_accum: heatmap accumulating
        % {3} heatmapI_Intensity: heatmap mean intensity
        % {4} heatmapI_PSpot_accum: heatmap mean intensity
    % heatmapI_std: standard deviation
        % {1} heatmapI_PSpot: heatmap PSpot
        % {2} heatmapI_PSpot_accum: heatmap accumulating
        % {3} heatmapI_Intensity: heatmap mean intensity
        % {4} heatmapI_PSpot_accum: heatmap mean intensity
    % heatmapI_indi: individual heatmap
        % {1,i} heatmapI_PSpot: heatmap PSpot of ith movie
        % {2,i} heatmapI_PSpot_accum: heatmap accumulating of ith movie
        % {3,i} heatmapI_Intensity: heatmap mean intensity of ith movie
        % {4,i} heatmapI_PSpot_accum: heatmap mean intensity of ith movie
mkdir('KymoFig');
blend_time=1;   % Blend time
% Check whether to save kymograph or not
if exist('DatasetName','var')
    DatasetName=DatasetName(1:end-4);
else
    DatasetName='';
end
if ~numel(DatasetName)
    show_movie=0;
else
    show_movie=1;
end
barplot=0;      % Show individual movies (0) or merged movie with error bar (1)
output_isIntensity = 1; % Make movie for which isIntensity: (1) Spot, (2) Average intensity
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
draw_time=[tmin:drawdt:tmax];        % When to show snapshot (for movie=1)
%% Initiating storage variables
snapshot_cellsize=[];
mean_snapshot={};
Irec_={};
pos_range=[-45:45];
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
% Process each trace separately
Irec__={};
Srec__={};
for tsidx=ts_spec
    tmp1=idselect_all{tsidx};
    inittime(tsidx)=100000;     % Time of beginning point of the cycle (or 1st spot appearance)
    endtime(tsidx)=0;
    if ~time_align
        % Option 1. Find mitosis of previous cell cycle
        for i=1:numel(tmp1)
            tmp3=datamat(tmp1(i)).Adjustedtime;    % Time record
            if inittime(tsidx)>tmp3(1)
                inittime(tsidx)=tmp3(1);
            end
        end
    else
        % Option 2. Find 1st spot appearance
        for i=1:numel(tmp1)
            tmp3=datamat(tmp1(i)).AdjustedIntensity;
            tmp4=find(tmp3,1,'first');
            tmp5=datamat(tmp1(i)).Adjustedtime;
            if inittime(tsidx)>tmp5(tmp4)
                inittime(tsidx)=tmp5(tmp4);
            end
        end
    end
    % Extract snapshot intensity from each embryos
    Irec=[];        % Holder for intensity - merged embryo
    Srec=[];        % Holder for spot appearance - merged embryo
    for i=1:numel(tmp1)
        % Find the time after the 1st spot appearance / mitosis
        alignedidx=find(datamat(tmp1(i)).Adjustedtime>=inittime(tsidx));
        if alignedidx(1)>1
            alignedidx=[alignedidx(1)-1 alignedidx];
        end
        % Get time and cell size
            % Load data
            tmp_t=datamat(tmp1(i)).Adjustedtime(alignedidx)-datamat(tmp1(i)).Adjustedtime(alignedidx(1));   % Load Time record
            tmp_s=datamat(tmp1(i)).sizerec(alignedidx);            % Load spot size record
            % Normalize the time if needed
            if time_normalize
                tmp_t=tmp_t/tphase(tsidx)*100;
                if any(tmp_t>=100)
                    'alert'
                end
            end
            % Padding signal for mitosis
            tmp_s=[tmp_s(1) tmp_s tmp_s(end)];
            tmp_t=[min(tmp_t)-1000 tmp_t max(tmp_t)+1000];
            Srec(i,:)=interp1(tmp_t,tmp_s,draw_time,'linear');
        % Treating intensity signal - snapshot or cummulative
            for hblendcnt=[1 2]
                % Load data
                tmp_I=datamat(tmp1(i)).AdjustedIntensity(alignedidx);          % Load Intensity record
                % Apply the blending time
                if hblendcnt==2
                    tmp_I=cumsum(tmp_I);
                    if tmp_t(end-1)<draw_time(end)
                        tmp_t(end-1)=draw_time(end);
                    end
                end
                % Padding limiter at mitosis to
                tmp_I=[-1000000 tmp_I -1000000];
                % Extract intensity at draw_time
                Irec{hblendcnt}(i,:)=interp1(tmp_t,tmp_I,draw_time,'linear');                
            end
    end
    for hblendcnt=[1 2]             % Blendtime: (1)snapshot or (2)cummulative
        for isIntensity=[1 2]       % isIntensity: (1)Spot or (2)Intensity
            Irec_=Irec{hblendcnt};
            if isIntensity==1
                Irec_(Irec_>0)=1;                   % Real intensity
                Irec_((Irec_<0)&(Irec_>-1))=-1;     % Non-existing cells
                Irec_(Irec_<=-1)=-1;                % Non-existing cells
            end
            Irec__{tsidx,hblendcnt,isIntensity}=Irec_;              % Record Spot existence/intensity
            Srec((Irec{hblendcnt}<0)&(Irec{hblendcnt}>-1))=NaN;     % Non-existing cells
            Srec(Irec{hblendcnt}<-1)=NaN;                           % Non-existing cells
            Srec_{tsidx}=Srec;                                      % Record Spot size
        end
    end
end
%% Make the kymograph - interpolating with common grid of position x time
    drawcnt=0;              % Drawing frame number
    draw_stop=0;
    snapshot_indi={};   % Holder for individual heatmap (ts_spec x time x hblend x isIntensity)
    mtphase=mean(tphase(ts_spec));
    
    
    if show_movie
        h=figure;
    end
    
    for hblendcnt=[1 2]
        for isIntensity=[1 2]
            for j=1:numel(draw_time)
                snapshot_interpolated={};   % Sub holder for individual heatmap
                for tsidx=ts_spec
                    Irec=Irec__{tsidx,hblendcnt,isIntensity};
                    Srec=Srec_{tsidx};
                    if (j<=size(Irec,2))
                        % Plot the mean Pon over time
                        pon=[];
                        Irectmp=Irec;
                        Irectmp(Irectmp<0)=0;
                        cnt_=0;
                        for pos=pos_range
                            cnt_=cnt_+1;
                            % Find cell with position position
                            tmp=(xaxis_all{tsidx}(:)-50-pos+binwidth/2).*(xaxis_all{tsidx}(:)-50-pos-binwidth/2)<=0;
                            % Calcylate mean Pon or intensity
                            pon=[pon nanmean(Irectmp(tmp,j))];
                            if (hblendcnt==1)&&(isIntensity==1)
                                Cellcount(j,cnt_)=Cellcount(j,cnt_)+sum(tmp);
                            end
                        end
                        tmp1=[-50 pos_range-xborder(tsidx) 50];     % Position
                        tmp2=[-1000000 pon -1000000];
                        snapshot_interpolated{tsidx}=interp1(tmp1,tmp2,pos_range-xborder(tsidx));
                        snapshot_interpolated{tsidx}(snapshot_interpolated{tsidx}<0)=NaN;
                        if show_movie&(hblendcnt==1)&(isIntensity==output_isIntensity)
                            if ~barplot % Show individual movies
                                plot(pos_range,snapshot_interpolated{tsidx},'LineWidth',1.5,'color',corder(tsidx),'DisplayName',[num2str(tsidx)]);hold on;
                            end
                        end
                        % Extract the cell size:
                        snapshot_cellsize(tsidx,j)=nanmean(Srec(:,j));
                    else
                        snapshot_interpolated{tsidx}=pos_range*0;
                    end
                    % Record individual heatmap:
                    snapshot_indi{tsidx,j,hblendcnt,isIntensity} = snapshot_interpolated{tsidx};
                    % Draw legend and plot beautify
                    if (tsidx==ts_spec(end))
                        if blend_time==1
                            % Plot mean snap_shot between embryos:
                            tmp=cell2mat(snapshot_interpolated(:));
                            mean_snapshot{j,hblendcnt,isIntensity}=nanmean(tmp,1);
                            var_snapshot{j,hblendcnt,isIntensity}=nanvar(tmp,[],1);
                            ticktmp=1:1:numel(var_snapshot{j,hblendcnt,isIntensity});
                            vartmp=var_snapshot{j,hblendcnt,isIntensity}*0;
                            vartmp(ticktmp)=var_snapshot{j,hblendcnt,isIntensity}(ticktmp);
                            if show_movie&(hblendcnt==1)&(isIntensity==output_isIntensity)
                                if barplot % Show merged movie with error bar
                                    h=errorbar(pos_range(ticktmp),mean_snapshot{j,hblendcnt,isIntensity}(ticktmp),vartmp(ticktmp)/numel(ts_spec),'b');
                                    set(h,'LineWidth',1,'LineStyle','-.','DisplayName',['Merged']);hold on;
                                else
                                    plot(pos_range,mean_snapshot{j,hblendcnt,isIntensity},'LineWidth',2,'color',[0.5 0.5 0.5],'LineStyle','-.','DisplayName',['Mean']);
                                end
                            end
                        end
                        % CHANGE isIntensity=1 if SPot, else 2 for Average intensity
                        if show_movie&(hblendcnt==1)&(isIntensity==output_isIntensity)
                            hold on;
                            switch timeunit
                                case '%'
                                    title([num2str(draw_time(j)*mtphase/100) ' ' timeunit ', nc' num2str(cycleno)]);
                                case 's'
                                    title([num2str(draw_time(j)*mtphase/100) ' ' timeunit ', nc' num2str(cycleno)]);
                            end
                            xlabel('AP axis (%EL)');
                            if output_isIntensity==2
                                ylabel('Spot intensity');
                            else
                                ylabel('P_{Spot}');
                                Ymax=1.1;
                            end
                            axis([AP 0 Ymax])
                            XTick=get(gca,'XTick');
                            XTickLabel=get(gca,'XTickLabel');
                            hold off;
                        end
                    end
                end
                if show_movie&(hblendcnt==1)&(isIntensity==output_isIntensity)
                    drawcnt=drawcnt+1;
                    legend('show');
                    if ~(draw_stop)
                        F(drawcnt)=getframe(gcf);
                        draw_time(j)
                    end
                end
            end
        end
    end
    
    if show_movie
        filename=['KymoFig/' DatasetName '_nc' num2str(cycleno) '_isIntensity' num2str(output_isIntensity) '.tif'];
        if exist(filename,'file')
            imwrite(F(1).cdata,filename,'writemode', 'overwrite');
            for i=2:numel(F)
                imwrite(F(i).cdata,filename,'writemode', 'append');
            end
        else
            for i=1:numel(F)
                imwrite(F(i).cdata,filename,'writemode', 'append');
            end
        end
        close(h);
    end
%% Plot the kymograph of snapshot vs time
nrowFISH=0;     % Set the width of the FISH stripev (not used)
numdraw=numel(nanmean(snapshot_cellsize,1));
heatmapI_output={};
heatmapI_std={};

if time_normalize
            draw_time=draw_time*mtphase/100;
end

for hblendcnt=[1 2]
    for isIntensity=[1 2]
        % mean heatmap
            heatmapI=cell2mat(mean_snapshot(:,hblendcnt,isIntensity));
            heatmapI(isnan(heatmapI))=0;
        % std heatmap
            heatmapI_=sqrt(cell2mat(var_snapshot(:,hblendcnt,isIntensity)));
            heatmapI_(isnan(heatmapI_))=0;
        % Indi heatmap:
            for tsidx = ts_spec
                tmp = cell2mat(snapshot_indi(tsidx,:,hblendcnt,isIntensity)');
                %tmp(isnan(tmp))=0;
                heatmapI_indi{tsidx,hblendcnt,isIntensity}=tmp;
            end
        % updating output
            heatmapI_output{hblendcnt,isIntensity}=heatmapI;
            heatmapI_std{hblendcnt,isIntensity}=heatmapI_;            
        % Plot kymograph if needed        
%             if show_movie
%                 figure;
%                 % Heat map
%                 subplot(222);
%                 surf(pos_range,draw_time,heatmapI,'EdgeColor','none');
%                 hold on;
%                 plot3(pos_range,pos_range*0+mtphase,pos_range*0+1,'LineWidth',2,'LineStyle','--');
%                 colorbar;
%                 if isIntensity==1
%                     caxis([0 1]);
%                 end
%                 xlabel('AP position');
%                 ylabel(['Time (s)']);
%                 colormap(hot);
%                 xlim(AP);
%                 ylim([0 maxT]);
%                 view([0,0,1]);% Adopt horizontal view
%                 time_axis=get(gca,'ylim');
%                 % Nuclei count:
%                 subplot(224);
%                 surf(pos_range,draw_time,Cellcount,'EdgeColor','none');
%                 hold on;
%                 plot3(pos_range,pos_range*0+mtphase,pos_range*0+max(Cellcount(:)),'LineWidth',2,'LineStyle','--');
%                 colorbar;
%                 %caxis([0 1]);
%                 xlabel('AP position');
%                 ylabel(['Time (s)']);
%                 ylim(time_axis);
%                 colormap(hot);
%                 xlim(AP);
%                 view([0,0,1]);% Adopt horizontal view
% 
% 
%                 % Cell size
%                 subplot(221);
%                 snapshot_cellsize(snapshot_cellsize==0)=NaN;
%                 plot(nanmean(snapshot_cellsize,1),draw_time(1:numdraw));
%                 ylim(time_axis);
%                 xlabel('Cell size (a.u)');
%                 filename=['KymoFig/' DatasetName 'sum_nc' num2str(cycleno)];
%                 title(['nc' num2str(cycleno)]);
%                 saveas(gcf,[filename '_relativetime' num2str(time_normalize) '.fig']);
%             end
    end
end
