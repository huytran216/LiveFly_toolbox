%% Re-extract all features in all strains before running the script
% THe data should be in Data_analysis/
function LiveFly_COMPARE_1x2x(compare_list,isBcd1X,fea_range,trimmed_trace,nc_range)
warning off;
close all;
%% Params
fld='../../Data analysis/';
load([fld 'feature_label.mat']);

Nsample_min=5;              % Minimum total nuclei per position
Nsample_per_embryo_min=5;   % Minimum nuclei per embryo per position
Nembryo_min = 0;            % Mininum number of embryo per position

plot_embryo_error=1;    % plot error based on embryo diversity (1) or nuclei error (merged, 0)
shaded_error_bar = 1;   % Plot shaded errorbar or normal errorbar


smooth_curve = 1;
normalize_intensity = 0;

fit_lambda = 1;
    impose_fit = 1;
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

if nargin==0
    %compare_list = [1 2 5 3 7];isBcd1X = [0 0 0 0 0]; % For B6-B9-B12 comparison
    %compare_list = [1 2 10]; isBcd1X = [0 0 0]; % For hb-B6-H6B6 comparison
    compare_list = [1 1];isBcd1X=[0 1]; % For hb-B6-H6B6 comparison, 1x2x
    %compare_list = [12 12];isBcd1X=[zeros(1,numel(compare_list)/2) ones(1,numel(compare_list)/2)]; % For hb-B6-H6B6 comparison, 1x2x
    %compare_list = [1 8 9]; isBcd1X =[0 0 0 ];% For vk33 vs random insertion
    %compare_list = [7 7];isBcd1X=[0 1];
    %compare_list = [1 2 10 12]; isBcd1X = compare_list*0;
    %compare_list = [3 3];isBcd1X=[0 1];
end

% Set folder containing mean data (contain dash)
folder={};
folder{1}='tmp/';
folder{2}='tmp_trimmed/';

avr = [600 750 1100];                % Mean nc13 duration
%% Cook label_list
DatasetLabel = {dtset(compare_list).label};
DatasetFile = {dtset(compare_list).filename};
for i=1:numel(compare_list)
    switch isBcd1X(i)
        case 0
            DatasetLabel{i}=[DatasetLabel{i} ''];
        case 1
            DatasetFile{i}=[DatasetFile{i} '-Bcd1X'];
            compare_1x2x = true;
        case 2
            DatasetFile{i}=[DatasetFile{i} '-dBcd'];
            compare_1x2x = true;
    end
end
%% Feature to plot, plot settings
if nargin==0
    fea_range=[1];
    trimmed_trace = 0;      % trimmed trace (1) or not (0)
    nc_range=[13];
end

AP_limit = [-32 20]; % for B6-B9-B12

scanwindow = [-25 20];
%% Plot stuffs
ymax=zeros(1,16);
ymin=zeros(1,16)+1e5;
rec_fea=[];
mI_rec = {};
sI_rec = {};
pos_rec = {};

h=[];
hplot = [];
xborder = []; hborder = []; vborder = []; noembryo = [];
cnt=0;
cnt2=0;

for nc=nc_range
    for fea=fea_range
        cnt=cnt+1;
        for i=1:numel(compare_list)
            % Check if 1x or 2x:
            if isBcd1X(i)
                original_i = i - numel(compare_list)/2;
            else
                original_i = i;
            end
            cnt2=cnt2+1;
            load(fullfile(fld,folder{trimmed_trace+1},DatasetFile{i}),'pos_range','mf_rec','sf_rec','nf_rec','nf_indi','mf_indi','sf_indi','FitRes');
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
                        ne_rec=ne_rec+(nf_indi{1,nc,j}>Nsample_per_embryo_min);
                    end
                end
            % Calculate mean of indi curve:
                for j=1:size(nf_indi,3)
                    if numel(nf_indi{1,nc,j})
                        tmp=mf_indi{fea,nc,j};
                        tmp(isnan(tmp))=0;
                        nf_indi{1,nc,j}(isnan(nf_indi{1,nc,j}))=0;
                        mf_indi_=mf_indi_ + tmp.*(nf_indi{1,nc,j}>Nsample_per_embryo_min);
                        sf_indi_=sf_indi_ + (tmp.*(nf_indi{1,nc,j}>Nsample_per_embryo_min)).^2;
                    end
                end
                mf_indi_=mf_indi_./ne_rec;
                sf_indi_=sqrt(sf_indi_./ne_rec - mf_indi_.^2);
            % Select if certain bin has enough samples
                flttmp = (nf_rec{fea,nc}>=Nsample_min)&(ne_rec>=Nembryo_min);
                
                if plot_embryo_error
                    mtmp = mf_indi_(flttmp);
                    stmp = sf_indi_(flttmp)./sqrt(ne_rec(flttmp)-1);
                else
                    mtmp = mf_rec{fea,nc}(flttmp);
                    stmp = sf_rec{fea,nc}(flttmp)./sqrt(nf_rec{fea,nc}(flttmp)-1);
                end
            % Smooth curve if specified
                if (~fit_lambda) & smooth_curve
                    mtmp = smooth(mtmp(end:-1:1));
                    mtmp = mtmp(end:-1:1);
                    %stmp = smooth(stmp(end:-1:1));
                    %stmp = stmp(end:-1:1);
                    %stmp(stmp_tmp)=nan;
                end
                stmp_tmp = isnan(stmp);
                mtmp(stmp_tmp)=nan;
                mtmp = mtmp(:);
            % Normalize by intensity
                if normalize_intensity
                    mtmp = mtmp/vborder(i,fea,nc,1);
                    stmp = stmp/vborder(i,fea,nc,1);
                end
            % Determine color based on 1x or 2x
                if isBcd1X(i)==0
                    color=corder(4);
                else
                    color=corder(2);
                end
            % Plot mean curve with error, in individual figures
                figure(20+original_i);
                set(gcf,'Position',[680   438   584   240]);
                title(DatasetLabel{original_i});
                if (numel(nc_range)==1)&&(numel(fea_range)==1)
                    % Calculating boundary shift
                    subplot(1,2,1);
                else
                    % Plot only one construct
                    subplot(numel(nc_range),numel(fea_range),cnt);
                end
                if shaded_error_bar
                    h20(cnt2)=shadedErrorBar(pos_range(flttmp),mtmp,stmp,{'color',color,'Display',DatasetLabel{i},'LineWidth',2},1,2);
                else
                    h20(cnt2)=errorbar(pos_range(flttmp),mtmp,stmp,'Display',DatasetLabel{i},'color',color);
                end
                xlabel('AP axis (%EL)');
                ylabel(feature_label{fea});
                xlim(AP_limit);
                hold on;
            % Plot mean curve with error, all in one figure
                figure(40);
                set(gcf,'Position',[680   438   numel(compare_list)*284/2   240]);
                subplot(1,numel(compare_list)/2,original_i);
                title(DatasetLabel{original_i});
                if shaded_error_bar                    
                    h40(cnt2)=shadedErrorBar(pos_range(flttmp),mtmp,stmp,{'color',color,'Display',DatasetLabel{i},'LineWidth',2},1,2);
                else
                    h40(cnt2)=errorbar(pos_range(flttmp),mtmp,stmp,'Display',DatasetLabel{i},'color',color,'LineWidth',2);
                end
                xlabel('AP axis (%EL)');
                %if original_i==1
                    ylabel(feature_label{fea});
                    ylim([0 1]);
                %else
                    %set(gca,'YTick',[]);
                %end
                xlim(AP_limit);
                hold on;
            % Plot prediction of position based on expression alone
                figure(80+original_i);
                subplot(1,2,1+(isBcd1X(i)>0));
                position_prediction_map(pos_range(flttmp),mtmp,stmp);
                title(DatasetLabel{i});
                
            % Saving for comparison between 1x 2x, unobservable point set
            % to either or plateau
                % Fill all NaN
                
                mIall = pos_range*0;mIall(flttmp)=fillmissing(mtmp,'nearest');
                sIall = pos_range*0;sIall(flttmp)=fillmissing(stmp,'nearest');
                tmpf = find(flttmp.*~isnan(mIall),1,'first');
                mIall(1:tmpf-1)=nanmean(mtmp(pos_range(flttmp)<scanwindow(1)));
                sIall(1:tmpf-1)=nanmean(stmp(pos_range(flttmp)<scanwindow(1)));
                
                mI_rec{original_i,1+(isBcd1X(i)>0)}=mIall;
                sI_rec{original_i,1+(isBcd1X(i)>0)}=sIall;
                pos_rec{original_i,1+(isBcd1X(i)>0)}=pos_range;
        end
        % Make a table out of it: constructs are put vertically
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
                    case 16
                        header{cnti,6}=[num2str(vborder(i,fea,nc,1)/100,'%.1f') ' ' char(177) ' ' ...
                            num2str((vborder(i,fea,nc,3) - vborder(i,fea,nc,2))/100,'%.1f')];
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
            htmp=figure(80);set(htmp,'Name',['nc ' num2str(nc),', Feature: ' feature_label{fea}],'Position',[100 100 800 500]);
            htmptable =   uitable('Parent',htmp,'Position',[0 0 800 500]);
            htmptable.Data = header;
            set(htmptable,'ColumnEditable',true(1,10))            
    end
end
%% Plot border shift vs position
fea = fea_range(1);
for nc=nc_range
    for i=1:numel(compare_list)/2
        xtmp(i,nc) = xborder(i,fea,nc,1);
        sxtmp(i,nc) = (xborder(i,fea,nc,3) - xborder(i,fea,nc,2))/2;
        dxtmp(i,nc) = xborder(i+numel(compare_list)/2,fea,nc,1) - xborder(i,fea,nc,1);
        sdxtmp(i,nc) = (xborder(i+numel(compare_list)/2,fea,nc,3) - xborder(i+numel(compare_list)/2,fea,nc,2))/2 + (xborder(i,fea,nc,3) - xborder(i,fea,nc,2))/2;
    end
    figure;
    errorbar(xtmp(:,nc),dxtmp(:,nc),sdxtmp(:,nc),sdxtmp(:,nc),sxtmp(:,nc),sxtmp(:,nc),'.');
    xlabel('Original position  (%EL)');
    ylabel('Mean boundary shift (%EL)');
end
%% Estimation of displace:
% Predict new position:
diff_plot={};
diff_grid={};
for i=1:numel(compare_list)/2
    % Ignore points w too little expression:
%     sI_rec{i,1}(mI_rec{i,1}<max(mI_rec{i,1})/20)=0;
%     mI_rec{i,1}(mI_rec{i,1}<max(mI_rec{i,1})/20)=0;
%     sI_rec{i,2}(mI_rec{i,2}<max(mI_rec{i,2})/20)=0;
%     mI_rec{i,2}(mI_rec{i,2}<max(mI_rec{i,2})/20)=0;
    pos_range1 = pos_rec{i,1}>=-30;
    pos_range2 = pos_rec{i,2}>=-30;
    figure(20+i);
    subplot(1,2,2);
    figure;
    [diff_plot{i},ax,ay,diff_grid{i},prange1,prange2]=shift_prediction_map(pos_rec{i,1}(pos_range1),mI_rec{i,1}(pos_range1),sI_rec{i,1}(pos_range1),...
        pos_rec{i,2}(pos_range2),mI_rec{i,2}(pos_range2),sI_rec{i,2}(pos_range2));
    
    %subplot(121);
    %xlim([-30 20]);
    %ylim([-30 20]);
    
    %subplot(122);
    
    set(gcf,'Position',[680   738   284   240]);
    xlim([-30 20]);
    ylim([-30 20]);
    caxis([0 0.5]);    
    title(DatasetLabel{i});
    colormap(flipud(gray));
    if i==1
        diff_all = diff_plot{1};
        diff_grid_all = diff_grid{1};
    else
        diff_all = diff_all*diff_plot{i};
        diff_grid_all = diff_grid_all.*diff_grid{i};
    end
%     figure(6);
%     subplot(1,numel(compare_list)/2,i);
%     mDX = [];
%     sDX = [];
%     for j=1:sum(pos_range1)
%         %mDX(j) = sum(diff_plot{i}(j,:)'.*(ay(:,j)-ax(:,j)));
%         %sDX(j) = sqrt(sum(diff_plot{i}(j,:)'.*((ay(:,j)-ax(:,j)).^2)) - mDX(j).^2);
%         mDX(j) = sum(diff_grid{i}(j,:).*pos_range(pos_range1));
%         sDX(j) = sqrt(sum(diff_grid{i}(j,:).*pos_range(pos_range1).^2) - mDX(j).^2);
%     end
%     errorbar(pos_range(pos_range1),mDX,sDX);
%     xlim([-30 20]);
%     title(DatasetLabel{i});
end

    % Plot combined probability map
        figure(7);
        HeatMap_(log(diff_grid_all'+1e-10),prange1,prange2,[log(1e-4) log(max(diff_grid_all(:)))]);
        plot3(prange1(:),prange1(:)*0,prange1(:)*0+log(max(diff_grid_all(:))),'LineStyle','-','color','k','LineWidth',1);
        hold on;
        %h = plot3(prange1(:),prange1(:)*0 - 14,prange1(:)*0+log(max(diff_grid_all(:))),'LineStyle','--','color','r','LineWidth',1,'Display','Theory');
        %legend(h,'Theory');
        
        set(gca,'Ydir','normal');
        colormap(flipud(gray));
        %xlabel('Original nuclei position (%EL)');
        %ylabel('Predicted shift (%EL)');
        xlim([-30 20]);
        ylim([-30 20]);
        zlim([log(1e-20) -1])
        set(gcf,'Position',[680   675   414   303]);
%% Summary
% figure(107);
% subplot(121);
% HeatMap_(diff_all',ax,ay-ax,[0 max(diff_all(:))]);
% set(gca,'YDir','normal');
% xlabel('Original position x (%EL)');
%     ylabel('Predicted shift x''-x (%EL)');
%     xlim([-30 10]);
%     ylim([-30 10]);
%     colormap(flipud(gray));
% 
% subplot(122);
% imshow(dualcolormap(diff_grid{1}',diff_grid{2}',[0 0.5],[0 0.5]),'XData',pos_range,'YData',pos_range);
% set(gca,'YDir','normal');
% xlabel('Original position x (%EL)');
%     ylabel('Predicted shift x''-x (%EL)');
%     xlim([-30 10]);
%     ylim([-30 10]);
    
% subplot(122);
%     mDX = [];
%     sDX = [];
%     for j=1:numel(pos_range)
%         % Cut the tail in the distribution: take 95% of mass distribution
%         [maxpro,midpos] = max(diff_all(j,:));
%         diff_all(j,diff_all(j,:)<maxpro/10)=0;
%         diff_all(j,:)=diff_all(j,:)/sum(diff_all(j,:));
%         
%         mDX(j) = sum(diff_all(j,:)'.*(ay(:,j)-ax(:,j)));
%         sDX(j) = sqrt(sum(diff_all(j,:)'.*((ay(:,j)-ax(:,j)).^2)) - mDX(j).^2);
%     end
%     errorbar(pos_range,mDX,sDX);
%     xlabel('Original position x (%EL)');
%     ylabel('Predicted shift x''-x (%EL)');
%     xlim([-30 20]);
%     ylim([-20 0]);
    %set(gcf,'Position',tmpfig);
%% Fit displacement with a curve:

        syms x;
        model_range = {'exp'};
        dat = struct;
        for modelidx = 1:numel(model_range)
            model = model_range{modelidx};
            switch model 
                case 'exp'
                    y = sym('y'); fun = exp(-x/y);  % Exponential gradient
                case 'hybrid'
                    y = sym('y',[1 3]); fun = exp(-x/y(1)) + y(2)*exp(-x/y(3));  % Two Exponential gradient
                case 'power'
                    y = sym('y', [1 2]);fun  = (x+y(1))^(-y(2)); % Power gradient
            end
            dat(modelidx).x = x;
            dat(modelidx).y = y;
            dat(modelidx).fun = fun;

            ax_ = prange1(1,:);
            ax_fit = [-30:10]; % Range of scan for shift in 2X, that is used in fitting:
            ratio = 0.5;      % change in Bcd level

            infsmall = 1e-10;
    %         newfun = @(y0) -sum(log(infsmall+diff_all(sub2ind(size(diff_all),...
    %             -ax_(1)+1+ax_fit,...
    %             max([ones(1,numel(ax_fit)); ...
    %             min([numel(ax_)*ones(1,numel(ax_fit)); ...
    %             (-ax_(1)+1+ax_fit+round(find_displacement(y0,fun,y,x,ax_fit,ratio)))
    %             ])...
    %             ])...
    %             ))));
    %          dat(modelidx).newfun = @(y0) -sum(log(1e-10+diff_all(sub2ind(size(diff_all),...
    %             -ax_(1)+1+ax_fit,...
    %             max([ones(1,numel(ax_fit)); ...
    %             min([numel(ax_)*ones(1,numel(ax_fit)); ...
    %             (-ax_(1)+1+ax_fit+round(find_displacement(y0,dat(modelidx).fun,dat(modelidx).y,dat(modelidx).x,ax_fit,ratio)))
    %             ])...
    %             ])...
    %             ))));

            newfun = @(y0) -sum(log(infsmall+diff_grid_all(sub2ind(size(diff_grid_all),...
                -ax_(1)+1+ax_fit,...    % Original position
                -ax_(1)+1+round(find_displacement(y0,fun,y,x,ax_fit,ratio))... % Shift position
                ))));

            dat(modelidx).newfun = @(y0) -sum(log(1e-10+diff_all(sub2ind(size(diff_all),...
                -ax_(1)+1+ax_fit,...            
                (-ax_(1)+1+round(find_displacement(y0,dat(modelidx).fun,dat(modelidx).y,dat(modelidx).x,ax_fit,ratio)))...
                ))));
            
            nitr = 10;
            if fit_lambda
                switch model 
                    case 'exp'
                        beta_best = cell(1,nitr);
                        f0 = zeros(1,nitr);
                        parfor cnt = 1:nitr
                            x0 = randbetween(1,30);
                            [beta_best_,f0_]=fminsearchbnd(newfun,x0,1,30);
                            beta_best{cnt} = beta_best_;
                            f0(cnt) = f0_;
                            display([cnt f0(cnt)]);
                        end
                        [f0,cnt]=min(f0);
                        beta_best = beta_best{cnt};
                    case 'hybrid'
                        beta_best = cell(1,nitr);
                        f0 = zeros(1,nitr);
                        parfor cnt = 1:nitr
                            x0 = randbetween([1 0 1],[30 10 30]);
                            [beta_best_,f0_]=fminsearchbnd(newfun,x0,[1 0 1],[30 10 30]);
                            beta_best{cnt} = beta_best_;
                            f0(cnt) = f0_;
                            display([cnt f0(cnt)]);
                        end
                        [f0,cnt]=min(f0);
                        beta_best = beta_best{cnt};
                    case 'power'
                        [beta_best,f0]=fminsearchbnd(newfun,[35 2],[-ax_fit(1)+1 0.1],[1000 100]);
                end

                dat(modelidx).beta_best = beta_best;
                dat(modelidx).f0 = f0;
            else
                % imposing lambda
            end
        end
    %% Plot the results

    for i=1:numel(compare_list)/2
        figure(120+i);
        plot(pos_rec{i,1},mI_rec{i,1},'-b','Display','2X');hold on;
        plot(pos_rec{i,2},mI_rec{i,2},'-r','Display','1X');
        title(DatasetLabel{i});
        xlim([AP_limit]);
    end

    figure(107);
        HeatMap_(log(diff_grid_all'+1e-10),prange1,prange2,[log(1e-10) log(max(diff_grid_all(:)))]);
        hold on;
        h1 = plot3(prange1(:),prange1(:)*0 - 14,prange1(:)*0+log(max(diff_grid_all(:))),'LineStyle','--','color','r','LineWidth',2);
        
        set(gca,'Ydir','normal');
        colormap(flipud(gray));
        %xlabel('Original nuclei position (%EL)');
        %ylabel('Predicted shift (%EL)');
        xlim([-30 20]);
        ylim([-30 20]);
        zlim([-30 -.1]);
        set(gcf,'Position',[680   675   414   303]);
        
    for modelidx = 1:numel(model_range)
        model = model_range{modelidx};
        beta_best_ = dat(modelidx).beta_best;
        dat(modelidx).newfun(beta_best_)
        [fun_dx, fun_eval] = find_displacement(beta_best_,dat(modelidx).fun,dat(modelidx).y,dat(modelidx).x,[-35:30],ratio);
        figure(108)
        %subplot(121);
        %plot(ax,fun_eval);
        subplot(1,3,1);
        hold on;
        plot3([-35:30],fun_dx,[-35:30]*0+1000,'LineStyle','--','LineWidth',2,'Display',model);
        [fun_dx, fun_eval] = find_displacement(beta_best_,dat(modelidx).fun,dat(modelidx).y,dat(modelidx).x,pos_range,ratio);
        subplot(1,3,2);
        hold on;
        semilogy(pos_range,fun_eval/fun_eval(1));
        subplot(1,3,3);
        hold on;
        plot(pos_range,fun_dx,'LineStyle','--','LineWidth',2,'Display',model);

        for i=1:numel(compare_list)/2
            figure(120+i);
            mItmp = mI_rec{i,1};
            postmp = pos_rec{i,1}+fun_dx;
            plot(postmp,mItmp,'LineStyle','--','Display',model);
        end
        figure(107);
        hold on;
        h2 = plot3(pos_range,fun_dx,pos_range*0-1,'LineStyle','--','LineWidth',2,'Display',model,'color','g');
    end
    plot3(pos_range,pos_range*0,pos_range*0-1,'LineStyle','-','color','k','LineWidth',1);
    caxis([-10 -1]);
    %legend([h2,h1],{'|\Delta(\it{x)|','Theory'},'box','off')
    box on;
    figure(108);
    subplot(1,3,1);legend show;
    subplot(1,3,2);legend show;
    subplot(1,3,3);legend show;
    for i=1:numel(compare_list)/2
       figure(120+i);
       legend show;
    end
%% If impose fit - shifted boundary
if impose_fit
    for i = 1:numel(compare_list)/2
        figure(40);
        subplot(1,numel(compare_list)/2,i);
        plot(pos_rec{i,1}-dat(modelidx).beta_best*log(2),mI_rec{i,1},'--k')        
    end
end

%% Save data:
if fit_lambda
    mkdir('shift_rec');
    filename = [num2str(compare_list) '_' num2str(isBcd1X) '_nc' num2str(nc_range)];
    save(['shift_rec/' filename]);
    figure(107);
    saveas(gcf,['shift_rec/logmap_' filename]);
    figure(40);
    saveas(gcf,['shift_rec/pattern_' filename]);
end