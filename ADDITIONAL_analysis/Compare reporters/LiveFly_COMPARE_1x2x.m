%% Re-extract all features in all strains before running the script
% THe data should be in Data_analysis/
warning off;
%% Params
fld='../../Data analysis/';
load([fld 'feature_label.mat']);

Nsample_min=2;              % Minimum total nuclei per position
Nsample_per_embryo_min=2;   % Minimum nuclei per embryo per position
Nembryo_min = 2;            % Mininum number of embryo per position

plot_embryo_error=1;    % plot error based on embryo diversity (1) or nuclei error (merged, 0)
shaded_error_bar = 1;   % Plot shaded errorbar or normal errorbar

trimmed_trace = 1;      % trimmed trace (1) or not (0)
smooth_curve = 1;
normalize_intensity = 0;
compare_1x2x = false;
%% Set up data list

dtset = struct('filename','','label','');
dtset(1).filename = 'hb-vk33';  dtset(1).label = 'hb-P2 vk33';

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
compare_list = [2 3 2 3];isBcd1X=[zeros(1,numel(compare_list)/2) ones(1,numel(compare_list)/2)]; % For hb-B6-H6B6 comparison, 1x2x
%compare_list = [1 8 9]; isBcd1X =[0 0 0 ];% For vk33 vs random insertion
%compare_list = [7 7];isBcd1X=[0 1];
%compare_list = [1 2 10 12]; isBcd1X = compare_list*0;
%compare_list = [3 3];isBcd1X=[0 1];

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
fea_range=[16];
nc_range=[13];

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
                        ne_rec=ne_rec+(nf_indi{1,nc,j}>=Nsample_per_embryo_min);
                    end
                end
            % Calculate mean of indi curve:
                for j=1:size(nf_indi,3)
                    if numel(nf_indi{1,nc,j})
                        tmp=mf_indi{fea,nc,j};
                        tmp(isnan(tmp))=0;
                        nf_indi{1,nc,j}(isnan(nf_indi{1,nc,j}))=0;
                        mf_indi_=mf_indi_ + tmp.*(nf_indi{1,nc,j}>=Nsample_per_embryo_min);
                        sf_indi_=sf_indi_ + (tmp.*(nf_indi{1,nc,j}>=Nsample_per_embryo_min)).^2;
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
                if smooth_curve
                    stmp_tmp = isnan(stmp);
                    mtmp = smooth(mtmp(end:-1:1));
                    mtmp = mtmp(end:-1:1);
                    %stmp = smooth(stmp(end:-1:1));
                    %stmp = stmp(end:-1:1);
                    %stmp(stmp_tmp)=nan;
                    mtmp(stmp_tmp)=nan;
                end
            % Normalize by intensity
                if normalize_intensity
                    mtmp = mtmp/vborder(i,fea,nc,1);
                    stmp = stmp/vborder(i,fea,nc,1);
                end
            % Determine color based on 1x or 2x
                if isBcd1X(i)==1
                    color=corder(2);
                else
                    color=corder(4);
                end
            % Plot mean curve with error
                figure(20+original_i);
                set(gcf,'Position',[680   738   284   240]);
                title(DatasetLabel{original_i});
                subplot(numel(nc_range),numel(fea_range),cnt);
                if shaded_error_bar                    
                    h40(cnt2)=shadedErrorBar(pos_range(flttmp),mtmp,stmp,{'color',color,'Display',DatasetLabel{i}},0.8,1);
                else
                    h40(cnt2)=errorbar(pos_range(flttmp),mtmp,stmp,'Display',DatasetLabel{i},'color',color);
                end
                xlabel('AP axis (%EL)');
                ylabel(feature_label{fea});
                xlim(AP_limit);
                hold on;
            
            % Plot prediction of position based on expression alone
                figure(80+original_i);
                subplot(1,2,1+isBcd1X(i));
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
                
                mI_rec{original_i,1+isBcd1X(i)}=mIall;
                sI_rec{original_i,1+isBcd1X(i)}=sIall;
                pos_rec{original_i,1+isBcd1X(i)}=pos_range;
                
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
    [diff_plot{i},ax,ay,diff_grid{i}]=shift_prediction_map(pos_rec{i,1}(pos_range1),mI_rec{i,1}(pos_range1),sI_rec{i,1}(pos_range1),...
        pos_rec{i,2}(pos_range2),mI_rec{i,2}(pos_range2),sI_rec{i,2}(pos_range2));
    
    %subplot(121);
    %xlim([-30 20]);
    %ylim([-30 20]);
    
    %subplot(122);
    
    set(gcf,'Position',[680   738   284   240]);
    xlim([-30 20]);
    ylim([-30 20]);
    caxis([0 0.5]);
    colormap(flipud(gray));
    if i==1
        diff_all = diff_plot{1};
        diff_grid_all = diff_grid{1};
    else
        diff_all = diff_all+diff_plot{i};
        diff_grid_all = diff_grid_all + diff_grid{1};
    end
    figure(6);
    subplot(1,numel(compare_list)/2,i);
    mDX = [];
    sDX = [];
    for j=1:sum(pos_range1)
        %mDX(j) = sum(diff_plot{i}(j,:)'.*(ay(:,j)-ax(:,j)));
        %sDX(j) = sqrt(sum(diff_plot{i}(j,:)'.*((ay(:,j)-ax(:,j)).^2)) - mDX(j).^2);
        mDX(j) = sum(diff_grid{i}(j,:).*pos_range(pos_range1));
        sDX(j) = sqrt(sum(diff_grid{i}(j,:).*pos_range(pos_range1).^2) - mDX(j).^2);
    end
    errorbar(pos_range(pos_range1),mDX,sDX);
    xlim([-30 20]);
    title(DatasetLabel{i});
end
diff_all = diff_all/i;
% %% Summary
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
%     
% % subplot(122);
% %     mDX = [];
% %     sDX = [];
% %     for j=1:numel(pos_range)
% %         % Cut the tail in the distribution: take 95% of mass distribution
% %         [maxpro,midpos] = max(diff_all(j,:));
% %         diff_all(j,diff_all(j,:)<maxpro/10)=0;
% %         diff_all(j,:)=diff_all(j,:)/sum(diff_all(j,:));
% %         
% %         mDX(j) = sum(diff_all(j,:)'.*(ay(:,j)-ax(:,j)));
% %         sDX(j) = sqrt(sum(diff_all(j,:)'.*((ay(:,j)-ax(:,j)).^2)) - mDX(j).^2);
% %     end
% %     errorbar(pos_range,mDX,sDX);
% %     xlabel('Original position x (%EL)');
% %     ylabel('Predicted shift x''-x (%EL)');
% %     xlim([-30 20]);
% %     ylim([-20 0]);
%     %set(gcf,'Position',tmpfig);
%% Fit displacement with a curve:
    syms x;
    model = 'hybrid';
    switch model 
        case 'exp'
            y = sym('y'); fun = exp(-x/y);  % Exponential gradient
        case 'hybrid'
            y = sym('y',[1 3]); fun = exp(-x/y(1)) + y(2)*exp(-x/y(3));  % Two Exponential gradient
        case 'power'
            y = sym('y', [1 2]);fun  = (x+y(1))^(-y(2)); % Power gradient
    end
    
    ax_ = ax(1,:);
    ax_fit = [-25:0];
    ratio = 0.5;
    
    newfun = @(y0) -sum(log(1e-10+diff_all(sub2ind(size(diff_all),...
        -ax_(1)+1+ax_fit,...
        max([ones(1,numel(ax_fit)); ...
        min([numel(ax_)*ones(1,numel(ax_fit)); ...
        (-ax_(1)+1+ax_fit+round(find_displacement(y0,fun,y,x,ax_fit,ratio)))
        ])...
        ])...
        ))));
    switch model 
        case 'exp'
            [beta_best,f0]=fminsearchbnd(newfun,10,1,100);
        case 'hybrid'
            [beta_best,f0]=fminsearchbnd(newfun,[5 3 20],[1 0 1],[100 10 100]);
        case 'power'
            [beta_best,f0]=fminsearchbnd(newfun,[35 2],[-ax_fit(1)+1 0.1],[1000 100]);
    end
%% Plot the results
beta_best_ = [beta_best];
newfun(beta_best_)
[fun_dx, fun_eval] = find_displacement(beta_best_,fun,y,x,[-35:30],ratio);
figure(107)
%subplot(121);
%plot(ax,fun_eval);
subplot(1,2,1);
hold on;
plot3([-35:30],fun_dx,[-35:30]*0+1000,'LineStyle','--','color','k','LineWidth',2);
%
figure(108);
[fun_dx, fun_eval] = find_displacement(beta_best_,fun,y,x,pos_range,ratio);
subplot(1,2,1);
hold on;
semilogy(pos_range,fun_eval/fun_eval(1));
subplot(1,2,2);
hold on;
plot(pos_range,fun_dx,'LineStyle','--','color','k','LineWidth',2);

%% Transform the axes:
for i=1:numel(compare_list)/2
    figure(120+i);
    mItmp = mI_rec{i,1};
    postmp = pos_rec{i,1}+fun_dx;
    plot(pos_rec{i,1},mI_rec{i,1},'-b');hold on;
    plot(pos_rec{i,2},mI_rec{i,2},'-r');
    plot(postmp,mItmp,'--k');
    xlim([AP_limit])
end

