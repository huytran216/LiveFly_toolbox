%% Plot stuffs
ymax=zeros(1,Nfea);
ymin=zeros(1,Nfea)+1e5;
rec_fea=[];
h70 = {};

header={};
cnti=2; % Line number
h=[];
hplot = [];
xborder = []; hborder = []; vborder = []; noembryo = [];
cnt=0;
cnt2=0;
sample_f = {};
sample_x = {};
for nc=nc_range
    for fea=fea_range
        cnt=cnt+1;
        RHO = [];
        PVAL = [];
        NSAMPLE = [];
        for i=1:numel(compare_list)
            cnt2=cnt2+1;
            load(fullfile(fld,folder{trimmed_trace+1},DatasetFile{i}),'pos_range','mf_rec','nf_rec','sf_rec','nf_indi','mf_indi','sf_indi','FitRes','samplef_rec','samplex_rec');
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
                subplot(numel(fea_range),numel(nc_range),cnt);
                color_order = get(gca,'colororder');
                set(gcf,'Position',[100   100   400*numel(nc_range)   250*numel(fea_range)]);
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
                    h40(cnt2)=shadedErrorBar(pos_range(flttmp),mtmp,stmp,{'color',color,'Display',DatasetLabel{i},'LineWidth',2},0.8,1);
                else
                    h40(cnt2)=errorbar(pos_range(flttmp),mtmp,stmp,'Display',DatasetLabel{i});
                end
                hold on;
            % Plot mean curve with error, normalized by maximum value
                figure(41);
                subplot(numel(fea_range),numel(nc_range),cnt);
                set(gcf,'Position',[100   100   400*numel(nc_range)   250*numel(fea_range)]);
                if shaded_error_bar
                    if compare_1x2x & numel(compare_list)==2
                        if cnt2==1
                            color=corder(2);
                        else
                            color=corder(4);
                        end
                    end
                    h41(cnt2)=shadedErrorBar(pos_range(flttmp),mtmp/vborder(i,fea,nc,2),stmp/vborder(i,fea,nc,2),{'color',color,'Display',dtset(i).label,'LineWidth',2},0.8,1);
                else
                    h41(cnt2)=errorbar(pos_range(flttmp),mtmp/vborder(i,fea,nc,2),stmp/vborder(i,fea,nc,2),'Display',dtset(i).label);
                end
                hold on;
            % Plot mean curve without error
                figure(42);
                subplot(numel(fea_range),numel(nc_range),cnt);
                set(gcf,'Position',[100   100   400*numel(nc_range)   250*numel(fea_range)]);
                plot(pos_range(flttmp),mtmp,'Display',DatasetLabel{i},'LineWidth',2);
                hold on;
            % Put all features into one (with noise)    
                figure(70+i);
                    subplot(1,numel(nc_range),find(nc_range==nc));
                    set(gcf,'Position',[300   100   400*numel(nc_range)   250]);
                    if shaded_error_bar
                        h70{i,fea}=shadedErrorBar(pos_range(flttmp),mtmp,stmp,{'Display',feature_label{fea},'LineWidth',2},0.8,1);
                    else
                        h70{i,fea}=errorbar(pos_range(flttmp),mtmp,stmp,'Display',feature_label{fea});
                    end
                    ylim([0 1]);
                    hold on;
                    title(DatasetLabel{i});
                    xlim(AP_limit);
            % Put all features into one (without noise)
                figure(60+i);
                    subplot(1,numel(nc_range),find(nc_range==nc));
                    set(gcf,'Position',[500   100   400*numel(nc_range)   250]);
                    plot(pos_range(flttmp),mtmp,'Display',feature_label{fea},'LineWidth',2);
                    hold on;
                    title(DatasetLabel{i});
                    xlim(AP_limit);
            % Calculate Spearman testz
                sample_f{i,fea} = samplef_rec{nc-8,fea};
                sample_x{i,fea} = samplex_rec{nc-8,fea};
                idselect = (sample_x{i,fea}>=AP_limit(1))&(sample_x{i,fea}<=AP_limit(2));
                sample_f{i,fea} = sample_f{i,fea}(idselect);
                sample_x{i,fea} = sample_x{i,fea}(idselect);
                if numel(idselect)>10
                    [RHO(i),PVAL(i)] = corr(sample_f{i,fea}',sample_x{i,fea}','Type','Spearman');
                    NSAMPLE(i) = sum(idselect);
                else
                    RHO(i)=0;
                    PVAL(i) = 1;
                    NSAMPLE(i) = 0;
                end
             
        end
        % Get legends and correct axis
            figure(40);
            subplot(numel(fea_range),numel(nc_range),cnt);
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
            subplot(numel(fea_range),numel(nc_range),cnt);
            if nc==nc_range(end)
                xlabel('AP axis (%EL)');
            end
            ylabel(['Normalized ' feature_label{fea}]);
            
            figure(42);
            subplot(numel(fea_range),numel(nc_range),cnt);
            if nc==nc_range(end)
                xlabel('AP axis (%EL)');
            end
            ylabel(feature_label{fea});
            if cnt==1
                legend('Show');
            end
            xlim(AP_limit);ylim([0 1.3]);
        % Make a table out of it: constructs are put vertically
            figure(30);
            header{1,1}='FEATURE:';
            header{1,2}=feature_label{fea};
            header{2,1}='Construct';
            header{2,2}='N';
            header{2,3}='xborder';
            header{2,4}='Hill';
            header{2,5}='Width';
            header{2,6}='HalfMax';
            header{2,7}='Rho';
            header{2,8}='p-val';
            header{2,9}='Nsample';
            for i=1:numel(compare_list)
                cnti=cnti+1;
                header{cnti,1}=DatasetLabel{i};
                header{cnti,2}=num2str(noembryo(i,nc));
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
                header{cnti,7}= RHO(i);
                header{cnti,8}= PVAL(i);
                header{cnti,9}= NSAMPLE(i);
            end
            % Make the table:
            htmp=figure(30);set(htmp,'Name',['nc ' num2str(nc),', Feature: ' feature_label{fea}],'Position',[100 100 800 500]);
            htmptable =   uitable('Parent',htmp,'Position',[0 0 800 500]);
            htmptable.Data = header;
            set(htmptable,'ColumnEditable',true(1,10))
    end
    for i=1:numel(compare_list)
        figure(70+i);
        subplot(1,numel(nc_range),find(nc_range==nc));
        legend([h70{i,fea_range}],feature_label(fea_range));
    end
end


% Apply similar axis to all subplots in Fig 11 and Fig 12
for i=1:cnt
    figure(40);
    subplot(numel(fea_range),numel(nc_range),i);
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
    subplot(numel(fea_range),numel(nc_range),i);
    xlim(AP_limit);
    if rec_fea(i)==1
        ylim([0 1.1])
    else
        ylim([ymin(rec_fea(i)) ymax(rec_fea(i))]);
    end
    set(gca,'XTick',round((AP_limit(1):5:AP_limit(2))/5)*5);
end