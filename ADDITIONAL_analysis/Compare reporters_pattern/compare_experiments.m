%% Plot stuffs
ymax=zeros(1,Nfea);
ymin=zeros(1,Nfea)+1e5;
rec_fea=[];

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
            load(fullfile(fld,folder{istrimed_range(fea_range == fea)+1},DatasetFile{i}),...
                'pos_range','mf_rec','nf_rec','sf_rec','nf_indi','mf_indi','sf_indi','FitRes','samplef_rec','samplex_rec');
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
                        ne_rec=ne_rec+(nf_indi{1,nc,j}>=Nsample_indi);
                    end
                end
            % Calculate mean of indi curve:
                for j=1:size(nf_indi,3)
                    if numel(nf_indi{1,nc,j})
                        tmp=mf_indi{fea,nc,j};
                        tmp(isnan(tmp))=0;
                        nf_indi{1,nc,j}(isnan(nf_indi{1,nc,j}))=0;
                        mf_indi_=mf_indi_ + tmp.*(nf_indi{1,nc,j}>=Nsample_indi);
                        sf_indi_=sf_indi_ + (tmp.*(nf_indi{1,nc,j}>=Nsample_indi)).^2;
                    end
                end
                mf_indi_=mf_indi_./ne_rec;
                sf_indi_=sqrt(sf_indi_./ne_rec - mf_indi_.^2);
            % Plot mean curve with error
            
                figure(40);
                subplot(numel(nc_range),numel(fea_range),cnt);
                color_order = get(gca,'colororder');
                set(gcf,'Position',[100   100   400*numel(fea_range)   250*numel(nc_range)]);
                flttmp = nf_rec{fea,nc}>Nsample_min;
                if plot_embryo_error
                    mtmp = mf_indi_(flttmp);
                    stmp = sf_indi_(flttmp)./sqrt(ne_rec(flttmp)-1);
                else
                    mtmp = mf_rec{fea,nc}(flttmp);
                    stmp = sf_rec{fea,nc}(flttmp)./sqrt(nf_rec{fea,nc}(flttmp)-1);
                end
                if smooth_curve
                    stmp_tmp = isnan(mtmp);
                    mtmp = smooth(mtmp(end:-1:1));
                    mtmp = mtmp(end:-1:1);
                    %stmp = smooth(stmp(end:-1:1));
                    %stmp = stmp(end:-1:1);
                    %stmp(stmp_tmp)=nan;
                    mtmp(stmp_tmp)=nan;
                end
                color = corder(compare_list(i));
                if ~any(flttmp)
                    flttmp(1:end) = true;
                    mtmp = zeros(1,numel(flttmp));
                    stmp = zeros(1,numel(flttmp));
                end
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
                subplot(numel(nc_range),numel(fea_range),cnt);
                set(gcf,'Position',[100   100   400*numel(fea_range)   250*numel(nc_range)]);
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
                subplot(numel(nc_range),numel(fea_range),cnt);
                set(gcf,'Position',[100   100   400*numel(fea_range)   250*numel(nc_range)]);
                plot(pos_range(flttmp),mtmp,'Display',DatasetLabel{i},'LineWidth',2);
                hold on;
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
            subplot(numel(nc_range),numel(fea_range),cnt);
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
            subplot(numel(nc_range),numel(fea_range),cnt);
            if nc==nc_range(end)
                xlabel('AP axis (%EL)');
            end
            ylabel(['Normalized ' feature_label{fea}]);
            
            figure(42);
            subplot(numel(nc_range),numel(fea_range),cnt);
            if nc==nc_range(end)
                xlabel('AP axis (%EL)');
            end
            ylabel(feature_label{fea});
            if cnt==1
                legend('Show');
            end
            xlim(AP_limit);ylim([0 1.3]);
        % Plot fitted params (Boundary position,Hill,maximum expression) with confidence error bar:
            if compare_1x2x
                barplot_fun = @barplot_bound_1x2x;
            else
                barplot_fun = @barplot_bound;
            end
            
            figure(43);            
            subplot(numel(nc_range),numel(fea_range),cnt);
            set(gcf,'Position',[700   100   400*numel(fea_range)   250*numel(nc_range)]);
            btmp=barplot_fun(xborder(:,fea,nc,1),xborder(:,fea,nc,2),xborder(:,fea,nc,3),DatasetLabel,AP_limit(1),'h');
            box on;
            xlim(AP_limit);
            xlabel('AP axis (%EL)');
            xtickangle(0);
            
            figure(44);
            subplot(numel(nc_range),numel(fea_range),cnt);
            set(gcf,'Position',[700   100   400*numel(fea_range)   250*numel(nc_range)]);
            btmp=barplot_fun(hborder(:,fea,nc,1),hborder(:,fea,nc,2),hborder(:,fea,nc,3),DatasetLabel,0);
            box on;
            xtickangle(0);
            ylabel('Hill coefficient');
            
            figure(45);
            v1=-log(1/19)./hborder(:,fea,nc,1)*20*2;
            v2=-log(1/19)./hborder(:,fea,nc,2)*20*2;
            v3=-log(1/19)./hborder(:,fea,nc,3)*20*2;
            subplot(numel(nc_range),numel(fea_range),cnt);
            set(gcf,'Position',[700   100   400*numel(fea_range)   250*numel(nc_range)]);
            btmp=barplot_fun(v1,v2,v3,DatasetLabel,0);
            box on;
            xtickangle(0);
            ylabel('Boundary width (%EL)');
            
            figure(46);
            subplot(numel(nc_range),numel(fea_range),cnt);
            set(gcf,'Position',[700   100   400*numel(fea_range)   250*numel(nc_range)]);
            btmp=barplot_fun(vborder(:,fea,nc,1),vborder(:,fea,nc,2),vborder(:,fea,nc,3),DatasetLabel,0);
            box on;
            xtickangle(0);
            ylabel('Maximum expression');
            % Full boundary position & Width
            
            figure(47);
            subplot(numel(nc_range),numel(fea_range),cnt);
            set(gcf,'Position',[700   100   400*numel(fea_range)   250*numel(nc_range)]);
            if compare_1x2x
                btmp=barplot_bound_width_1x2x(xborder(:,fea,nc,1),xborder(:,fea,nc,2),xborder(:,fea,nc,3),v1,v2,v3,DatasetLabel,[-50 20],'h');
            else
                btmp=barplot_bound_width(xborder(:,fea,nc,1),xborder(:,fea,nc,2),xborder(:,fea,nc,3),v1,v2,v3,DatasetLabel,[-50 20],'h');
            end
            box on;
            xtickangle(45);
            xlabel('Boundary position (%EL)');
        % Make a table out of it: constructs are put vertically
            figure(30);
            header={};
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
            cnti=2; % Line number
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
                        header{cnti,6}=[num2str(vborder(i,fea,nc,1),'%.1f') ' ' char(177) ' ' ...
                            num2str((vborder(i,fea,nc,3) - vborder(i,fea,nc,2)),'%.1f')];
                    case 5
                        header{cnti,6}=[num2str(vborder(i,fea,nc,1),'%.2f') ' ' char(177) ' ' ...
                            num2str((vborder(i,fea,nc,3) - vborder(i,fea,nc,2)),'%.2f')];
                    case 4
                        header{cnti,6}=[num2str(vborder(i,fea,nc,1),'%.2f') ' ' char(177) ' ' ...
                            num2str((vborder(i,fea,nc,3) - vborder(i,fea,nc,2)),'%.2f')];
                    otherwise
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
end

% Apply similar axis to all subplots in Fig 11 and Fig 12
for i=1:cnt
    figure(40);
    subplot(numel(nc_range),numel(fea_range),i);
    xlim(AP_limit);
    if rec_fea(i)==1
        ylim([0 1.1])
    else
        ylim([ymin(rec_fea(i)) ymax(rec_fea(i))]);
    end
    h=legend(h40(1:numel({dtset(compare_list).label})),DatasetLabel);
    legend boxoff
    set(gca,'XTick',round((AP_limit(1):10:AP_limit(2))/5)*5);    

    figure(41);
    subplot(numel(nc_range),numel(fea_range),i);
    xlim(AP_limit);
    ylim([0 1.1]);
    h=legend(h41(1:numel({dtset(compare_list).label})),DatasetLabel);
    legend boxoff
    set(gca,'XTick',round((AP_limit(1):10:AP_limit(2))/5)*5);   
    
    figure(42);
    subplot(numel(nc_range),numel(fea_range),i);
    xlim(AP_limit);
    if rec_fea(i)==1
        ylim([0 1.1])
    else
        ylim([ymin(rec_fea(i)) ymax(rec_fea(i))]);
    end
    legend boxoff
    set(gca,'XTick',round((AP_limit(1):5:AP_limit(2))/5)*5);
end