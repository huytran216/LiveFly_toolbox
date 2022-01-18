function [pos_range,tau,pval,Nsample,Nsample_ts,F,tau_mix_within,pval_mix_within,tau_mix_all,pval_mix_all]=probe_relationship(x,y,fea1,fea2,wd,tscnt,label1,label2,pos_range,plot_sliding_window,Nmix_within,Nmix_all)
% Probing for relationship between features of cell 1 (mother) and features
% of cell 2 (daughter) given the same position
% Input:
    % x: x position of cell2
    % y: y position of cell2
    % fea1: feature of cell1
    % fea2: feature of cell2
    % wd: width of the scanning window
    % tscnt (optional) time series separation
    % toplot: plot the result per sliding window or not
    % Nmix_within: mix the daughters within individual embryos
    % Nmix_all: mix the daughters for all embryos
% Output
    % pos_range: scanned position along AP axis
    % tau: [3 x pos_range] contain the spearman rank correlation between
    %   fea1 and fea2
    %   x and fea2
    %   y and fea2
    % pval: [3 x pos_range] contain the p-value of spearman test between
    %   fea1 and fea2
    %   x and fea2
    %   y and fea2
    % NSample: number of samples in each position
    % F: Figure images plotting Fea_m vs Fea_d
    % Hm: Entropy of mother features
        % Require the binning of data
    % Hd: Entropy of daughter features
        % Require the binning of data
    % Imd: Mutual information between mother and daughter features
        % Require the binning of data
%% Defining pos_range and Nmix 
    if ~exist('pos_range','var')
        pos_range=[-30:1:0];    % Unit in %EL
    end
    if ~exist('Nmix_within','var')
        Nmix_within=0;
    end
    if ~exist('Nmix_all','var')
        Nmix_all=0;
    end
%% Create data holder
    tau=zeros(3,numel(pos_range));      % tau - Spearman correlation test
    pval=zeros(3,numel(pos_range));     % pval - Spearman correlation test
    Nsample=zeros(1,numel(pos_range));  % Number of samples
    ts_range = sort(unique(tscnt));
    Nsample_ts=zeros(max(ts_range),numel(pos_range));   % Number of samples per embryos
    tau_mix_within = zeros(Nmix_within,numel(pos_range));
    pval_mix_within = tau_mix_within;
    tau_mix_all = zeros(Nmix_all,numel(pos_range));
    pval_mix_all = tau_mix_all;
    if ~exist('tscnt','var')
        tscnt=fea1*0+1;
    end
%% Moving the scan window and test the indepence between cells' features
    cnt=0;
    hold off;    
    
    for pos=pos_range
        cnt=cnt+1;
        tmp=find((x(:)-pos+wd/2).*(x(:)-pos-wd/2)<=0);
        Nsample(1,cnt)=numel(tmp);
        Nsample_ts(ts_range,cnt)= histc(tscnt(tmp),ts_range);
        % Probe for relationship between mother and daughter features
        if Nsample(1,cnt)>0
            [r,p]=corr(fea1(tmp)',fea2(tmp)','type','Spearman');
        else
            r=NaN;
            p=NaN;
        end
        tau(1,cnt)=r;
        pval(1,cnt)=p;
        
        % Probe for relationship between mother and daughter features for
        % mixed samples
        if Nsample(1,cnt)>0
            for i=1:Nmix_within
                newidx = randperm_subset(tscnt(tmp));
                [r,p]=corr(fea1(tmp)',fea2(tmp(newidx))','type','Spearman');
                tau_mix_within(i,cnt)=r;
                pval_mix_within(i,cnt)=p;
            end
            for i=1:Nmix_all
                newidx = randperm(Nsample(1,cnt));
                [r,p]=corr(fea1(tmp)',fea2(tmp(newidx))','type','Spearman');
                tau_mix_all(i,cnt)=r;
                pval_mix_all(i,cnt)=p;
            end
        else
            tau_mix_within(:,cnt)=NaN;
            pval_mix_within(:,cnt)=NaN;
        end
        % Probe for relationship between daughter features and x position (test if the wd is small enough)
        if Nsample(1,cnt)>0
            [r,p]=corr(x(tmp)',fea2(tmp)','type','Spearman');
        else
            r=NaN;
            p=NaN;
        end
        tau(2,cnt)=r;
        pval(2,cnt)=p;
        
        % Probe for relationship between daughter features and y position
        if Nsample(1,cnt)>0
            [r,p]=corr(y(tmp)',fea2(tmp)','type','Spearman');
        else
            r=NaN;
            p=NaN;
        end
        tau(3,cnt)=r;
        pval(3,cnt)=p;
        
        % Plot of needed
        if plot_sliding_window
            
            subplot(131);
            if p<0.05
                plot_color(fea1(tmp),fea2(tmp),tscnt(tmp),'x');
            else
                plot_color(fea1(tmp),fea2(tmp),tscnt(tmp),'o');
            end
            title(['pos: ' num2str(pos) ', tau: ' num2str(tau(1,cnt),'%2.2g') ', p-val: ' num2str(pval(1,cnt),'%2.3g') ]);
            xlabel(label1);
            ylabel(label2);
            if numel(fea1)*numel(fea2)>0
                axis([0 max(fea1)*1.1+0.1 0 max(fea2)*1.1+0.1]);
            end

            subplot(132);
            if p<0.05
                plot_color(x(tmp),fea2(tmp),tscnt(tmp),'x');
            else 
                plot_color(x(tmp),fea2(tmp),tscnt(tmp),'o');
            end
            title(['pos: ' num2str(pos) ', tau: ' num2str(tau(2,cnt),'%2.2g') ', p-val: ' num2str(pval(2,cnt),'%2.3g') ]);
            xlabel('x');
            ylabel(label2);
            if numel(fea1)*numel(fea2)>0
                axis([pos-wd/2 pos+wd/2 0 max(fea2)*1.1+0.1]);
            end        
            legend(['N=' num2str(sum(tmp))]);


            subplot(133);
            if p<0.05
                plot_color(y(tmp),fea2(tmp),tscnt(tmp),'x');
            else 
                plot_color(y(tmp),fea2(tmp),tscnt(tmp),'o');
            end
            title(['pos: ' num2str(pos) ', tau: ' num2str(tau(3,cnt),'%2.2g') ', p-val: ' num2str(pval(3,cnt),'%2.3g') ]);
            xlabel('y');
            ylabel(label2);
            if numel(fea1)*numel(fea2)>0
                axis([0 max(y)*1.1 0 max(fea2)*1.1+0.1]);
            end
            
            % Check if ON - OFF is tested
            if strcmp(label1,'ON mother')&strcmp(label2,'ON daughter')
                subplot(131);
                fea_=fea1(tmp)+fea2(tmp)*2;
                hist(fea_,[0 1 2 3]);
                set(gca,'XTickLabel',{'00','10','10','11'});
                title('Mother Daughter count');
            end
            
            % Record the images
            tightfig;
            set(gcf, 'Position', [0 300 1200 300]);
            F(cnt)=getframe(gcf);
        else
            F = [];
        end
    end