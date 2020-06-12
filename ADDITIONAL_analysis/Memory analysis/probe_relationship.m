function [pos_range,tau,pval,Nsample,F,H1,H2,I12,Ibg,pI12]=probe_relationship(x,y,fea1,fea2,wd,tscnt,label1,label2,pos_range)
% Probing for relationship between features of cell 1 (mother) and features
% of cell 2 (daughter) given the same position
% Input:
    % x: x position of cell2
    % y: y position of cell2
    % fea1: feature of cell1
    % fea2: feature of cell2
    % wd: width of the scanning window
    % tscnt (optional) time series separation
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
%% Defining pos_range:
    if ~exist('pos_range','var')
        pos_range=[-30:1:0];    % Unit in %EL
    end
%% Create data holder
    tau=zeros(3,numel(pos_range));      % tau - Spearman correlation test
    pval=zeros(3,numel(pos_range));     % pval - Spearman correlation test
    Nsample=zeros(1,numel(pos_range));  % Number of samples
    if ~exist('tscnt','var')
        tscnt=fea1*0+1;
    end
    mkdir('img_rec');
%% Binning the features: For mutual information only
    [~,ax1]=histcounts(fea1,0:0.2:1);
    fea1_=interp1(ax1,ax1,fea1,'nearest','extrap');
    [~,ax2]=hist(fea2,5);
    fea2_=interp1(ax2,ax2,fea2,'nearest','extrap');
    H1=zeros(1,numel(pos_range));
    H2=zeros(1,numel(pos_range));
    I12=zeros(1,numel(pos_range));
    Ibg=zeros(1,numel(pos_range));
    pI12=ones(1,numel(pos_range));
%% Moving the scan window and test the indepence between cells' features
    cnt=0;
    hold off;
    for pos=pos_range
        cnt=cnt+1;
        tmp=(x(:)-pos+wd/2).*(x(:)-pos-wd/2)<=0;
        % Probe for relationship between mother and daughter features
        if sum(tmp)>0
            [r,p]=corr(fea1(tmp)',fea2(tmp)','type','Spearman');
        else
            r=NaN;
            p=NaN;
        end
        tau(1,cnt)=r;
        pval(1,cnt)=p;
        Nsample(1,cnt)=sum(tmp);
        subplot(131);
        if p<0.05
            plot_marker(fea1(tmp),fea2(tmp),tscnt(tmp),'r');
        else
            plot_marker(fea1(tmp),fea2(tmp),tscnt(tmp),'b');
        end
        title(['pos: ' num2str(pos) ', tau: ' num2str(r,'%2.2g') ', p-val: ' num2str(p,'%2.3g') ]);
        xlabel(label1);
        ylabel(label2);
        if numel(fea1)*numel(fea2)>0
            axis([0 max(fea1)*1.1+0.1 0 max(fea2)*1.1+0.1]);
        end
        % Probe for relationship between daughter features and x position (test if the wd is small enough)
        if sum(tmp)>0
            [r,p]=corr(x(tmp)',fea2(tmp)','type','Spearman');
        else
            r=NaN;
            p=NaN;
        end
        tau(2,cnt)=r;
        pval(2,cnt)=p;
        subplot(132);
        if p<0.05
            plot_marker(x(tmp),fea2(tmp),tscnt(tmp),'r');
        else 
            plot_marker(x(tmp),fea2(tmp),tscnt(tmp),'b');
        end
        title(['pos: ' num2str(pos) ', tau: ' num2str(r,'%2.2g') ', p-val: ' num2str(p,'%2.3g') ]);
        xlabel('x');
        ylabel(label2);
        if numel(fea1)*numel(fea2)>0
            axis([pos-wd/2 pos+wd/2 0 max(fea2)*1.1+0.1]);
        end        
        legend(['N=' num2str(sum(tmp))]);
        
        % Probe for relationship between daughter features and y position
        if sum(tmp)>0
            [r,p]=corr(y(tmp)',fea2(tmp)','type','Spearman');
        else
            r=NaN;
            p=NaN;
        end
        tau(3,cnt)=r;
        pval(3,cnt)=p;
        subplot(133);
        if p<0.05
            plot_marker(y(tmp),fea2(tmp),tscnt(tmp),'r');
        else 
            plot_marker(y(tmp),fea2(tmp),tscnt(tmp),'b');
        end
        title(['pos: ' num2str(pos) ', tau: ' num2str(r,'%2.2g') ', p-val: ' num2str(p,'%2.3g') ]);
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
        
        %% Calculate the entropy of fea1 and fea2 and mutual information
        H1(cnt)=entropy_(fea1_(tmp));
        H2(cnt)=entropy_(fea2_(tmp));
        [I12(cnt),Ibg(cnt),pI12(cnt)]=randmi(fea1_(tmp),fea2_(tmp));
    end