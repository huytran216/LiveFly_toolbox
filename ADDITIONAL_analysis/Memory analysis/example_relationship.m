function [pos_range,tau,pval,Nsample,F,H1,H2,I12,Ibg,pI12]=probe_relationship(x,y,fea1,fea2,wd,tscnt,label1,label2)
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
    pos_range=[-10 -2 3];
    color='rgb';
    tau=zeros(3,numel(pos_range));
    pval=zeros(3,numel(pos_range));
    Nsample=zeros(1,numel(pos_range));
    if ~exist('tscnt','var')
        tscnt=fea1*0+1;
    end
    mkdir('img_rec');
%% Binning the features:
    [~,ax1]=hist(fea1,5);
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
    figure;
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
        
        plot_marker(fea1(tmp),fea2(tmp),tscnt(tmp),color(cnt));
        leg{cnt}=(['position: ' num2str(pos) '%, p-value ~ ' num2str(p,'%2.3g') ]);
        xlabel(label1);
        ylabel(label2);
        if numel(fea1)*numel(fea2)>0
            axis([0 max(fea1)*1.1+0.1 0 max(fea2)*1.1+0.1]);
        end
hold on;
        % Record the images
        tightfig;
        F(cnt)=getframe(gcf);
        
        %% Calculate the entropy of fea1 and fea2 and mutual information
        H1(cnt)=entropy_(fea1_(tmp));
        H2(cnt)=entropy_(fea2_(tmp));
        [I12(cnt),Ibg(cnt),pI12(cnt)]=randmi(fea1_(tmp),fea2_(tmp));
    end
    box off;
    h=legend(leg);
    set(h,'Location','NorthWest');