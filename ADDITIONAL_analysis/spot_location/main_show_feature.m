%% Plot features:
        % First spot position
        % Last spot position
        % Plot all spot position
        % Anisothropy plot
datafolder=['..\Data analysis\'];
load(fullfile(datafolder,'data\WT_CP.mat'));
load(fullfile(datafolder,'feature_label.mat'));
%% Define the plot parameters:
cycleno=13;
ts_spec=find(arrayfun(@(x) getfield(DatasetList,{x},['nc' num2str(cycleno)]),1:Nmov));

color='crgb';
%% Initiating the feature storage
    xrec=[];
    fxpos=[];
    fypos=[];
    lxpos=[];
    lypos=[];
    fea=[];
%% Plot spot position
    for tsidx=ts_spec
        % Take cell id
        idselect=find(([datamat(:).cycle]==cycleno)&(([datamat(:).tscnt]==tsidx)));
        % Take cell position
        xaxis=[datamat(idselect).x]*100-50; % No embryo alignment required
        % Begining probing
        cnt=0;
        for id=idselect
            cnt=cnt+1;
            % Calculate first spot position
            first=find(datamat(id).Intensity,1,'first');
            last=find(datamat(id).Intensity,1,'last');
            % If spots are found
            if numel(first)
                fxpos=[fxpos;datamat(id).xspotrec(first)/datamat(id).sizerec(first)];
                fypos=[fypos;datamat(id).yspotrec(first)/datamat(id).sizerec(first)];
                lxpos=[lxpos;datamat(id).xspotrec(last)/datamat(id).sizerec(last)];
                lypos=[lypos;datamat(id).yspotrec(last)/datamat(id).sizerec(last)];
                fea=[fea;datamat(id).Feature];
            else
                fxpos=[fxpos;NaN];
                fypos=[fypos;NaN];
                lxpos=[lxpos;NaN];
                lypos=[lypos;NaN];
                fea=[fea;datamat(id).Feature];
            end
        end
        % Record cell position
        xrec=[xrec xaxis];
    end
    %% Correct for drift?
    fxpos=fxpos-nanmean(fxpos);
    fypos=fypos-nanmean(fypos);
    lxpos=lxpos-nanmean(lxpos);
    lypos=lypos-nanmean(lypos);
    %% Plot the results
    figure;
    subplot(221);
    plot(fxpos,fypos,'x');hold on;
    plot(sin(linspace(0,2*pi,100)),cos(linspace(0,2*pi,100)));
    plot([-1 1],[0 0]);
    plot([0 0],[-1 1]);
    axis([-1.1 1.1 -1.1 1.1]);
    title('First spot position');
    subplot(222);
    plot(lxpos,lypos,'x');hold on;    
    plot(sin(linspace(0,2*pi,100)),cos(linspace(0,2*pi,100)));
    plot([-1 1],[0 0]);
    plot([0 0],[-1 1]);
    axis([-1.1 1.1 -1.1 1.1]);
    title('Last spot position');
    subplot(223)
    ax=[0:0.1:1];
    fdist=[sqrt(fxpos.^2+fypos.^2)];
    ldist=[sqrt(lxpos.^2+lypos.^2)];
    bf=hist(fdist,ax);
    bl=hist(ldist,ax);
    plot(ax,bf);hold on;plot(ax,bl);
    h=legend('first','last');
    set(h,'Location','NorthWest');
    xlabel('distance toward the center');
    ylabel('count ');
    title('Dist. spot distance to center');
    subplot(224)
    ddist=fdist-ldist;
    hist(ddist,-1:0.1:1);
    [~,p]=(ttest(ddist));
    title(['t-test. p-val ' num2str(p)]);
    xlabel('spot displacement');
    %% Plot the spot position vs Features
    figure;
    fea_id=2;
    plot(fdist,fea(:,[fea_id]),'x');
    xlabel('Spot Displacement');
    ylabel(feature_label(fea_id));
    [rho,val]=corr(ddist(~isnan(ddist)),fea(~isnan(ddist),[4]),'type','Spearman');
    rho
    val