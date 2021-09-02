function [pos_prediction_map,ax,ay,pos_prediction_map_,pos_range1,pos_range2] = shift_prediction_map(pos_range1,mIin1, sIin1,pos_range2,mIin2, sIin2)
    % Expand pos_range by patching data with maximum and mininum:
    pos_range1_ = [pos_range1-numel(pos_range1) pos_range1() pos_range1+numel(pos_range1)];
    pos_range2_ = [pos_range2-numel(pos_range2) pos_range2() pos_range2+numel(pos_range2)];
    [~,maxpt1] = max(mIin1); 
    mIin1(1:find(~isnan(mIin1),1,'first')-1)=mIin1(maxpt1);
    sIin1(1:find(~isnan(mIin1),1,'first')-1)=sIin1(maxpt1);
    mIin1(find(~isnan(mIin1),1,'last')+1:end)=0;
    sIin1(find(~isnan(mIin1),1,'last')+1:end)=0;
    mIin1 = [mIin1(maxpt1)*ones(1,numel(pos_range1)) mIin1 0*ones(1,numel(pos_range1))];
    sIin1 = [sIin1(maxpt1)*ones(1,numel(pos_range1)) sIin1 0*ones(1,numel(pos_range1))];
    [~,maxpt2] = max(mIin2);
    mIin2(1:find(~isnan(mIin2),1,'first')-1)=mIin2(maxpt2);
    sIin2(1:find(~isnan(mIin2),1,'first')-1)=sIin2(maxpt2);
    mIin2(find(~isnan(mIin2),1,'last')+1:end)=0;
    sIin2(find(~isnan(mIin2),1,'last')+1:end)=0;
    mIin2 = [mIin2(maxpt2)*ones(1,numel(pos_range2)) mIin2 0*ones(1,numel(pos_range2))];
    sIin2 = [sIin2(maxpt2)*ones(1,numel(pos_range2)) sIin2 0*ones(1,numel(pos_range2))];
    

    sIin1_ = sIin1(:);sIin1_(sIin1_<max(sIin1_)/10)=max(sIin1_)/10;
    mIin1_ = mIin1(:);
    
    sIin2_ = sIin2(:);sIin2_(sIin2_<max(sIin2_)/10)=max(sIin2_)/10;
    mIin2_ = mIin2(:);

    % Predicting expression probabilty: discretize expression level first
    ax_I = linspace(-max([sIin1_*3;sIin2_*3]),max([mIin1_+sIin1_*3;mIin2_+sIin2_*3]),1000);
    % Prepare matrix G:
    matI2 = zeros(numel(ax_I)+1,numel(mIin2_));
    for j=1:numel(mIin2_)
        matI2(1:numel(ax_I),j) = normpdf(ax_I,mIin2_(j),sIin2_(j));
    end
    matI2(end,:)=ones(1,numel(mIin2_));
    % Prepare the map
    pos_prediction_map = zeros(numel(mIin1_),numel(mIin2_));
    for i = 1:numel(mIin1_)
        for j=1:numel(mIin2_)
            tmp_sum = (sIin1_(i)^2+sIin2_(j)^2);
            tmp_coeff = sIin1_(i)^2*sIin2_(j)^2/tmp_sum;
            
            pos_prediction_map(i,j) = exp(-(tmp_coeff*(mIin1_(i)-mIin2_(j))^2/tmp_sum + ...
                tmp_coeff*log(2*pi*(tmp_sum)))/2/tmp_coeff);
        end
        pos_prediction_map(i,:)=pos_prediction_map(i,:)/sum(pos_prediction_map(i,:));
    end
    [ax,ay] = meshgrid(pos_range1_,pos_range2_);
    subplot(121);
    HeatMap_(pos_prediction_map',ax,ay,[0 max(pos_prediction_map(:))]);
    hold on;
    plot3(pos_range1_(:),pos_range1_(:),pos_range1_(:)*0+1,'LineStyle','--','color','w','LineWidth',2);
    set(gca,'Ydir','normal');
    xlabel('Original position x (%EL)');
    ylabel('Predicted position x'' (%EL)');
    subplot(122);
    
    pos_prediction_map_ = zeros(numel(pos_range1),numel(pos_range2));
    mapped = zeros(numel(pos_range1),numel(pos_range2));
    for i=1:numel(pos_range1)
        for j=1:numel(pos_range2)
            postmp = (ax'==pos_range1(i))&((ay-ax)'==pos_range2(j));
            if sum(postmp(:))==1
                pos_prediction_map_(i,j) = pos_prediction_map(postmp);
                mapped(i,j)=1;
            end
        end
        pos_prediction_map_(i,:) = pos_prediction_map_(i,:)/sum(pos_prediction_map_(i,:));
    end
    
    % Plot the results:
    %HeatMap_(pos_prediction_map',ax,ay-ax,[0 max(pos_prediction_map(:))]);
    HeatMap_(pos_prediction_map_',pos_range1,pos_range2,[0 max(pos_prediction_map(:))]);
    hold on;
    plot3(pos_range1(:),pos_range1(:)*0,pos_range1(:)*0+1,'LineStyle','--','color','k','LineWidth',1);
    set(gca,'Ydir','normal');
    colormap(flipud(gray))
    xlabel('Original position (%EL)');
    ylabel('Predicted shift (%EL)');
    xlim([pos_range1(1) pos_range1(end)]);
    ylim([pos_range2(1) pos_range2(end)]);
end