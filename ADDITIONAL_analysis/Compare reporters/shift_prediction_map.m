function [pos_prediction_map,ax,ay] = shift_prediction_map(pos_range1,mIin1, sIin1,pos_range2,mIin2, sIin2)
    
    pos_range1_ = pos_range1(:);pos_range1_ = pos_range1_(~isnan(mIin1));
    sIin1_ = sIin1(:);sIin1_ = sIin1_(~isnan(mIin1));sIin1_(sIin1_==0)=max(sIin1_)/20;
    mIin1_ = mIin1(:);mIin1_ = mIin1_(~isnan(mIin1));
    
    pos_range2_ = pos_range2(:);pos_range2_ = pos_range2_(~isnan(mIin2));
    sIin2_ = sIin2(:);sIin2_ = sIin2_(~isnan(mIin2));sIin2_(sIin2_==0)=max(sIin2_)/20;
    mIin2_ = mIin2(:);mIin2_ = mIin2_(~isnan(mIin2));

    % Predicting expression probabilty: discretize expression level first
    ax_I = linspace(-max([sIin1_*3;sIin2_*3]),max([mIin1_+sIin1_*3;mIin2_+sIin2_*3]),1000);
    pIin_all = ax_I*0;
    for i = 1:numel(mIin2_)
        pIin_all = pIin_all + normpdf(ax_I,mIin2_(i),sIin2_(i));
    end
    % Prepare the map
    pos_prediction_map = zeros(numel(mIin1_),numel(mIin2_));
    for i = 1:numel(mIin1_)
        for j=1:numel(mIin2_)
            tmp = ...
                    normpdf(ax_I,mIin1_(i),sIin1_(i)).*...
                    normpdf(ax_I,mIin2_(j),sIin2_(j))./...
                    pIin_all;
            pos_prediction_map(i,j) = nansum(tmp);
        end
        pos_prediction_map(i,:)=pos_prediction_map(i,:)/sum(pos_prediction_map(i,:));
    end
    [ax,ay] = meshgrid(pos_range1_,pos_range2_);
    figure;
    subplot(121);
    HeatMap_(pos_prediction_map',ax,ay,[0 max(pos_prediction_map(:))]);
    hold on;
    plot3(pos_range1_(:),pos_range1_(:),pos_range1_(:)*0+1,'LineStyle','--','color','w','LineWidth',2);
    set(gca,'Ydir','normal');
    xlabel('Original position x (%EL)');
    ylabel('Predicted position x'' (%EL)');
    subplot(122);
    HeatMap_(pos_prediction_map',ax,ay-ax,[0 max(pos_prediction_map(:))]);
    hold on;
    plot3(pos_range1_(:),pos_range1_(:)*0,pos_range1_(:)*0+1,'LineStyle','--','color','w','LineWidth',2);
    set(gca,'Ydir','normal');
    xlabel('Original position x (%EL)');
    ylabel('Predicted shift x''-x (%EL)');
end