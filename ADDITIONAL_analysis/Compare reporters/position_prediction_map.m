function [pos_prediction_map] = position_prediction_map(pos_range,mIin, sIin)
    
    pos_range_ = pos_range(:);pos_range_ = pos_range_(~isnan(mIin));
    sIin_ = sIin(:);sIin_ = sIin_(~isnan(mIin));sIin_(sIin_==0)=max(sIin_)/20;
    mIin_ = mIin(:);mIin_ = mIin_(~isnan(mIin));

    % Predicting expression probabilty: discretize expression level first
    ax_I = linspace(-max(sIin_*3),max(mIin_+sIin_*3),1000);
    pIin_all = ax_I*0;
    for i = 1:numel(mIin_)
        pIin_all = pIin_all + normpdf(ax_I,mIin_(i),sIin_(i));
    end
    % Prepare the map
    pos_prediction_map = zeros(numel(mIin_));
    for i = 1:numel(mIin_)
        for j=1:numel(mIin_)
            tmp = ...
                    normpdf(ax_I,mIin_(i),sIin_(i)).*...
                    normpdf(ax_I,mIin_(j),sIin_(j))./...
                    pIin_all;
            pos_prediction_map(i,j) = nansum(tmp);
        end
        pos_prediction_map(i,:)=pos_prediction_map(i,:)/sum(pos_prediction_map(i,:));
    end
    HeatMap_(pos_prediction_map,pos_range_,pos_range_,[0 max(pos_prediction_map(:))]);
    set(gca,'Ydir','normal');
    xlabel('Real position x (%EL)');
    ylabel('Predicted position x'' (%EL)');
end