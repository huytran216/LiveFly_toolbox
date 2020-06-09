function [b] = barplot_bound_1x2x(x,lb,ub,label,baseline,orientation)
    plot_hor=0;
    if exist('orientation','var')
        if strcmp(orientation,'h')
            plot_hor=1;
        end
    end
    x= x(:);
    lb = lb(:);
    ub = ub(:);
    ngroups = numel(x)/2;
    nbars = 2;
    % Calculating the width for each bar group
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    if ~plot_hor
        x=reshape(x,[numel(x)/2 2]);
        lb=reshape(lb,[numel(x)/2 2]);
        ub=reshape(ub,[numel(x)/2 2]);
        
        b=bar(x,'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'BarWidth',0.8);
        hold on;
        for i = 1:nbars
            x_axis = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
            errorbar(x_axis, (lb(:,i)+ub(:,i))/2, abs(ub(:,i)-lb(:,i))/2, abs(ub(:,i)-lb(:,i))/2 ,'LineStyle','none','color','k','LineWidth',2);
        end
        if exist('label','var')
            set(gca,'XTick',1:numel(label(1:2:end)),'XTickLabel',label(1:end/2));
        end
        if exist('baseline','var')
            b(1).BaseValue = baseline;
            yl = get(gca,'ylim');
            ylim([baseline yl(2)]);
        end

    else
        x = x(end:-1:1);
        lb = lb(end:-1:1);
        ub = ub(end:-1:1);
        
        x=reshape(x,[numel(x)/2 2]);
        lb=reshape(lb,[numel(x)/2 2]);
        ub=reshape(ub,[numel(x)/2 2]);
        
        label = label(end:-1:1);
        b=barh(x,'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'BarWidth',0.8);
        hold on;
        for i = 1:nbars
            x_axis = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
            errorbar((lb(:,i)+ub(:,i))/2, x_axis, abs(ub(:,i)-lb(:,i))/2, abs(ub(:,i)-lb(:,i))/2 ,'horizontal','LineStyle','none','color','k','LineWidth',2);
        end
        
        if exist('label','var')
            set(gca,'YTick',1:numel(label(end/2+1:end)),'YTickLabel',label(end/2:-1:1));
        end
        if exist('baseline','var')
            b(1).BaseValue = baseline;
            xl = get(gca,'xlim');
            xlim([baseline xl(2)]);
        end
    end
    if numel(x)>2
        set(b(1),'FaceColor',corder(2));
        set(b(2),'FaceColor',corder(4));
    end
    