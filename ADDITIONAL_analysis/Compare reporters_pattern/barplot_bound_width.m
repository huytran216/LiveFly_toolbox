function [b] = barplot_bound_width(x,x_lb,x_ub,y,y_lb,y_ub,label,baseline,orientation)
    % x is the boundary position
    % y is the boudary width
    plot_hor=0;
    if exist('orientation','var')
        if strcmp(orientation,'h')
            plot_hor=1;
        end
    end
    if ~plot_hor
        b=bar(x,'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'BarWidth',0.8);
        hold on;
        errorbar(1:numel(x),(x_lb+x_ub)/2,abs(x_ub-x_lb)/2,abs(x_ub-x_lb)/2,'LineStyle','none','color','k','LineWidth',2);
        if exist('label','var')
            set(gca,'XTick',1:numel(label),'XTickLabel',label);
        end
        if exist('baseline','var')
            b(1).BaseValue = baseline(1);
            ylim(baseline);
        end
        
    else
        x = x(end:-1:1);
        x_lb = x_lb(end:-1:1);
        x_ub = x_ub(end:-1:1);
        y = y(end:-1:1);
        y_lb = y_lb(end:-1:1);
        y_ub = y_ub(end:-1:1);        
        label = label(end:-1:1);
        x1 = x-y/2;
        x1_lb = x-y_ub/2;
        x1_ub = x-y_lb/2;
        x2 = x+y/2;
        x2_ub = x+y_ub/2;
        x2_lb = x+y_lb/2;
        b=barh(x2,'FaceColor',[.7 .7 .7],'EdgeColor',[0 0 0],'BarWidth',0.8);hold on;
        c=barh(x,'FaceColor',[.7 .7 .7],'EdgeColor',[0 0 0],'BarWidth',0.8);hold on;
        d=barh(x1,'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'BarWidth',0.8);
        hold on;
        errorbar((x_lb+x_ub)/2,1:numel(x),abs(x_ub-x_lb)/2,abs(x_ub-x_lb)/2,'horizontal','LineStyle','none','color','k','LineWidth',2);
        errorbar((x1_lb+x1_ub)/2,1:numel(x),abs(x1_ub-x1_lb)/2,abs(x1_ub-x1_lb)/2,'horizontal','LineStyle','none','color','k','LineWidth',1);
        errorbar((x2_lb+x2_ub)/2,1:numel(x),abs(x2_ub-x2_lb)/2,abs(x2_ub-x2_lb)/2,'horizontal','LineStyle','none','color','k','LineWidth',1);
        if exist('label','var')
            set(gca,'YTick',1:numel(label),'YTickLabel',label);
        end
        if exist('baseline','var')
            b(1).BaseValue = baseline(1);
            xl = get(gca,'xlim');
            xlim(baseline);
        end
        h=gca; h.YAxis.TickLength = [0 0];
    end
    
    