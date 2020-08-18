function [b] = barplot_bound(x,lb,ub,label,baseline,orientation)
    plot_hor=0;
    if exist('orientation','var')
        if strcmp(orientation,'h')
            plot_hor=1;
        end
    end
    if ~plot_hor
        b=bar(x,'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'BarWidth',0.8);
        hold on;
        errorbar(1:numel(x),(lb+ub)/2,abs(ub-lb)/2,abs(ub-lb)/2,'LineStyle','none','color','k','LineWidth',2);
        if exist('label','var')
            set(gca,'XTick',1:numel(label),'XTickLabel',label);
        end
        if exist('baseline','var')
            b(1).BaseValue = baseline(1);
            yl = get(gca,'ylim');
            ylim([baseline(1) yl(2)]);
        end

    else
        x = x(end:-1:1);
        lb = lb(end:-1:1);
        ub = ub(end:-1:1);
        label = label(end:-1:1);
        b=barh(x,'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'BarWidth',0.8);
        hold on;
        errorbar((lb+ub)/2,1:numel(x),abs(ub-lb)/2,abs(ub-lb)/2,'horizontal','LineStyle','none','color','k','LineWidth',2);
        if exist('label','var')
            set(gca,'YTick',1:numel(label),'YTickLabel',label);
        end
        if exist('baseline','var')
            b(1).BaseValue = baseline(1);
            xl = get(gca,'xlim');
            xlim([baseline(1) xl(2)]);
        end
    end
    
    