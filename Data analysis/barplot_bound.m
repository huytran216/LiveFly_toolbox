function [] = barplot_bound(x,lb,ub,label)
    bar(x,'FaceColor',[1 1 1],'EdgeColor',[0 0 0]);
    hold on;
    errorbar(1:numel(x),(lb+ub)/2,abs(ub-lb)/2,abs(ub-lb)/2,'LineStyle','none','color','k','LineWidth',2);
    if exist('label','var')
        set(gca,'XTick',1:numel(label),'XTickLabel',label);
    end