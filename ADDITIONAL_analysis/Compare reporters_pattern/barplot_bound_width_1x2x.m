function [b] = barplot_bound(x,x_lb,x_ub,y,y_lb,y_ub,label,baseline,orientation)
    plot_hor=0;
    if exist('orientation','var')
        if strcmp(orientation,'h')
            plot_hor=1;
        end
    end
    x= x(:);
    x_lb = x_lb(:);
    x_ub = x_ub(:);
    ngroups = numel(x)/2;
    nbars = 2;
    % Calculating the width for each bar group
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    if ~plot_hor
        x=reshape(x,[numel(x)/2 2]);
        x_lb=reshape(x_lb,[numel(x)/2 2]);
        x_ub=reshape(x_ub,[numel(x)/2 2]);
        
        b=bar(x,'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'BarWidth',0.8);
        hold on;
        for i = 1:nbars
            x_axis = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
            errorbar(x_axis, (x_lb(:,i)+x_ub(:,i))/2, abs(x_ub(:,i)-x_lb(:,i))/2, abs(x_ub(:,i)-x_lb(:,i))/2 ,'LineStyle','none','color','k','LineWidth',2);
        end
        if exist('label','var')
            set(gca,'XTick',1:numel(label(1:2:end)),'XTickLabel',label(1:end/2));
        end
        if exist('baseline','var')
            b(1).BaseValue = baseline(1);
            yl = get(gca,'ylim');
            ylim([baseline(1) yl(2)]);
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
        
        x=reshape(x,[numel(x)/2 2]);
        x_lb=reshape(x_lb,[numel(x)/2 2]);
        x_ub=reshape(x_ub,[numel(x)/2 2]);
        x1=reshape(x1,[numel(x)/2 2]);
        x1_lb=reshape(x1_lb,[numel(x)/2 2]);
        x1_ub=reshape(x1_ub,[numel(x)/2 2]);
        x2=reshape(x2,[numel(x)/2 2]);
        x2_lb=reshape(x2_lb,[numel(x)/2 2]);
        x2_ub=reshape(x2_ub,[numel(x)/2 2]);
        
        label = label(end:-1:1);
        b=barh(x2,'FaceColor',[.7 .7 .7],'EdgeColor',[0 0 0],'BarWidth',0.8);hold on;
        c=barh(x,'FaceColor',[.7 .7 .7],'EdgeColor',[0 0 0],'BarWidth',0.8);hold on;
        d=barh(x1,'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'BarWidth',0.8);hold on;
        hold on;
        for i = 1:nbars
            x_axis = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
            errorbar((x_lb(:,i)+x_ub(:,i))/2, x_axis, abs(x_ub(:,i)-x_lb(:,i))/2, abs(x_ub(:,i)-x_lb(:,i))/2 ,'horizontal','LineStyle','none','color','k','LineWidth',2);
            errorbar((x1_lb(:,i)+x1_ub(:,i))/2, x_axis, abs(x1_ub(:,i)-x1_lb(:,i))/2, abs(x1_ub(:,i)-x1_lb(:,i))/2 ,'horizontal','LineStyle','none','color','k','LineWidth',2);
            errorbar((x2_lb(:,i)+x2_ub(:,i))/2, x_axis, abs(x2_ub(:,i)-x2_lb(:,i))/2, abs(x2_ub(:,i)-x2_lb(:,i))/2 ,'horizontal','LineStyle','none','color','k','LineWidth',2);
        end
        
        if exist('label','var')
            set(gca,'YTick',1:numel(label(end/2+1:end)),'YTickLabel',label(end/2:-1:1));
        end
        if exist('baseline','var')
            b(1).BaseValue = baseline(1);
            xl = get(gca,'xlim');
            xlim([baseline(1) xl(2)]);
        end
    end
    if numel(x)>2
        set(b(1),'FaceColor',corder(2));
        set(b(2),'FaceColor',corder(4));
        set(c(1),'FaceColor',corder(2));
        set(c(2),'FaceColor',corder(4));
    end
    