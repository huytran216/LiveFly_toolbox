function parent = HeatMap_(Map,ax,ay,cscale,parent)
    if ~exist('parent','var')
        parent = gca;
    end
    if ~exist('ax','var')
        [ax,ay]=size(Map');
        ax=1:ax;
        ay=1:ay;
    end
    if ~exist('cscale','var')
        cscale=[0 1];
    end
    surf(ax,ay,Map,'Parent',parent,'EdgeColor','none');
    xlim([min(ax) max(ax)]);
    ylim([min(ay) max(ay)]);
    colormap(hot);
    view([0 0 1]);
    set(parent,'Ydir','reverse');
    colorbar;
    caxis(cscale);