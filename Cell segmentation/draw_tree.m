function draw_tree(tree_parent,tree_time,cc)
    % Draw the tree on figure h
    % get tree layout
    % Draw the cell cycle if any
    [x,y]=treelayout(tree_parent);
    z=1:numel(tree_parent);
    y=tree_time;
    hold on;
    for i=1:numel(x)
        if tree_parent(i)
            plot3([x(i) x(tree_parent(i))],[y(i) y(tree_parent(i))],[z(i) z(tree_parent(i))]);
        end
    end
    plot3(x,y,z,'o');
    set (gca,'Ydir','reverse');
    ylabel('<< frame <<');
    title('Lineage tree');
    dcm_obj = datacursormode(figure(gcf));
    set(dcm_obj,'UpdateFcn',@myupdatefcn);
    % Draw the nuclear cycle:
    [cctype,~,pos]=unique(cc);
    newcolor=cell(numel(cctype),1);
    cchue=0;
        % Get a specific color
        for i=1:numel(cctype)
            cchue=(cchue+0.87);
            cchue=cchue-floor(cchue);
            newcolor{i}=hsv2rgb([cchue 0.5 0.5]);
            % Add label on the nuclear cycle
            if cctype(i)>0
                draw_pt=find(cc==cctype(i),1,'first');
                text(1,draw_pt+5,['nc' num2str(cctype(i))]);
            end
        end
        % Draw the color strip
        for i=1:numel(cc)
            plot([1.1 1.2],[i i],'color',newcolor{pos(i)},'LineWidth',3);
        end
end