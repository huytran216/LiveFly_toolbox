
i=1;
pos_range = -85:85;
lambda = 20;
lambdaH = 5;

mI_rec{1,1} = 1./(1+exp(pos_range/lambdaH));
sI_rec{1,1} = pos_range*0+0.01;
pos_rec{1,1} = pos_range;
mI_rec{1,2} = 1./(1+exp((pos_range+20)/lambdaH));
sI_rec{1,2} = pos_range*0+0.01;
pos_rec{1,2} = pos_range;
figure(1);
subplot(121);
color2x = corder(2);
color1x = corder(4);
shadedErrorBar(pos_range,mI_rec{1,1},sI_rec{1,1},{'color',color2x,'Display','2x','LineWidth',2},0.8,2);hold on;
shadedErrorBar(pos_range,mI_rec{1,2},sI_rec{1,2},{'color',color1x,'Display','1x','LineWidth',2},0.8,2);
xlim([-45 45]);
ylim([0 1.1]);
set(gca,'YTick',[],'XTick',[]);
xlabel('Nuclei position');
ylabel('P_{Spot}');
%% Calculate shift in boundary
pos_range1 = pos_rec{i,1}>=-85;
pos_range2 = pos_rec{i,2}>=-85;
    figure(2);
    [diff_plot{i},ax,ay,diff_grid{i},prange1,prange2]=shift_prediction_map(pos_rec{i,1}(pos_range1),mI_rec{i,1}(pos_range1),sI_rec{i,1}(pos_range1),...
        pos_rec{i,2}(pos_range2),mI_rec{i,2}(pos_range2),sI_rec{i,2}(pos_range2));
    
    %subplot(121);
    %xlim([-30 20]);
    %ylim([-30 20]);
    
    %subplot(122);
    
    colormap(flipud(gray));
    if i==1
        diff_all = diff_plot{1};
        diff_grid_all = diff_grid{1};
    else
        diff_all = diff_all*diff_plot{i};
        diff_grid_all = diff_grid_all.*diff_grid{i};
    end
    figure(1);
    subplot(122);
        HeatMap_(log(diff_grid_all'+1e-10)-5,prange1,prange2,[log(1e-10) log(max(diff_grid_all(:)))]);
        hold on;
        plot3(prange1(:),prange1(:)*0,prange1(:)*0+log(max(diff_grid_all(:))),'LineStyle','-','color','k','LineWidth',1);
        maxlog = log(max(diff_grid_all(:)));
        set(gca,'Ydir','normal');
        h= colormap(flipud(gray));
        xlabel('Original nuclei position (%EL)');
        ylabel('Predicted shift (%EL)');
        xlim([-45 45]);
        ylim([-45 45]);
        zlim([log(1e-10) maxlog]);
        set(gca,'XTick',[],'YTick',[]);