alpha_range = 0:0.05:1;
rT_range = 10.^(-0.5:0.05:0.5);
rTeff = zeros(numel(alpha_range),numel(rT_range));
rerr = zeros(numel(alpha_range),numel(rT_range));

for i=1:numel(alpha_range)
    for j=1:numel(rT_range)
        x = alpha_range(i);
        y = rT_range(j);
        rTeff(i,j) = 1/(x^2/y + (1-x)^2);
        rerr(i,j)=(x^2/y + (1-x)^2);
    end
end
contour(log10(rT_range),alpha_range,rerr,[0:0.1:4],'ShowText','off','LineWidth',1);hold on;
contour(log10(rT_range),alpha_range,rerr,[0:0.2:1 1.5:0.5:4],'ShowText','on','LineWidth',3);
%set(gca,'XScale','log');
xlabel('(\partial Y_{real}/\partial Y)^2');
ylabel('\alpha')

hold on;
plot3(alpha_range*0+log10(750/1100),alpha_range,alpha_range*0+10,'LineStyle','--','LineWidth',2,'color','k');
plot3(alpha_range*0+log10(600/750),alpha_range,alpha_range*0+10,'LineStyle','--','LineWidth',2,'color','k');
plot3(alpha_range*0+log10(1100/2700),alpha_range,alpha_range*0+10,'LineStyle','--','LineWidth',2,'color','k');

rT_real = 1100/2700;

rTeff_best = (0.2.^2/rT_real + (1-0.2).^2);
rTeff_best