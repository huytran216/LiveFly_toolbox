function []=plot_marker(fea1,fea2,tscnt,color)
   [C]=unique(tscnt);
   markerlist='+o*.xsd^v<>ph';
   
   for i=1:numel(C)
       tmp=tscnt==C(i);
       plot(fea1(tmp),fea2(tmp),[markerlist(i) color]);hold on;
   end
   hold off;