function show_info(x)
    for i=1:numel(x)/4
        kon = x(i*2-1);
        koff = x(i*2);
        display(num2str(i));
        display(['   PON: ' num2str(kon/(kon+koff))]);
        display(['   kon: ' num2str(kon)]);
        display(['   koff: ' num2str(koff)]);
        display(['   tau: ' num2str(1/(kon+koff))]);
        display(['  MaxI: ' num2str(x(end-numel(x)/2+i))]);
        display([' Delay: ' num2str(x(end-numel(x)/4+i)) ' s']);
    end
end