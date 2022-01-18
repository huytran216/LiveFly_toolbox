function newind = randperm_subset(tscnt)
    newind = 1:numel(tscnt);
    for i=unique(tscnt)
        tmp = find(tscnt==i);
        newind(tmp) = tmp(randperm(numel(tmp)));
    end
end