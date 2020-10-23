function x=randbetween(a,b)
    % assuming size(a)==size(b)
    x = a+(b-a).*rand(size(a));