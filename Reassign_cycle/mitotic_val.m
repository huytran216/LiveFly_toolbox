function  y=mitotic_val(p,x)
    % p: function parameters
        % p(1): max position for y
        % p(2): x that maximize y
        % p(3): speed that y is decreasing
    % x: function input
        y=p(1)-abs(p(2)-x)*p(3);
end