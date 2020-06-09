function Iout = medfilt2_roi(I,mn,mask)

% Find true pixel in mask
[sx,sy]=find(mask);
m=mn(1);
n=mn(2);
mn=size(I);

% mask
mask=conv2(double(mask),ones(m,n),'same');
I(~mask)=-1;

%fun=@(x) median(x.data(:));
Iout=nlfilter(I,[m n],@fun);

function fval = fun(x)
        if any(x==-1)
            fval=1e10;  % Very high background - remove spots here
        else
            fval=median(x(:));
        end
    end
end
