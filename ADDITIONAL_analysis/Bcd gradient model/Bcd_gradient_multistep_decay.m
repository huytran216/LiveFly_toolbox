function [c0_grid] = Bcd_gradient_multistep_decay(xax,gridt,D,L)
% xax: position axes (in %EL)
% tax: time axes (in s)
% D: diffusion coefficient (unit in micron^2/s)
% L: Protein life time 
    % convert x %EL to micron
    xax = (xax+50)*5;
    % Define function:
    syms x;
    syms t;
    
    % Grid
    dx = 0.1;
    dt =0.001;
    gridx = 0:dx:500;
    gridt = 0:dt:1000;
    % Predicted lambda:
    lambda = sqrt(D*L);
    % set initial condition:
    cax_1 = exp(-gridx/lambda);
    dcax_1 = gridx*0;
    d2cax_1 = gridx*0;
    
    cax_2 = exp(-gridx/lambda);
    dcax_2 = gridx*0;
    d2cax_2 = gridx*0;
    
    
    semilogy(xax,exp(-xax/lambda),'--');
    hold on;
    
    for tidx = 2:numel(gridt)
        % Calculate decay:
        for xidx = 1:numel(gridx)
            dcax_1(xidx) = -cax_1(xidx)/L;
            dcax_2(xidx) = +cax_1(xidx)/L;
        end
        % Birth only at anterior pole
        dcax_1(1) = dcax_1(1)+1;
        % Calculate 2nd derivatives over x:
        for xidx = 2:numel(gridx)-1
            d2cax_1(xidx) = D*(cax_1(xidx-1)+cax_1(xidx+1)-2*cax_1(xidx))/2/dx.^2;
            d2cax_2(xidx) = D*(cax_2(xidx-1)+cax_2(xidx+1)-2*cax_2(xidx))/2/dx.^2;
        end
        d2cax_1(1) = D*(+cax_1(2)-cax_1(1))/dx.^2;
        d2cax_1(end) = D*(-cax_1(end)+cax_1(end-1))/dx.^2;
        d2cax_2(1) = D*(+cax_2(2)-cax_2(1))/dx.^2;
        d2cax_2(end) = D*(-cax_2(end)+cax_2(end-1))/dx.^2;
        % Final: 
        cax_1 = cax_1 + dt*(dcax_1+d2cax_1);
        cax_2 = cax_2 + dt*(dcax_2+d2cax_2);
        % Record:
        %cax_rec = [cax_rec;cax];
        if mod(tidx,10000)==2
            semilogy(gridx,cax_2+1e-10);
            hold on;
        end
    end
    
    % Final: 
    c0_grid = cax_1;