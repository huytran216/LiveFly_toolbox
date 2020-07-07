function [c0_grid] = Bcd_gradient(xax,gridt,D,L)
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
    % set initial condition:
    cax = gridx*0; cax(1) = 1;
    %cax_rec = zeros(numel(gridt),numel(gridx));
    dcax = gridx*0;
    d2cax = gridx*0;
    % Predicted lambda:
    lambda = sqrt(D*L);
    semilogy(xax,exp(-xax/lambda),'--');
    hold on;
    
    for tidx = 2:numel(gridt)
        % Calculate decay:
        for xidx = 1:numel(xax)
            dcax(xidx) = -cax(xidx)/L;
        end
        % Birth only at anterior pole
        dcax(1) = dcax(1)+0.1/dt;
        % Calculate 2nd derivatives over x:
        for xidx = 2:numel(xax)-1
            d2cax(xidx) = D*(cax(xidx-1)+cax(xidx+1)-2*cax(xidx))/2/dx.^2;
        end
        d2cax(1) = D*(+cax(2)-cax(1))/dx.^2;
        d2cax(end) = D*(-cax(end)+cax(end-1))/dx.^2;
        % Final: 
        cax = cax + dt*(dcax+d2cax);
        % Record:
        %cax_rec = [cax_rec;cax];
        if mod(tidx,1000)==2
            semilogy(gridx,cax+1e-5);
            hold on;
        end
    end
    
    % Final: 
    c0_grid = cax;