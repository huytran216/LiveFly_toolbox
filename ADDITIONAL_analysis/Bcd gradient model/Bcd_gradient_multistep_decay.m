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
    gridt = 0:dt:10000;
    % Predicted lambda:
    lambda = sqrt(D*L);
    % set initial condition:
    cax_1 = exp(-gridx/lambda);
    %cax_1 = gridx*0;
    dcax_1 = gridx*0;
    d2cax_1 = gridx*0;
    
    %cax_2 = exp(-gridx/lambda);
    cax_2 = gridx*0;
    dcax_2 = gridx*0;
    d2cax_2 = gridx*0;
    
    figure;
    semilogy(xax,dt*L*exp(-xax/lambda),'LineStyle','--','LineWidth',2,'color','k');
    hold on;
    for tidx = 2:numel(gridt)
        % Calculate decay:
        dcax_1 = -cax_1/L;
        dcax_2 = +cax_1/L - cax_2/L/2;
        % Birth only at anterior pole
        dcax_1(1) = dcax_1(1) +1;
        % Diffusion term
        d2cax_1(2:end-1) = diff(diff(cax_1));
        d2cax_2(2:end-1) = diff(diff(cax_2));
        d2cax_1(1) = D*(+cax_1(2)-cax_1(1));
        d2cax_1(end) = D*(-cax_1(end)+cax_1(end-1));
        d2cax_2(1) = D*(+cax_2(2)-cax_2(1));
        d2cax_2(end) = D*(-cax_2(end)+cax_2(end-1));
        % Final: 
        cax_1 = cax_1 + dt*(dcax_1+d2cax_1/dx^2);
        cax_2 = cax_2 + dt*(dcax_2+d2cax_2/dx^2);
        % Record:
        %cax_rec = [cax_rec;cax];
        if mod(tidx*dt,100)==0
            %semilogy(gridx,cax_2+1e-10); hold on;
            semilogy(gridx,cax_1+1e-10); hold on;
        end
    end
    semilogy(xax,dt*L*exp(-xax/lambda),'LineStyle','--','LineWidth',2,'color','k');
    hold on;
    % Final: 
    c0_grid = cax_1;