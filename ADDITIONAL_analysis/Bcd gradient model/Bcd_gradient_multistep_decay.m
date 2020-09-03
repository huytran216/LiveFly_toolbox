function [c_grid,gridx] = Bcd_gradient_multistep_decay(xax,gridt,D,L)
% xax: position axes (in %EL)
% tax: time axes (in s)
% D: diffusion coefficient (unit in micron^2/s)
% L: Protein life time (s)
    % convert x %EL to micron
    xax = (xax+50)*5;
    % Define function:
    syms x;
    syms t;
    
    % Grid
    dx = 0.5;
    dt =0.01;
    gridx = 0:dx:500;
    gridt = 0:dt:10000;
    % Predicted lambda:
    lambda = sqrt(D*L);
    % set initial condition:
    cax_1 = L*exp(-gridx/lambda);
    %cax_1 = gridx*0;
    
    dcax_1 = gridx*0;
    d2cax_1 = gridx*0;
    
    %cax_2 = exp(-gridx/lambda);
    stability = 3;
    cax_2 = exp(-gridx/lambda/stability);
    dcax_2 = gridx*0;
    d2cax_2 = gridx*0;
    
    figure;
    semilogy(xax,exp(-xax/lambda),'LineStyle','--','LineWidth',2,'color','k');
    hold on;
    for tidx = 2:numel(gridt)
        % Calculate decay:
        dcax_1 = -cax_1/L;
        dcax_2 = +cax_1/L - cax_2/L/stability;
        % Birth only at anterior pole
        dcax_1(1) = dcax_1(1) +1;
        % Diffusion term
        d2cax_1(2:end-1) = diff(diff(cax_1));
        d2cax_1(1) = (+cax_1(2)-cax_1(1));
        d2cax_1(end) = (-cax_1(end)+cax_1(end-1));
        d2cax_2(2:end-1) = diff(diff(cax_2));
        d2cax_2(1) = (+cax_2(2)-cax_2(1));
        d2cax_2(end) = (-cax_2(end)+cax_2(end-1));
        % Final: 
        cax_1 = cax_1 + dt*(dcax_1+D*d2cax_1/dx^2);
        cax_2 = cax_2 + dt*(dcax_2+D*d2cax_2/dx^2);
        % Record:
        %cax_rec = [cax_rec;cax];
        if mod(tidx*dt,100)==0
            derivative12 = [tidx*dt sum(abs(dcax_2+D*d2cax_2/dx^2)) sum(abs(dcax_1+D*d2cax_1/dx^2))]
            semilogy(gridx,cax_2+1e-10); hold on;
            %semilogy(gridx,cax_1+1e-10); hold on;
            if sum(derivative12(2:3))<1e-4
                break;
            end
        end
    end
    semilogy(xax,exp(-xax/lambda),'LineStyle','--','LineWidth',2,'color','k');
    hold on;
    % Final: 
    figure;
    semilogy(gridx,cax_1+1e-10,'Display','c1 at SS'); hold on;
    semilogy(gridx,cax_2+1e-10,'Display','c2 at SS'); hold on;
    semilogy(gridx,5*cax_1+cax_2+1e-10,'Display','c1+c2 at SS'); hold on;
    semilogy(xax,exp(-xax/lambda),'LineStyle','--','LineWidth',2,'color','k','Display','c1 prediction');
    xlabel('position');
    ylabel('c');
    legend show
    
    c_grid = cax_2;