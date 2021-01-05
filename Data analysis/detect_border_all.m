function [xhat,hhat,what,vhat,CIxhat,CIhhat,CIwhat,CIvhat] = detect_border_all(x,y,limit,isplot,fitoption,erroroption,isfitt0)
    % Detect the border position based on the cell position (x), cell
    % feature (y).
    % If feature =10 (interphase duration), then Hill coeff = 10;
    % Input:
    %   x: 1xN embryos of cell position
    %   y: 1xN embryos of feature value
    %   limit: define the value of vborder*2
    %   isplot: plot the figure or not
    %   fitoption: options to fit the function
            % fitoption(1): share maximum level
            % fitoption(2): share Hill coefficient
            % fitoption(3): use relative time (not used)
    %   erroroption: calculate confidance level
            % erroroption(1): xhat
            % erroroption(2): yhat
            % erroroption(3): what
            % erroroption(4): vhat
    % Output: 
    %   xhat: position of the border
    %   hhat: hill coefficient
    %   what: width of the border - from 5% to 95% maximum expression
    %   vborder: value of intensity at the border (half maximum)
    N=numel(x);
    % Find the range of y, and calculate the middle point (border_value)
    
    if ~exist('limit','var')
        limit=0;
    end
    if ~exist('erroroption','var')
        erroroption=zeros(1,4);
    else
        if nargout<=4
            erroroption=zeros(1,4);
        end
    end
        
    % Fitting the border line
    function v=fitline(x,beta)
        % beta1: amplitude
        % beta2: border position
        % beta3: order/hill coeff
        v=sigmf(x,beta([3 2]))*beta(1)*2;
    end

    % Case010: limit==0, fitoption = [0 1 0];
        % Error function
        function fval=obj010(beta)
            % beta: size of 2xN+1, containing [(N)border_value (N)border_position (1)hill_function]
            fval=0;
            for k=1:N
                fval=fval+sum(abs(y{k}-fitline(x{k},beta([k k+N end]))).^2);
            end
        end
    
    % Case011: limit==0, fitoption = [0 1 1];
        % Error function
        function fval=obj011(beta)
            % beta: size of N+2, containing [(N)border_value (1)border_position (1)hill_function]
            fval=0;
            for k=1:N
                fval=fval+sum(abs(y{k}-fitline(x{k},beta([k end-1 end]))).^2);
            end
        end
    
    % Case100: limit==0, fitoption = [1 0 0];
        % Error function
        function fval=obj100(beta)
            % beta: size of 2xN+1, containing [(1)border_value (N)border_position (N)hill_function]
            fval=0;
            for k=1:N
                fval=fval+sum(abs(y{k}-fitline(x{k},beta([1 1+k 1+N+k]))).^2);
            end
        end
    
    % Case101: limit>0, fitoption = [1 0 1];
        % Error function
        function fval=obj101(beta)
            % beta: size of N+2, containing [(1)border_value (1)border_position (N)hill_function]
            fval=0;
            for k=1:N
                fval=fval+sum(abs(y{k}-fitline(x{k},beta([1 2 2+k]))).^2);
            end
        end
    
    % Case110: limit==0, fitoption = [1 1 0];
        % Error function
        function fval=obj110(beta)
            % beta: size of 2xN+1, containing [(1)border_value (N)border_position (1)hill_function]
            fval=0;
            for k=1:N
                fval=fval+sum(abs(y{k}-fitline(x{k},beta([1 1+k end]))).^2);
            end
        end
    % Case111: limit==0, fitoption = [1 1 1];
        % Error function
        function fval=obj111(beta)
            % beta: size of 2xN+1, containing [(1)border_value (1)border_position (1)hill_function]
            fval=0;
            for k=1:N
                fval=fval+sum(abs(y{k}-fitline(x{k},beta([1 2 3]))).^2);
            end
        end
    
    % Case000: limit==0, fitoption = [0 0 0];
        function fval=obj000(beta)
            % beta: size of 3xN, containing [(N)border_value (N)border_position (N)hill_function]
            fval=0;
            for k=1:N
                fval=fval+sum(abs(y{k}-fitline(x{k},beta([k k+N k+2*N]))).^2);
            end
        end
    
    % Case001: limit==0, fitoption = [0 0 1];
        function fval=obj001(beta)
            % beta: size of 3xN, containing [(N)border_value (1)border_position (N)hill_function]
            fval=0;
            for k=1:N
                fval=fval+sum(abs(y{k}-fitline(x{k},beta([k k+N k+N+1]))).^2);
            end
        end
    
    %% Create initial value
    beta0=[];
    beta_lb=[];
    beta_ub=[];
    betaidx=[];
    if ~fitoption(1)
        betaidx=[betaidx numel(beta0)+(1:N)];
        beta0=[beta0 1*ones(1,N)];
        beta_ub=[beta_ub 1e10*ones(1,N)];
        beta_lb=[beta_lb 1e-10*ones(1,N)];
    else
        betaidx=[betaidx numel(beta0)+ones(1,N)*1];
        beta0=[beta0 1];
        beta_ub=[beta_ub 1e10];
        beta_lb=[beta_lb 1e-10];
    end
    if ~fitoption(3)
        betaidx=[betaidx numel(beta0)+(1:N)];
        beta0=[beta0 -5*ones(1,N)];
        beta_ub=[beta_ub 50*ones(1,N)];
        beta_lb=[beta_lb -50*ones(1,N)];
    else
        betaidx=[betaidx numel(beta0)+ones(1,N)*1];
        beta0=[beta0 -5];
        beta_ub=[beta_ub 50];
        beta_lb=[beta_lb -50];
    end
    if ~isfitt0
        H_ = [-5 -20 -0.1];
    else
        H_ = [-20 -21 -19];
    end
    if ~fitoption(2)
        betaidx=[betaidx numel(beta0)+(1:N)];
        beta0=[beta0 H_(1)*ones(1,N)];
        beta_ub=[beta_ub H_(3)*ones(1,N)];
        beta_lb=[beta_lb H_(2)*ones(1,N)];
    else
        betaidx=[betaidx numel(beta0)+ones(1,N)*1];
        beta0=[beta0 H_(1)];
        beta_ub=[beta_ub H_(3)];
        beta_lb=[beta_lb H_(2)];
    end
    %% Get function
    fun=@mean;
    setfunname=['fun=@obj' num2str(fitoption(1)) num2str(fitoption(2)) num2str(fitoption(3)) ';'];
    eval(setfunname);
    %% Eval function:
    [betahat,fhat]=fminsearchbnd(fun,beta0,beta_lb,beta_ub);
    %% Extract from betahat:
    vhat=betahat(betaidx(1:N));
    xhat=betahat(betaidx(N+(1:N)));
    hhat=-betahat(betaidx(2*N+(1:N)))*20;
    what=log(19)./hhat*20;
    %% Calculate the errorbar
    p=0.05;
    nsample=numel(cell2mat(x));
    % Calculate tolerance objective
    fval=fhat/nsample;
    diffloglikelihood = chi2inv(1-p,1);
    sigma=sqrt(fval);
    floglike = -nsample*fval/sigma^2;
    
    % Error bar for xhat
    CIxhat=[0 0];
    if erroroption(1)        
        if fitoption(3)
            xhat_range=[0.1:0.1:20];
            lb=0;
            ub=0;
            for xhat_=1:numel(xhat_range)
                % xhat fixed, hhat, vhat fixed
                if ~ub
                    if -nsample*obj111([vhat(1) xhat(1)+xhat_range(xhat_) -hhat(1)/20])/sigma^2/nsample<(floglike-diffloglikelihood)
                        ub=1;
                        CIxhat(2)=xhat(1)+xhat_range(xhat_);
                    end
                end
                if ~lb
                    if -nsample*obj111([vhat(1) xhat(1)-xhat_range(xhat_) -hhat(1)/20])/sigma^2/nsample<(floglike-diffloglikelihood)
                        lb=1;
                        CIxhat(1)=xhat(1)-xhat_range(xhat_);
                    end
                end                    
            end            
        end 
    end
    
    % Error bar for hhat
    CIhhat=[1 1];
    if erroroption(2)
        if fitoption(2)
            hhat_range=[0.1:0.1:20];
            lb=0;
            ub=0;
            for hhat_=1:numel(hhat_range)
                % xhat fixed, hhat, vhat fixed
                if ~ub
                    if -nsample*obj111([vhat(1) xhat(1) -(hhat(1)+hhat_range(hhat_))/20])/sigma^2/nsample<(floglike-diffloglikelihood)
                        ub=1;
                        CIhhat(2)=hhat(1)+hhat_range(hhat_);
                    end
                end
                if ~lb
                    if -nsample*obj111([vhat(1) xhat(1) -(hhat(1)-hhat_range(hhat_))/20])/sigma^2/nsample<(floglike-diffloglikelihood)
                        lb=1;
                        CIhhat(1)=hhat(1)-hhat_range(hhat_);
                    end
                end                    
            end            
        end
    end
    % Error bar for what
    if erroroption(3)
        CIwhat=log(19)./CIhhat([2 1])*20;
    end
    % Error bar for vhat
    CIvhat=[0 0];
    if erroroption(4)
        if fitoption(1)
            vhat_range=10.^[0.01:0.01:0.5];
            lb=0;
            ub=0;
            for vhat_=1:numel(vhat_range)
                % vhat not fixed, hhat, xhat fixed
                if ~ub
                    if -nsample*obj111([vhat(1)*vhat_range(vhat_) xhat(1) -hhat(1)/20])/sigma^2/nsample<(floglike-diffloglikelihood)
                        ub=1;
                        CIvhat(2)=vhat(1)*vhat_range(vhat_);
                    end
                end
                if ~lb
                    if -nsample*obj111([vhat(1)/vhat_range(vhat_) xhat(1) -hhat(1)/20])/sigma^2/nsample<(floglike-diffloglikelihood)
                        lb=1;
                        CIvhat(1)=vhat(1)/vhat_range(vhat_);
                    end
                end
            end
        end
    end
end