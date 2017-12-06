function [xhat,hhat,what,vhat] = detect_border_all(x,y,limit,isplot,fitoption)
    % Detect the border position based on the cell position (x), cell
    % feature (y).
    % Input:
    %   x: 1xN embryos of cell position
    %   y: 1xN embryos of feature value
    %   limit: define the value of vborder*2
    %   isplot: plot the figure or not
    %   fitoption: options to fit the function
            % fitoption(1): share maximum level
            % fitoption(2): share Hill coefficient
            % fitoption(3): use relative time (not used)
    % Output: 
    %   xhat: position of the border
    %   hhat: hill coefficient
    %   what: width of the border
    %   vborder: value of intensity at the border
    N=numel(x);
    % Find the range of y, and calculate the middle point (border_value)
    
    if ~exist('limit','var')
        limit=0;
    end
    % Case1: limit>0, fitoption = [1 1];
        % Fitting the border line
        function v=fitline1(x,beta)
            v=sigmf(x,beta([2 1]))*limit;
        end
        % Error functions
        function fval=obj1(beta)
            % beta: size of N+1, containing [(N)border_position (1)hill_function]
            fval=0;
            for k=1:N
                fval=fval+sum(abs(y{k}-fitline1(x{k},beta([k end]))).^2);
            end
        end
    % Case2: limit==0, fitoption = [0 1];
        % Fitting the border line
        function v=fitline2(x,beta)
            v=sigmf(x,beta([3 2]))*beta(1)*2;
        end
        % Error function
        function fval=obj2(beta)
            % beta: size of 2xN+1, containing [(N)border_value (N)border_position (1)hill_function]
            fval=0;
            for k=1:N
                fval=fval+sum(abs(y{k}-fitline2(x{k},beta([k k+N end]))).^2);
            end
        end
    % Case3: limit==0, fitoption = [1 0];
        % Fitting the border line
        function v=fitline3(x,beta)
            v=sigmf(x,beta([3 2]))*beta(1)*2;
        end
        % Error function
        function fval=obj3(beta)
            % beta: size of 2xN+1, containing [(1)border_value (N)border_position (N)hill_function]
            fval=0;
            for k=1:N
                fval=fval+sum(abs(y{k}-fitline3(x{k},beta([1 1+k 1+N+k]))).^2);
            end
        end
    
    % Case4: limit==0, fitoption = [1 1];
        % Fitting the border line
        function v=fitline4(x,beta)
            v=sigmf(x,beta([3 2]))*beta(1)*2;
        end
        % Error function
        function fval=obj4(beta)
            % beta: size of 2xN+1, containing [(1)border_value (N)border_position (1)hill_function]
            fval=0;
            for k=1:N
                fval=fval+sum(abs(y{k}-fitline4(x{k},beta([1 1+k end]))).^2);
            end
        end
    % Case5: limit==0, fitoption = [0 0]; (fit embryo separately)
        % Fitting the border line
        function v=fitline5(x,beta)
            v=sigmf(x,beta([3 2]))*beta(1)*2;
        end
        % Error function
        function fval=obj5(beta)
            % beta: size of 2xN+1, containing [(1)border_value (1)border_position (1)hill_function]
            fval=0;
            fval=fval+sum(abs(y{cnt}-fitline5(x{cnt},beta([1 2 3]))).^2);
        end
    
    
    if limit
        % If vborder is known
        vborder_value=ones(1,N)*limit/2;

        % Perform the fitting
        beta0=[-5+zeros(1,N) -1];

        problem=opti('fun',@obj1,'ydata',0,'x0',beta0);
        betahat=solve(problem);
        betahat=[vborder_value(:);betahat(:)];
    else
        % If vborder is not known - make some educated guess
        for i=1:N
            [ncnt,ax]=hist(y{i}(x{i}<-10),10);
            [~,pos]=max(ncnt(2:end));
            vborder_value(i)=ax(pos+1)/2;
        end
        % Case 2. [N N 1]
        if (fitoption(1)==0)&&(fitoption(2)==1)
            % Perform the fitting
            beta0=[vborder_value -5+zeros(1,N) -1];
            problem=opti('fun',@obj2,'ydata',0,'x0',beta0);
            betahat=solve(problem);
            % Converting to the output    
            vhat=betahat(1:N);
            xhat=betahat(N+(1:N));
            hhat=-betahat(end)*20*ones(1,N);
            what=-1/2/betahat(end)*4*ones(1,N);
        end
        % Case 3. [1 N N]
        if (fitoption(1)==1)&&(fitoption(2)==0)
            % Perform the fitting
            beta0=[mean(vborder_value) -5+zeros(1,N) -ones(1,N)];
            problem=opti('fun',@obj3,'ydata',0,'x0',beta0);
            betahat=solve(problem);
            % Converting to the output    
            vhat=betahat(1)*ones(1,N);
            xhat=betahat(1+(1:N));
            hhat=-betahat(2+N:end)*20;
            what=-1/2/betahat(2+N:end)*4;
        end
        % Case 4. [1 N 1]
        if (fitoption(1)==1)&&(fitoption(2)==1)
            % Perform the fitting
            beta0=[mean(vborder_value) -5+zeros(1,N) -1];
            problem=opti('fun',@obj4,'ydata',0,'x0',beta0);
            betahat=solve(problem);
            % Converting to the output
            vhat=betahat(1)*ones(1,N);
            xhat=betahat(1+(1:N));
            hhat=-betahat(end)*20;
            what=-1/2/betahat(end)*4;
        end
        % Case 5. [1 1 1]xN
        if (fitoption(1)==0)&&(fitoption(2)==0)
            % Perform fitting separately
            for cnt=1:N
                beta0=[vborder_value(cnt) -5 -1];
                problem=opti('fun',@obj5,'ydata',0,'x0',beta0);
                betahat=solve(problem);
                % Converting to the output
                vhat(cnt)=betahat(1);
                xhat(cnt)=betahat(2);
                hhat(cnt)=-betahat(end)*20;
                what(cnt)=-1/2/betahat(end)*4;
            end
        end
    end
end