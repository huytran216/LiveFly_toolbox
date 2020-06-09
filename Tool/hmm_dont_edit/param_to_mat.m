function [A,R] = param_to_mat(varargin)
% Generate transition matrix based on the model and parameters
% INPUT:
%   vargin
    % 1. model: model type
    % 2. DT: sampling time
    % 3. following parameter vector for the model
% OUTPUT:
    % if nargin==2: generate random rate matrix R and transition matrix A
    % if vargin>2: generate parameter dependent rate matrix R and transition matrix A
    
    model=varargin{1};
    DT=varargin{2};

    
    switch model
        case 'twostate'
            if nargin==2
                sumk=10^(3*rand(1)-4);
                pon=rand()*0.95+0.05;
            else
                sumk=varargin{3}(1);
                pon=varargin{3}(2);
            end
            kon=sumk*pon;
            koff=sumk*(1-pon);
            R = [-kon koff;...
                kon -koff];
            A = rate_to_prob(R,DT);
        case 'threestatecircle'
            if nargin==2
                sumk=10^(3*rand(1)-4);
                pon=rand()*0.95+0.05;
                k1k2=rand();
            else
                sumk=varargin{3}(1);
                pon=varargin{3}(2);
                k1k2=varargin{3}(3);
            end
            k1=sumk*pon*(k1k2+1);
            k2=sumk*pon*(k1k2+1)/k1k2;
            koff=sumk*(1-pon);
            R = [-k1 0 koff ; k1 -k2  0; 0 k2 -koff];
            A = rate_to_prob(R,DT);
        case 'twocopycorr'
            if nargin==2
                sumk=10^(3*rand(1)-4);
                pon=rand()*0.95+0.05;
            else
                sumk=varargin{3}(1);
                pon=varargin{3}(2);
            end
            kon=sumk*pon;
            koff=sumk*(1-pon);
            R = [-kon koff;...
                kon -koff];
            A = rate_to_prob(R,DT);
        case 'twocopyuncorr'
            if nargin==2
                sumk=10^(3*rand(1)-4);
                pon=rand()*0.95+0.05;
            else
                sumk=varargin{3}(1);
                pon=varargin{3}(2);
            end
            kon=sumk*pon;
            koff=sumk*(1-pon);
            R = [-2*kon koff 0; ...
                2*kon -kon-koff 2*koff;...
                0 kon -2*koff];
            A = rate_to_prob(R,DT);
        case 'twocopyalpha'
            if nargin==2
                sumk=10^(3*rand(1)-4);
                pon=rand()*0.95+0.05;
                alpha=rand()*20;
            else
                sumk=varargin{3}(1);
                pon=varargin{3}(2);
                alpha=varargin{3}(3);
            end
            kon=sumk*pon;
            koff=sumk*(1-pon);
            R = [-2*kon alpha*koff 0; ...
                2*kon -alpha*kon-alpha*koff 2*koff;...
                0 alpha*kon -2*koff];
            A = rate_to_prob(R,DT);
        case 'threestate'   % TO BE UPDATED
            if nargin==2
                A_=rand(2,3);
                A(1,:)=A_(1,:).*A_(2,:);
                A(2,:)=(1-A_(1,:)).*A_(2,:);
                A(3,:)=1-A(1,:)-A(2,:);
                R = prob_to_rate(A,varargin{2});
            else
                
            end
    end