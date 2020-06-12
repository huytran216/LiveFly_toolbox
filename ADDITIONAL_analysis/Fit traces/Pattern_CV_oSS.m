function [CV,S2,S1]=Pattern_CV_oSS(A,XPDF,XON,trange)
    % A: transition rate matrix
    % XPDF: steady state distribution
    % XON: emission vector
    % trange: trange

    A=double(A);
    
    XON=XON(:);
    XPDF=XPDF(:);
    % Do stuffs numerically
    NT=100; % number of time bins
    
    S1=[];
    S2=[];
    for j=1:numel(trange)
        expm_rec=cell(1,NT);
        T=trange(j);
        % Find second moment
        tmp=0;
        urange=linspace(0,T,NT);
        for u=1:numel(urange)
            % Find exponential matrix
            expm_rec{u}=expm(A*urange(u));
        end
        S1_=zeros(1,NT);
        S2_=zeros(NT,NT);
        
        for u=1:numel(urange)
            % Find first moment
            S1_(u)=XON'*expm(A*urange(u))*XPDF*T/NT;
            % Find second moment
            for v=1:u-1
                p=expm_rec{v}*XPDF;
                m=expm_rec{u-v};
                % Find mean value <n(t)*n(s)>
                S2_(u,v)=sum(XON*XON'.*m*p)*T^2/NT^2;
            end
        end
        S1(j)=sum(S1_)/T;
        S2(j)=2*sum(S2_(:))/T^2;
    end
    % Calculate noise:
    CV=(S2./S1.^2)-1;
   