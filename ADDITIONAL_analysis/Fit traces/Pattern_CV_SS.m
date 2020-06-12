function [CV,S2,S1]=Pattern_CV_SS(A,XPDF,XON,trange)
    % A: transition rate matrix
    % XPDF: steady state distribution
    % XON: emission vector
    % trange: trange

    A=double(A);
    
    % Get xfire
    XPDF=XPDF(:).*XON(:);
    % Diagonalize A:
    [V,D] = eig(A);
    % Calculate lambda
    L=cell(1,numel(trange));
    for j=1:numel(trange)
        L{j}=zeros(numel(XON));
        for i=1:size(D,1)
            if abs(D(i,i))>=1e-5
                L{j}(i,i)=-trange(j)/D(i,i)+1/D(i,i)^2*(exp(D(i,i)*trange(j))-1);
            else
                L{j}(i,i)=trange(j)^2/2;
            end
        end
    end
    % Calculate <S^2>
    S2=cellfun(@(x) 2*XON(:)'*V*x*(V^(-1))*XPDF,L);
    % Calculate <S>
    S1=sum(XPDF)*trange;
    % Calculate noise:
    CV=(S2./S1.^2)-1;
   