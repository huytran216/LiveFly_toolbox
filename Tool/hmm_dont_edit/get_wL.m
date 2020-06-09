function wL=get_wL(Nbs,padsize,DT)
% Generate the filtering window wL based on the reporter gene configuration:
% INPUT:
%   Nbs: number of binding sites
%   padsize: size of the subsequent gene at the 5" of the reporter gene.
% OUTPUT:
%   wL: wL(i)=size of the nascent RNA after time DT*i
    if (Nbs==1)&&(padsize==0)
        wL=1;
    else
        ms2spacing=52;  % spacing between ms2 binding sites
        ke=45;          % Elongation rate: bp/s
        ms=zeros(1,Nbs*ms2spacing+padsize)+Nbs;
        for i=1:Nbs
            ms((i-1)*ms2spacing+1:i*ms2spacing)=i;
        end
        sizeDT=DT*ke;    % Polymerase size
        for i=1:numel(ms)/sizeDT
            wL(i)=sum(ms((sizeDT*(i-1)+1):(sizeDT*i)))/sizeDT;
        end
        if i<numel(ms)/sizeDT
            wL(i+1)=sum(ms((sizeDT*(i)+1):end))/sizeDT;
        end
    end
    wL=wL/wL(end);
end
