function [error,Iout] = fitall_(xin,xcrit,Irec_,time_rec_,fun_,L)
    xout = xin(xcrit);
    error =0;
    if ~exist('L','var')
        L=1;
    end
    Iout = cell(1,numel(Irec_));
    if (size(Irec_{1},1)>1)&&(size(Irec_{1},2)>1)
        % Fit based on individual curves
        for i=1:numel(Irec_)
            if any(isnan(Irec_{i}))
                dt = time_rec_{i}(2)-time_rec_{i}(1);
                tmpfun = fun_(xout([i*2+[-1 0] end-numel(Irec_)*2+i end-numel(Irec_)+i]),[time_rec_{i} time_rec_{i}(end)+(1:numel(L))*dt]);

                tmpfun = conv(tmpfun,L,'full');
                tmpfun = tmpfun(1:end-2*numel(L)+1);
                for j=1:size(Irec_{i},1)
                    error = error + sum((Irec_{i}(j,:)-tmpfun).^2);
                end
                Iout{i}=tmpfun;
            end
        end
    else
        % Fit based on mean curve
        for i=1:numel(Irec_)
            if any(isnan(Irec_{i}))
                error = error + sum((Irec_{i}-fun_(xout([i*2+[-1 0] end-numel(Irec_)*2+i end-numel(Irec_)+i]),time_rec_{i})).^2);
            end
        end
    end
end