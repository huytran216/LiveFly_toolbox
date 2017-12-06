function [Iout,rd]=imwarp_poly(Iin,tform)
    if strcmp(class(tform),'projective2d')
        % The input is actually the displacement matrix
        [Iout,rd]=imwarp(Iin,tform,'nearest');
    else
        % If input is actually a polynomial
        [x,y]=meshgrid(1:size(Iin,2),1:size(Iin,1));
        xout=reshape(tform.x([x(:),y(:)]),size(Iin));
        yout=reshape(tform.y([x(:),y(:)]),size(Iin));
        Iout=cat(3,xout,yout);
        [Iout,rd]=imwarp(Iin,Iout,'nearest');
    end

    