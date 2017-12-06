function Inew=resetframe(Iold,rd,Osize)
    % Convert image from world coordinate to intrinsic coordinate
    % Input:
        % Iold: Old image in world coordinate
        % rd: image world coordinate reference
        % Osize: Size of original image
    % Real world coodinate
    rdx=round(max(1,1-rd.XWorldLimits(1)):min(rd.ImageSize(2),Osize(2)-rd.XWorldLimits(1)));
    rdy=round(max(1,1-rd.YWorldLimits(1)):min(rd.ImageSize(1),Osize(1)-rd.YWorldLimits(1)));
    % Intrinsic coodinate
    rdx_=round(max(1,1+rd.XWorldLimits(1)):min(rd.ImageSize(2)+rd.XWorldLimits(1),Osize(2)));
    rdy_=round(max(1,1+rd.YWorldLimits(1)):min(rd.ImageSize(1)+rd.YWorldLimits(1),Osize(1)));
    Inew=zeros(Osize);Inew(rdy_,rdx_)=Iold(rdy,rdx);
end