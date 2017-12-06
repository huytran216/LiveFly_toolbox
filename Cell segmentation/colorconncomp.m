function CC = colorconncomp(BW)
%BWCONNCOMP Find connected components in color image.
%   CC = COLORCONNCOMP(BW) returns the connected components CC found in BW. 
%   BW is a non-binary image that can have any dimension. CC is a structure 
%   with four fields:
%
%      Connectivity   Connectivity of the connected components (objects).
%
%      ImageSize      Size of BW.
%
%      NumObjects     Number of connected components (objects) in BW.
%
%      PixelIdxList   1-by-NumObjects cell array where the kth element
%                     in the cell array is a vector containing the linear
%                     indices of the pixels in the kth object.
%    
%   Default connectivity =8

CC = struct(...
    'Connectivity', 0, ...
    'ImageSize', size(BW), ...
    'NumObjects', [], ...
    'PixelIdxList', []);

idxlist=unique(BW(:));
cnt=0;
for i=1:numel(idxlist)
    if idxlist(i)~=0
        cnt=cnt+1;
        CC.PixelIdxList{cnt}=find(BW==idxlist(i));
    end
end
CC.Connectivity=0;
CC.NumObjects=sum(idxlist>0);
CC.ImageSize=size(BW);