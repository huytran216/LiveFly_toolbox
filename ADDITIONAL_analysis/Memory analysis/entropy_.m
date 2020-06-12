function E = entropy_(x)
%ENTROPY Entropy of intensity image.    
%   E = ENTROPY(I) returns E, a scalar value representing the entropy of an
%   intensity image.  Entropy is a statistical measure of randomness that can be
%   used to characterize the texture of the input image.  Entropy is defined as
%   -sum(p.*log2(p)) where p contains the histogram counts returned from IMHIST.

%   Modified from the original entropy to accomodate binned data - Huy Tran (22-4-2016)
if numel(unique(x))<1
    E=NaN;
    return;
end
if numel(unique(x))==1
    E=0;
    return;
end

% calculate histogram counts
p = hist(x,unique(x));

% normalize p so that sum(p) is one.
p = p ./ sum(p);

E = -sum(p.*log2(p));
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function I = ParseInputs(varargin)

narginchk(1,1);

validateattributes(varargin{1},{'uint8','uint16', 'double', 'logical'},...
              {'real', 'nonempty', 'nonsparse'},mfilename, 'I',1);

I = varargin{1};
