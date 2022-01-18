function [I,Ibg,pval]=randmi(A,B,varargin)
%RANDMI Determines the mutual information of two images or signals
%
%   I=mi(A,B)   Mutual information of A and B.
%   Ibg: Mutual information of A and B, if A is permuted randomly
%   pval: p-val: probablity that I is generated from a shuffled set
%
%   Assumption: 0*log(0)=0
%
%   See also ENTROPY.

%   jfd, 15-11-2006
%        01-09-2009, added case of non-double images
%        24-08-2011, speed improvements by Andrew Hill
itrmax=1; % Set to 200 for better bootstrap. Set to 1 to ignore
I=mi(A,B);
for itr=1:itrmax
    Ibg(itr)=mi(A(randperm(numel(A))),B);
end
if ~isnan(I)
    pval=sum(Ibg>=I)/itrmax;
else
    pval=1;
end
Ibg=mean(Ibg);
