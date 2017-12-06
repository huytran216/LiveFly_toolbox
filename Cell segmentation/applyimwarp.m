function [Inew]=applyimwarp(Iold,tform,alignratio)
% Apply scale up and down with alignratio and imwarp in between
    Itmp_=imresize(Iold,alignratio,'nearest');
    [Itmp,rd]=imwarp_poly(Itmp_,tform);
    Itmp=resetframe(Itmp,rd,size(Itmp_));
    Inew=imresize(Itmp,size(Iold),'nearest');