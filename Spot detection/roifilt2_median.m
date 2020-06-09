function I = roifilt2_median(J,BW,mn)
    %Determine rectangle that encloses
    % the non-zero elements in BW.  The rectangle vector is chosen so that
    % no non-zero element in J is considered to be a boundary pixel by
    % imfilter.  In other words, the row and column padding should be equal
    % to the row size and column size of H, respectively.  Also, rectangle
    % cannot be bigger than size of original image.
    
    [row, col] = find(BW==1);
    colpad = ceil(mn(2) / 2);
    rowpad = ceil(mn(1) / 2);
    mincol = max(1, min(col(:)) - colpad);
    minrow = max(1, min(row(:)) - rowpad); 
    maxcol = min(size(J, 2), max(col(:)) + colpad);
    maxrow = min(size(J, 1), max(row(:)) + rowpad);
    
    % perform filtering on y that is cropped to the rectangle.
    I = J;
    J = J(minrow:maxrow, mincol:maxcol);
    BW(minrow:maxrow, mincol:maxcol)=true;
    filtI = medfilt2(J,mn);

    I(BW) = filtI;
