function Iout = dualcolormap(I1,I2,r1,r2)
    % I1: image 1   % 2D grayscale
    % I2: image 2   % 2D grayscale
    % Iout: image combined   % 2D color
    % r1: range1 (min max)
    % r2: range2 (min max)
    
    I1 = (I1-r1(1))/(r1(2)-r1(1));
    I2 = (I2-r2(1))/(r2(2)-r2(1));
    
    c1 = rgb2hsv(double([1,0,0]));
    c2 = rgb2hsv(double([0,1,0]));
    
    
    I1_ = hsv2rgb(cat(3,I1*0+c1(1),I1,I1*0+1));
    I2_ = hsv2rgb(cat(3,I2*0+c2(1),I2,I2*0+1));
    
    Iout = I1_/2 + I2_/2;
%     subplot(131);
%     imshow(I1_);
%     subplot(132);
%     imshow(I2_);

%     subplot(133);
%     imshow(Iout);