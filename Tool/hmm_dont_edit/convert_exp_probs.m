%Convert Probability Matrices from General 3 State HMM to Rates
%Manually Transferred prob matrices
dT = 40.8;
path = 'C:\Users\Nicholas\OneDrive\Berkeley Biophysics\Garcia_Fall 2016\Transcriptional Bursting\HMM Model\HMM Output\';
filename = 'AP 40-46 2016_10_17';

A(:,:,1) = [0.7682,0.2326,0.1807; 0.2318,0.6502,0.5441;1.622e-13,0.1173,0.2751]; %40
A(:,:,2) = [0.7595,0.2533,0.06471; 0.2405,0.6207,0.6579;1.123e-15,0.126,0.2774]; %41
A(:,:,3) = [0.7208,0.2779,0.02346;0.2792,0.5924,0.7575;2.885e-40,0.1297,0.2191]; %42
A(:,:,4) = [0.7552,0.2981,0.1567;0.2445,0.6698,0.4589;0.0003535,0.03209,0.3843]; %43
A(:,:,5) = [0.8126,0.1634,0.3037;5.092e-21,0.4054,0.1128;0.1874,0.4312,0.5835]; %44
A(:,:,6) = [0.7793,0.3276,0.2963;0.1677,0.5843,0.3576;0.05299,0.08807,0.3461]; %45

R = zeros(6,10);
R_const = zeros(6,10);
A_out = zeros(6,10);
for i = 1:6
    R(i,1) = i+39;
    R(i,2:end) = reshape(rate_fit(A(:,:,i),dT,0,.005),1,[]);
    R_const(i,1) = i+39;
    R_const(i,2:end) = reshape(rate_fit(A(:,:,i),dT,1,.005),1,[]);
    A_out(i,1) = i+39;
    A_out(i,2:end) = reshape(A(:,:,i),1,[]);
end;

write = strcat(path,filename,'.xlsx');

xlswrite(write,R, 'Rate Fits General');
xlswrite(write,R_const, 'Rate Fits 2-State');
xlswrite(write,A_out, 'Probs');







