function [llh,std]=get_loglikelihood(fout,totalsample)
    % Mean square end:
    std=sqrt(fout/totalsample);  % Mean square of difference of totalsample sample
    % Likelihood
    llh = -totalsample/2 - totalsample*1/2*log(2*pi*std^2);
end