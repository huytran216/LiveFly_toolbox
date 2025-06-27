function [betahat_,fhat_] = fminsearch_global(fun,beta_lb,beta_ub,nitr)
    fhat_ = 1e10;
    betahat_ = beta_lb;
    for itr = 1:nitr
        beta0 = rand(size(beta_lb)).*(beta_ub-beta_lb)+beta_lb;
        [betahat,fhat]=fminsearchbnd(fun,beta0,beta_lb,beta_ub);
        if fhat<fhat_
            fhat_ = fhat;
            betahat_ = betahat;
        end
    end
end