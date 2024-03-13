function c = CCG(x1,x2,lag,dt)
    if(length(x1)~=length(x2))
        error('Spike trains must be the same length');
    end

    lag = lag/dt;
    lambda1 = sum(x1)/(dt*length(x1));
    lambda2 = sum(x2)/(dt*length(x2));
    c = xcorr(x1,x2,lag,'unbiased');
    c = c/sqrt(lambda1*lambda2);
end