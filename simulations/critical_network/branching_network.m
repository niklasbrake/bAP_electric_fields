function varargout = simulatespikes(m,nE,nI,tmax)
    if(m>1)
        error('Branching index must be less than 1.')
    end
    tmax = tmax*1e-3;
    dt = 4e-3;
    t = 0:dt:tmax;
    N = length(t);

    eiFraction = nE/(nE+nI);

    k = 4;
    h = dt*nE*(1-m);

    % Generate mean firing rate using critical branching process
    B = zeros(N,1);
    exN = poissrnd(h,N,1);
    B(1) = poissrnd(dt*nE);
    for i = 1:N-1
        if(B(i))
            count = binornd(B(i)*k,m/k);
        else
            count = 0;
        end
        B(i+1) = count+exN(i);
    end
    if(nargout==2)
        varargout = {t,B};
        return;
    end

    numspikes = sum(B);
    spikeTime = zeros(5*numspikes,1);
    neuronIDs = zeros(5*numspikes,1);
    eiRatio = 5;
    j = 0;
    for i = 1:length(t)
        nEx = poissrnd(B(i));
        spikeTime(j+1:j+nEx) = t(i) + rand(1,nEx)*dt - dt/2;
        neuronIDs(j+1:j+nEx) = randperm(nE,nEx);
        j = j+nEx;

        nIn = poissrnd(B(i)*eiRatio*(1-eiFraction));
        spikeTime(j+1:j+nIn) = t(i) + rand(1,nIn)*dt - dt/2;
        neuronIDs(j+1:j+nIn) = nE+randperm(nI,nIn);
        j = j+nIn;
    end
    spikeTime(j+1:end) = [];
    neuronIDs(j+1:end) = [];
    ei = [zeros(1,nE),ones(1,nI)];
    [neuronIDs,I] = sort(neuronIDs);
    spikeTime = 1e3*spikeTime(I);

    varargout = {neuronIDs,spikeTime,ei,t,B};
end