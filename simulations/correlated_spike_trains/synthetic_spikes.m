function sp = synthetic_spikes(rho,dt,N)
    tmax = 10;
    t = 0:1/16e3:tmax;
    sp = zeros(length(t)-1,N);

    [L0,gam] = compute_DG(N,tmax,dt,rho);
    [ids,ts] = synthesize_spikes(N,tmax,dt,L0,gam);
    for i = 1:N
        sp(:,i) = histcounts(ts(ids==i),t);
    end
end

function [L,gam] = compute_DG(N,tmax,dt,rho)
    % Firing rate (1 Hz)
    p = dt*ones(N,1);

    % Max correlation
    R = p(1)*(1-p(1))*rho;

    KxxCell = cell(N,3); % Store sparse cov matrix (i,j,v) for every synapse
    parfor i = 1:N
        idcs = (i+1:N)';
        KxxCell(i,:) = {idcs,i+0*zeros(size(idcs)),R*ones(size(idcs))};
    end
    i0 = cat(1,KxxCell{:,1});
    j0 = cat(1,KxxCell{:,2});
    S = cat(1,KxxCell{:,3});

    % Find mu and sig of n-d Gaussian with which to sample spikes
    gam = icdf('normal',p,0,1);
    fun = @(x,kk,rr,gg) kk+rr-mvncdf(gg,0,[1,x;x,1]);
    A = zeros(length(i0),1);
    ibnds = [-1+1e-4,1-1e-4];
    parfor k = 1:length(i0)
        i = i0(k);
        j = j0(k);
        A(k) = fzero(@(x) S(k)+p(i)*p(j)-mvncdf(gam([i,j]),0,[1,x;x,1]),ibnds);
    end
    L = sparse(i0,j0,A,N,N);
    L = L+L'+eye(N);
end

function [ids,ts] = synthesize_spikes(N,tmax,dt,L0,gam)
    % Sample spikes using a discretized Gaussian
    M = ceil(tmax/dt);
    Lc = chol(nearcorr(L0),'lower');
    xCell = cell(M,2);
    for i = 1:M
        ids0 = find((Lc*randn(N,1))<gam);
        ts0 = (i*ones(size(ids0))+rand(size(ids0)))*dt;
        xCell(i,:) = {ids0,ts0};
    end
    ids = cat(1,xCell{:,1});
    ts = cat(1,xCell{:,2});
end