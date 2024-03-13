function sp = synthetic_spikes(rho,dt,m)
% rho = 0.5;
% rho = linspace(0,1,21);
% for i = 1:length(rho)
%     [M1(:,i),M3(:,i)] = get_spike_correlation(rho(i));
% end

% figureNB;
%     plot(M1,M3,'.k','MarkerSize',10)
%     FT = fitlm(M1(:),M3(:),'intercept',false)
%     t = linspace(0,1,1e3);
%     hold on;
%     plot(t,FT.predict(t(:)),'r','LineWidth',1)
%     xlabel('Spike time correlation')
%     ylabel('apEEG correlation')
%     gcaformat

% function [M1,M3] = get_spike_correlation(rho)
    N = 30;
    tmax = 10;
    % dt = 1/4e3;
    % rho = 0.2;

    % [uAP,mtype] = getUnitaryAP;
    % uAP = resample(uAP,1/dt,16e3);
    [L0,gam] = compute_DG(N,tmax,dt,rho);

    t = 0:1/16e3:tmax;
    sp = zeros(length(t)-1,N,m);
    for rep = 1:m
        [ids,ts] = synthesize_spikes(N,tmax,dt,L0,gam);
        % ids = randi(N,5*tmax*N,1);
        % ts = tmax*rand(5*tmax*N,1);
        % uAP_sampled = sampleAPs(uAP,mtype,N);

        % apEEG = zeros(length(t)-1,N);
        for i = 1:N
            sp(:,i,rep) = histcounts(ts(ids==i),t);
            % apEEG(:,i) = filter(uAP_sampled(:,i),1,sp(:,i));
        end

        % R1 = 1-pdist(sp','correlation');
        % R3 = 1-pdist(apEEG','correlation');
        % idcs = find(R1>=0);
        % idcs = 1:length(R1);
        % M1(rep) = mean(R1(idcs));
        % M3(rep) = mean(R3(idcs));
    end
end

function [L,gam] = compute_DG(N,tmax,dt,rho)
    % Firing rate (1 Hz)
    p = dt*ones(N,1);

    % Max correlation
    R = p(1)*(1-p(1))*rho;

    % Use haversine distance among synapses to construct covariance matrix
    KxxCell = cell(N,3); % Store sparse cov matrix (i,j,v) for every synapse
    for i = 1:N
        % d = hav_d(X(i,:),X(i+1:end,:));
        % idcs = find(d<0.25);
        % R = R*exp(-10*d(idcs)).*sqrt(rVar(i)*rVar(i+idcs));
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
    for k = 1:length(i0)
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
function [uAP,mtype] = getUnitaryAP
    load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\unitaryAP_all.mat');
    [sa,X] = network_simulation_beluga.getHeadModel;
    uAP = network_simulation_beluga.getEEG(savedUnitaryAP,sa,43e3);
end
function uAP_sampled = sampleAPs(uAP,mtype,N)
    load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\mtype_abundance.mat','mtype_abundance')
    [J,ID] = findgroups(mtype);
    count = splitapply(@(x) sum(x),ones(size(J)),J)';
    scaledAbundance = mtype_abundance.Abundance./count;
    prop = scaledAbundance(J);
    F = cumsum(prop)/sum(prop);
    [~,idcs] = unique(F);
    R = rand(N,1);
    cellID = interp1(F(idcs),idcs,R,'next','extrap');

    uAP_sampled = uAP(:,cellID);
end