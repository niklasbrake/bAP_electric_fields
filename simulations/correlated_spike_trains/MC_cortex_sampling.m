%{
% Load precomputed cross spectra
load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\spike_synchrony\cross_spectra.mat')

% Calculate new coordinates uiAi
[sa,X] = network_simulation_beluga.getHeadModel;
M = size(X.vertices,1);
Lxyz = zeros(M,3);
for idx = 1:M
    L0 = squeeze(sa.cortex75K.V_fem(49,idx,:))'; % czIDX = 49
    vz = sa.cortex75K.normals(idx,:);
    [vx,vy] = getOrthBasis(vz);
    A = [vx(:),vy(:),vz(:)];
    Lxyz(idx,:) = 1e-6*L0*A;
end
% Set up a KD tree to efficiently search for nearby vertices
Mdl = KDTreeSearcher(X.vertices);
%}



% Result variables
Sxy = zeros(100,1); Pxy = zeros(100,1);
Sxx = zeros(100,1); Pxx = zeros(100,1);
% dValues = [0,10.^linspace(-3,2,1e3)];
dValues = linspace(0,100,1e3);
Pxy_D = zeros(100,length(dValues));
SSExy_D = zeros(100,length(dValues));
SSE_xx = zeros(100,1);
count_D = zeros(length(dValues),1);

% Convergence variable
iConverge = 1;
maxLag = 1e3;
maxMeanChange = zeros(maxLag,1);

% Sample neurons proportional to abundance
iSampling = 1;
samplingInterval = 1e4;
neuronIdcs = reshape(sample_blue_neurons(2*samplingInterval),samplingInterval,2);
mapIJ = @(i,j) (i-1)*1035-(i-1)*i/2+j-i; % Map neuron ids to Rij dimension

N = 0;


figureNB;
f = 10.^linspace(0,log10(3e3),100);
x = [f,flip(f)];
subplot(1,2,1);
    SE = nan+sqrt(Sxx/(N-1))/sqrt(N);
    mu = nan+Pxx/N;
    y = [mu-1.96*SE;flip(mu)+1.96*flip(SE)];
    plt2 = fill(x,y,'k','FaceAlpha',0.2);
    set(gca,'xscale','log');
    ttl2 = title('');
subplot(1,2,2);
    SE = nan+sqrt(sum(SSExy_D,2)/(N-1))/sqrt(N);
    mu = nan+sum(Pxy_D,2)/N;
    y = [mu-1.96*SE;flip(mu)+1.96*flip(SE)];
    plt1 = fill(x,y,'k','FaceAlpha',0.2);
    set(gca,'xscale','log');
    ttl1 = title('');
while true
    N = N+1;

    % Randomly sample unitary AP profiles
    i = min(neuronIdcs(iSampling,:));
    j = max(neuronIdcs(iSampling,:));
    Rij = Sij(:,:,:,mapIJ(i,j));

    % Randomly sample vertex from cortex
    iX = randi(M);

    % Sample second vertex with higher probability for nearby neurons
    d = exprnd(8);
    idcs = rangesearch(Mdl,X.vertices(iX,:),d);
    jX = randsample(idcs{1},1);


    uiAi = Lxyz(iX,:);
    ujAj = Lxyz(jX,:);

    for iF = 1:100
        Sxy(iF) = uiAi*Rij(:,:,iF)*ujAj';
        Sxx(iF) = 0.5*(uiAi*Sii(:,:,iF,i)*uiAi'+ujAj*Sii(:,:,iF,j)*ujAj');
    end
    % Scale by distance (spike time correlation)
    Sxy = exp(-d.^2./6)*Sxy;

    % Update mean
    et = Sxy-Pxy;
    Pxy = Pxy + et/N;
    etxx = (Sxx-Pxx);
    Pxx = Pxx + etxx/N;
    SSE_xx = SSE_xx + etxx.*(Sxx-Pxx);

    % Update d dependence
    iD = interp1(dValues,1:length(dValues),d,'nearest','extrap');
    Pxy_D(:,iD) = Pxy_D(:,iD) + Sxy;
    SSExy_D(:,iD) = SSExy_D(:,iD) + et.*(Sxy-Pxy);
    count_D(iD) = count_D(iD)+1;

    % Convergence check
    if(iConverge==maxLag)
        SE = sqrt(sum(SSExy_D,2)/(N-1))/sqrt(N);
        mu = sum(Pxy_D,2)/N;
        plt1.YData = [mu-1.96*SE;flip(mu)+1.96*flip(SE)];
        ttl1.String = sprintf('%d (%f%c)',N,100*N/1e15,char(37));

        SE = sqrt(SSE_xx/(N-1))/sqrt(N);
        mu = Pxx;
        plt2.YData = [mu-1.96*SE;flip(mu)+1.96*flip(SE)];
        ttl2.String = sprintf('%d (%f%c)',N,100*N/76985370,char(37));
        drawnow;
        if(max(1.96*2*SE)<1e-17);
            break;
        end
        % fprintf('Convergence index = %.3e (N = %d, %.5f%c of possible pairs)\n',mean(maxMeanChange),N,N/1e15*100,char(37));
        iConverge = 0;
    end
    iConverge = iConverge + 1;

    % Resample check
    if(iSampling==samplingInterval)
        neuronIdcs = reshape(sample_blue_neurons(2*samplingInterval),samplingInterval,2);
        iSampling = 0;
    end
    iSampling = iSampling+1;
end







nCount = 16e9;

figureNB;
plot(f,8.5*nCount*Pxx + 8.5*nCount*(nCount-1)*Pxy)
