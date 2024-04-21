function MC_cortex_sampling_Rxy(chain_number)

saveFile = sprintf('/lustre04/scratch/nbrake/data/simulation_analyzed/cross_spectra/chain%s.mat',chain_number);

% Load anatomy information
load('/lustre04/scratch/nbrake/data/simulation_analyzed/cortical_area/MC_data.mat','A','dValues');
% Load savedUnitaryAP
load('/lustre04/scratch/nbrake/data/simulation_analyzed/unitaryAP/unitaryAPNew.mat','savedUnitaryAP')
savedUnitaryAP = permute(savedUnitaryAP,[2,1,3]);
% Calculate new coordinates uiAi
load('/lustre04/scratch/nbrake/data/anatomy_nyhead_model.mat');
X = struct();
X.vertices = sa.cortex75K.vc;
X.faces= sa.cortex75K.tri;
Cz = find(strcmp(sa.clab_electrodes,'Cz'));

% Calculate normal lead fields
M = size(X.vertices,1);
Lxyz = zeros(M,3);
for idx = 1:M
    L0 = squeeze(sa.cortex75K.V_fem(Cz,idx,:))'; % czIDX = 49
    vz = sa.cortex75K.normals(idx,:);
    [vx,vy] = getOrthBasis(vz);
    Lxyz(idx,:) = 1e-6*L0*[vx(:),vy(:),vz(:)];
end

% Set up a KD tree to efficiently search for nearby vertices
Mdl = KDTreeSearcher(X.vertices);

% Result variables
L = 16e3;
fs = 16e3;
f = fs/L:fs/L:fs/2;
Sxy = zeros(length(f),1); Pxy = zeros(length(f),1);
Pxy_D = zeros(length(f),length(dValues));
SSExy_D = zeros(length(f),length(dValues));
count_D = zeros(1,length(dValues));

% Convergence variable
iLargeLoop = 1;
loopSize = 1e3;
d_abs = 1e-4/(16e9)^2*sqrt(40);

% Sample neurons proportional to abundance
iSampling = 1;
samplingInterval = 1e4;
neuronIdcs = reshape(sample_blue_neurons(2*samplingInterval),samplingInterval,2);
mapIJ = @(i,j) (i-1)*1035-(i-1)*i/2+j-i; % Map neuron ids to Rij dimension

figureNB;
x = [f,flip(f)];
    SE = nan+f(:);
    mu = nan+f(:);
    y = [mu-1.96*SE;flip(mu)+1.96*flip(SE)];
    plt1 = fill(x,y,'k','FaceAlpha',0.2,'EdgeColor','none');
    hold on;
    plt0 = plot(f,mu*nan,'k','LineWidth',1);
    set(gca,'xscale','log');
    ttl1 = title('');

N = 0;
tic
while true
    N = N+1;

    % Randomly sample unitary AP profiles
    i = neuronIdcs(iSampling,1);
    j = neuronIdcs(iSampling,2);

    % Randomly sample vertex from cortex
    iX = randi(M);
    % Sample second vertex with higher probability for nearby neurons
    d = exprnd(6);
    idcs = rangesearch(Mdl,X.vertices(iX,:),d);
    if(length(idcs{1})>1)
        jX = randsample(idcs{1},1);
    else
        jX = idcs{1};
    end
    d = norm(X.vertices(iX,:)-X.vertices(jX,:));
    uiAi = Lxyz(iX,:);
    ujAj = Lxyz(jX,:);

    % Compute eeg cross spectrum
    eeg1 = uiAi*savedUnitaryAP(:,:,i);
    eeg2 = ujAj*savedUnitaryAP(:,:,j);
    R_eeg = xcorr(eeg2,eeg1,8e3,'unbiased');
    R_eeg = detrend(R_eeg(2:end),'constant');
    S12 = fft(R_eeg)/fs^2*L; % units = per spike rather than per Hztime
    Sxy = 2*S12(2:8001).*exp(2*pi*sqrt(-1)*f*7999/fs);
    Sxy = real(Sxy(:));

    % Update d dependence
    et = Sxy-Pxy;
    Pxy = Pxy + et/N;
    iD = interp1(dValues,1:length(dValues),d,'nearest','extrap');
    Pxy_D(:,iD) = Pxy_D(:,iD) + Sxy;
    SSExy_D(:,iD) = SSExy_D(:,iD) + et.*(Sxy-Pxy);
    count_D(iD) = count_D(iD)+1;

    % Convergence check
    if(iLargeLoop==loopSize)
        SIG = SSExy_D./(count_D-1);
        SIG(isinf(SIG)) = nan;
        SIGN = nansum(A.*SIG,2);
        STOP = min(0.99-SIGN/N/d_abs^2);

        % mu = Pxy_D./count_D;
        % Rxy = nansum(A.*mu,2);


        % CI = 1.96*sqrt(SIGN)./sqrt(N);
        % Rxy_CI = [Rxy-CI;flip(Rxy+CI)];

        % plt0.YData = Rxy;
        % plt1.YData = Rxy_CI;
        % ttl1.String = sprintf('%d (%f)',N,STOP);
        % drawnow;
        iLargeLoop = 0;
        if(STOP>=0 || N>15e3)
            break;
        end
    end
    iLargeLoop = iLargeLoop + 1;

    % Resample check
    if(iSampling==samplingInterval)
        neuronIdcs = reshape(sample_blue_neurons(2*samplingInterval),samplingInterval,2);
        iSampling = 0;
    end
    iSampling = iSampling+1;
end


save(saveFile,'dValues','Pxy_D','SSExy_D','count_D','f');