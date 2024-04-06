addpath 'C:\Users\brake\Documents\GitHub\bAP_electric_fields\simulations\calculate_unitary_bAP'

% Load precomputed cross spectra
load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\unitaryAPNew.mat')
savedUnitaryAP = permute(savedUnitaryAP,[2,1,3]);

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


% Result variables
L = 2e3;
fs = 16e3;
f = fs/L:fs/L:fs/2;
Sxy = zeros(length(f),1); Pxy = zeros(length(f),1);
Sxx = zeros(length(f),1); Pxx = zeros(length(f),1);
% dValues = [0,10.^linspace(-3,2,1e3)];
dValues = linspace(0,50,1e3);
Pxy_D = zeros(length(f),length(dValues));
SSExy_D = zeros(length(f),length(dValues));
SSE_xx = zeros(length(f),1);
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
x = [f,flip(f)];
subplot(1,2,1);
    SE = nan+sqrt(Sxx/(N-1))/sqrt(N);
    mu = nan+Pxx/N;
    y = [mu-1.96*SE;flip(mu)+1.96*flip(SE)];
    hold on;
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
    uiAi = Lxyz(iX,:);
    ujAj = Lxyz(jX,:);

    % Compute eeg cross spectrum
    eeg1 = uiAi*savedUnitaryAP(:,:,i);
    eeg2 = ujAj*savedUnitaryAP(:,:,j);
    R_eeg = xcorr(eeg2,eeg1,1e3,'unbiased');
    R_eeg = detrend(R_eeg(2:end),'constant');
    S12 = fft(R_eeg)/fs^2*L; % units = per spike rather than per Hztime
    Sxy = 2*S12(2:1001).*exp(2*pi*sqrt(-1)*f*999/fs);
    % Scale by distance (spike time correlation)
    Sxy = real(Sxy(:));

    % Compute eeg unitary spectrum
    R_eeg = xcorr(eeg1,eeg1,1e3,'unbiased');
    R_eeg = detrend(R_eeg(2:end),'constant');
    S11 = fft(R_eeg)/fs^2*L; % units = per spike rather than per Hztime
    Sxx = 2*S11(2:1001).*exp(2*pi*sqrt(-1)*f*999/fs);
    Sxx = real(Sxx(:));

    % Update mean
    et = Sxy-Pxy;
    Pxy = Pxy + et/N;
    etxx = (Sxx-Pxx/N);
    Pxx = Pxx + Sxx;
    SSE_xx = SSE_xx + etxx.*(Sxx-Pxx/N);

    % Update d dependence
    iD = interp1(dValues,1:length(dValues),d,'nearest','extrap');
    Pxy_D(:,iD) = Pxy_D(:,iD) + Sxy;
    SSExy_D(:,iD) = SSExy_D(:,iD) + et.*(Sxy-Pxy);
    count_D(iD) = count_D(iD)+1;

    % Convergence check
    if(iConverge==maxLag)
        SE = sqrt(SSE_xx/(N-1))/sqrt(N);
        mu = Pxx/N;
        plt2.YData = [mu-1.96*SE;flip(mu)+1.96*flip(SE)];
        ttl2.String = sprintf('%d (%f%c)',N,100*N/76985370,char(37));

        SE = sqrt(sum(SSExy_D,2)/(N-1))/sqrt(N);
        mu = sum(Pxy_D,2)/N;
        plt1.YData = [mu-1.96*SE;flip(mu)+1.96*flip(SE)];
        ttl1.String = sprintf('%d (%f%c)',N,100*N/1e15,char(37));
        drawnow;
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



N2 = 16e9;
dataFolder = 'E:\Research_Projects\004_Propofol\manuscript\Nature Communications\final_submission\manuscript_source_data';
load(fullfile(dataFolder,'cortex_anatomy','anatomy_cortical_pairwise_distance_distribution.mat'));
area = B;
dMids = 0.5*(rValues(2:end)+rValues(1:end-1));
cumN = interp1(rValues,area,[dValues,dValues(end)+1],'spline','extrap');
diffN = mean(diff(cumN),2)*200000/N2;
diffN(end) = nan;

Pxy_D2 = Pxy_D./count_D'.* exp(-dValues.^2./6);

Pxy2 = nansum(Pxy_D2.*diffN',2);

Pxx2 = Pxx/N;

dS = [0,10.^linspace(-3,-1.5,6)];

low_noise = (16/1e3)^2;;
lam = 8.5;

figureNB(7,7);
    % plotwitherror(freq(freq<300),pre_CAF(freq<300,:),'M','LineWidth',1,'color','b');
    % plotwitherror(freq(freq<300),post_CAF(freq<300,:),'M','LineWidth',1,'color','r');
    plot(f,f*0+low_noise,'--k','LineWidth',1);
    hold on;
    for i = 1:length(dS)
        S = dS(i);
        B = exp(-(2*pi*f(:)*S).^2/2);
        if(i==1)
            plot(f,lam*N2*Pxx2 + 0.2*lam*N2*(N2-1)*B.*Pxy2,'color',[0,0,0],'LineWidth',1)
        else
            plot(f,lam*N2*Pxx2 + 0.2*lam*N2*(N2-1)*B.*Pxy2,'color',[0.6,0.6,0.6],'LineWidth',1)
        end
    end
    B = exp(-(2*pi*f(:)*8e-3*sqrt(2)).^2/2);
    % plot(f,lam*N2*Pxx2 + 0.2*lam*N2*(N2-1)*B.*Pxy2,'color','k','LineWidth',1)


    xlim([10,1e3])
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlabel('Frequency (Hz)')
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    gcaformat;