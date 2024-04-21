addpath 'C:\Users\brake\Documents\GitHub\bAP_electric_fields\simulations\calculate_unitary_bAP'
% Load precomputed cross spectra
load('C:\Users\brake\Desktop\unitaryAPNew.mat')
savedUnitaryAP = permute(savedUnitaryAP,[2,1,3]);

% Calculate new coordinates uiAi
[sa,X] = network_simulation_beluga.getHeadModel;
electrode = find(strcmp(sa.clab_electrodes,'Cz'));

% Restrict simulation to area under POz
% electrode = find(strcmp(sa.clab_electrodes,'C2'));
% x0 = sa.locs_3D_orig(electrode,:);
% d = vecnorm(X.vertices-x0,2,2);
% idcs = find(d<50);
% X.vertices = X.vertices(idcs,:);
% sa.cortex75K.V_fem = sa.cortex75K.V_fem(:,idcs,:);
% sa.cortex75K.normals = sa.cortex75K.normals(idcs,:);

M = size(X.vertices,1);
Lxyz = zeros(M,3);
for idx = 1:M
    L0 = squeeze(sa.cortex75K.V_fem(electrode,idx,:))'; % czIDX = 49
    vz = sa.cortex75K.normals(idx,:);
    [vx,vy] = getOrthBasis(vz);
    A = [vx(:),vy(:),vz(:)];
    Lxyz(idx,:) = 1e-6*L0*A;
end
% Set up a KD tree to efficiently search for nearby vertices
Mdl = KDTreeSearcher(X.vertices);
[f0,S] = import_paralytic_data;

% Result variables
L = 16e3;
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

R = 0.2;
lam = 10;
dF = 5*2.^[0:5];
sig = 8e-3*sqrt(2);
% N2 = 1151350505; % caluclated from total surface area times 70000
N2 = 16e9;

low_noise = 1e-3;
blue = [17,82,185]/255;

figureNB(15,7.6);
for i = 1:length(dF)
    subplot(2,3,i);
    plot(f0,S(:,1:8),'k','Linewidth',1);
    hold on;

    plot([1,3e3],low_noise*[1,1],'r','LineWidth',1.5,'LineStyle','-')

    F0 = dF(i);
    B = exp(-2*(pi*(f(:)-F0)*sig).^2);
    plt(i) = plot(f,lam*N2*Pxx*nan + R*lam*N2*(N2-1)*B.*Pxy*nan,'color',blue,'LineWidth',1);
        % plot(f,lam*N2*Rxx+R*lam*N2*(N2-1)*Rxy,'color',[0.6,0.6,0.6],'LineWidth',1,'LineStyle',':')
        % plot(f,lam*N2*Rxx+0*R*lam*N2*(N2-1)*Rxy,'color',[0.6,0.6,0.6],'LineWidth',1,'LineStyle',':')

    set(gca,'xscale','log')
    set(gca,'yscale','log');
    xlabel('Frequency (Hz)');
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    gcaformat;
    set(gca,'FontSize',8)

    xticks([1,10,100])
    xticklabels([1,10,100])
    xlim([1,400])

    ylim([1e-6,1e2])
    yticks([1e-6,1e-2,1e2])
end


badIdcs = find(isnan(squeeze(savedUnitaryAP(1,1,:))))';

while true
    N = N+1;

    % Randomly sample unitary AP profiles
    i = neuronIdcs(iSampling,1);
    j = neuronIdcs(iSampling,2);

    if(sum(badIdcs==i) || sum(badIdcs==j))
        N=N-1;
        iSampling = iSampling+1;
        continue;
    end

    % Randomly sample vertex from cortex
    iX = randi(M);
    jX = randi(M);
    % Sample second vertex with higher probability for nearby neurons
    % d = exprnd(6);
    % idcs = rangesearch(Mdl,X.vertices(iX,:),d);
    % if(length(idcs{1})>1)
    %     jX = randsample(idcs{1},1);
    % else
    %     jX = idcs{1};
    % end
    % d = norm(X.vertices(iX,:)-X.vertices(jX,:));
    uiAi = Lxyz(iX,:);
    ujAj = Lxyz(jX,:);

    % Compute eeg cross spectrum
    eeg1 = uiAi*savedUnitaryAP(:,:,i);
    eeg2 = ujAj*savedUnitaryAP(:,:,j);
    R_eeg = xcorr(eeg2,eeg1,4e3,'unbiased');
    R_eeg = detrend(R_eeg(2:end),'constant');
    S12 = fft(R_eeg)/fs^2*L; % units = per spike rather than per Hztime
    Sxy = 2*S12(2:4001).*exp(2*pi*sqrt(-1)*f*3999/fs);
    % Scale by distance (spike time correlation)
    Sxy = real(Sxy(:));
    % Sxy(Sxy<=0) = 1e-30;
    % Sxy = log10(Sxy);

    % Compute eeg unitary spectrum
    R_eeg = xcorr(eeg1,eeg1,4e3,'unbiased');
    R_eeg = detrend(R_eeg(2:end),'constant');
    S11 = fft(R_eeg)/fs^2*L; % units = per spike rather than per Hztime
    Sxx = 2*S11(2:4001).*exp(2*pi*sqrt(-1)*f*3999/fs);
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
        for i = 1:length(dF)
            Rxx = Pxx/N;
            Rxy = Pxy/N;
            F0 = dF(i);
            B = exp(-2*(pi*(f(:)-F0)*sig).^2);
            plt(i).YData = lam*N2*Rxx + R*lam*N2*(N2-1)*B.*Rxy;
        end
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




Rxx = Pxx/N;
Rxy = Pxy/N;




% N2 = 16e9;
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

