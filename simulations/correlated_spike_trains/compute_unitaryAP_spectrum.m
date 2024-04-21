load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\unitaryAPNew.mat')

savedUnitaryAP = savedUnitaryAP-nanmedian(savedUnitaryAP(1:1500,:,:));
savedUnitaryAP = gpuArray(savedUnitaryAP(2:end,:,:));

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
Lxyz = Lxyz';
Lxyz = gpuArray(Lxyz(:,sa.cortex10K.in_from_cortex75K));

% Compute unitary spectrum
fs = 16e3;
m = size(savedUnitaryAP,3);
psd_unit = zeros(1e3,m);
n = size(savedUnitaryAP,1);
freq = fs/n:fs/n:fs/2;

psd = gpuArray(zeros(n/2,m));
h = waitbar(0);
for j = 1:m
    update_waitbar(h,j,m);
    eeg1 = savedUnitaryAP(:,:,j)*Lxyz;
    xdft = fft(eeg1);
    xdft = xdft(2:n/2+1,:);
    psd(:,j) = 2*mean(abs(xdft/fs).^2,2);
end
