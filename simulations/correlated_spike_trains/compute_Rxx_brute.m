load('/lustre04/scratch/nbrake/data/simulation_analyzed/unitaryAP/unitaryAP.mat');
savedUnitaryAP = savedUnitaryAP(2:end,:,:);

load('/lustre04/scratch/nbrake/data/anatomy_nyhead_model.mat');
M = size(sa.cortex75K.vc,1);
Lxyz = zeros(M,3);
for idx = 1:M
    L0 = squeeze(sa.cortex75K.V_fem(49,idx,:))'; % czIDX = 49
    vz = sa.cortex75K.normals(idx,:);
    [vx,vy] = getOrthBasis(vz);
    A = [vx(:),vy(:),vz(:)];
    Lxyz(idx,:) = 1e-6*L0*A;
end
Lxyz = Lxyz';

% Compute unitary spectrum
fs = 16e3;
m = size(savedUnitaryAP,3);
psd_unit = zeros(1e3,m);
n = size(savedUnitaryAP,1);
freq = fs/n:fs/n:fs/2;

psd = zeros(n/2,m);
parfor j = 1:m
    eeg1 = savedUnitaryAP(:,:,j)*Lxyz;
    xdft = fft(eeg1);
    xdft = xdft(2:n/2+1,:);
    psd(:,j) = 2*mean(abs(xdft/fs).^2,2);
end

save('/lustre04/scratch/nbrake/data/simulation_analyzed/unitaryAP/unitarySpectrum.mat','freq','psd')

function [x1,x2] = getOrthBasis(x0)
    x0 = x0/norm(x0,2);
    if(prod(size(x0))==3)
        if(x0>0.9)
            x1 = [0,1,0];
        else
            x1 = [1,0,0];
        end
        x1 = x1-x0*sum(x1.*x0);
        x1 = x1/norm(x1,2);
        x2 = cross(x1,x0);
    else
        idcs = x0(:,1)>0.9;
        x1 = zeros(size(x0));
        x1(idcs,2) = 1;
        x1(~idcs,1) = 1;

        x1 = x1-x0.*sum(x1.*x0,2);
        x1 = x1./vecnorm(x1,2,2);
        x2 = cross(x1,x0);
    end
end