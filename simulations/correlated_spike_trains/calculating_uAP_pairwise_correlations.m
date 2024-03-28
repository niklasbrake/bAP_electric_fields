
load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\unitaryAP_all.mat');

N = m*(m-1)/2;
m = size(savedUnitaryAP,3);
Sij = zeros(3,3,100,N);
Sii = zeros(3,3,100,m);

fs = 16e3;
L = 4e3;
f0 = fs/L:fs/L:fs/2;
f = 10.^linspace(0,log10(3e3),100);

h = waitbar(0); tic;
count = 1;
for i=1:m
    for x = 1:3
        for y = 1:3
            R_eeg = xcorr(savedUnitaryAP(:,x,i),savedUnitaryAP(:,y,i));
            R_eeg = detrend(R_eeg(2:end),'constant');
            S12 = fft(R_eeg)/fs;
            S12 = 2*S12(2:2001).*exp(2*pi*sqrt(-1)*f0(:)*1999/fs);
            Sii(x,y,:,i) = interp1(f0,real(S12),f,'linear','extrap');
        end
    end

    for j = i+1:m
        for x = 1:3
            for y = 1:3
                R_eeg = xcorr(savedUnitaryAP(:,x,i),savedUnitaryAP(:,y,j));
                R_eeg = detrend(R_eeg(2:end),'constant');
                S12 = fft(R_eeg)/fs;
                S12 = 2*S12(2:2001).*exp(2*pi*sqrt(-1)*f0(:)*1999/fs);
                Sij(x,y,:,count) = interp1(f0,real(S12),f,'linear','extrap');
            end
        end
        count = count + 1;
    end
    update_waitbar(h,count,N);
end