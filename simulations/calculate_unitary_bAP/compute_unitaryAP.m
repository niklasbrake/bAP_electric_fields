warning('off','signal:findpeaks:largeMinPeakHeight');
% [sa,X] = network_simulation_beluga.getHeadModel;

% folder = 'E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\neuron_models';
folder = '/lustre04/scratch/nbrake/data/simulations/unitary_AP';

% Calculate spike triggered average
% [mtype,ei_type,layer,morph] = getMtypes(folder);
% [savedUnitaryAP,N,EI,files] = getUnitaryAP(folder);
% save('/lustre04/scratch/nbrake/data/simulation_analyzed/unitaryAP/unitaryAP.mat');
% return;
%{
save('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\unitaryAP_all.mat','savedUnitaryAP','mtype','ei_type','layer','morph');
%}

% Compute unitary spectrum
% load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\unitaryAP_all.mat')
load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\unitaryAPNew.mat')
psd_unit = getUnitarySpectrum(savedUnitaryAP)
% save('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\unitary_AP_PSD.mat','freq','psd_unit','mtype');

% Bootstrap average

bs_sample = reshape(sample_blue_neurons(1e6),1e3,1e3);
for i = 1:1e3
    waitbar(i/1e3)
    bootstrap_Est(:,i) = nanmean(psd_unit(:,bs_sample(:,i)),2);
end


avgEnergy = sum(bootstrap_Est*mean(diff(freq)));

mean(avgEnergy)
std(avgEnergy)

figureNB(6,5)
    plotwitherror(freq,bootstrap_Est,'SE','color','k','LineWidth',1)
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    ylabel(['Unitary AP energy density (' char(956) 'V^2/Hz)'])
    xlabel('Frequency (Hz)')
    xlim([1,2e3])
    xticks([1,10,100,1000])
    xticklabels([1,10,100,1000])
    gcaformat;




function [mtype,ei_type,layer,morph] = getMtypes(folder)
    F = dir(folder);
    F = F(3:end);
    M = length(F);

    for i = 1:M
        items = split(F(i).name,'_');
        layer{i} = items{1};
        morph{i} = items{2};
        if(strcmp(items{3},'L1') || strcmp(items{3},'L4'))
            morph{i} = [items{2} '_' items{3}];
        end
        mtype{i} = [layer{i} '_' morph{i}];
    end
    e_mtypes = {'PC','SS','SP','STPC','UTPC','TTPC1','TTPC2','TPC','BPC','IPC','TPC_L1','TPC_L4'};
    ei_type = ismember(morph,e_mtypes);
end
function [savedUnitaryAP,N,EI,files] = getUnitaryAP(folder)

    F = dir(folder);
    F = F(3:end);
    M = length(F);

    savedUnitaryAP = nan(2001,3,M);
    N = nan(M,1);

    fs = 16e3; % Hz

    h = waitbar(0);
    files = cell(M,1);
    tic
    for i = 1:M
        update_waitbar(h,i,M);
        files{i} = F(i).name;
        folder0 = fullfile(folder,F(i).name,'matlab_recordings');

        fid = fopen(fullfile(folder,F(i).name,'EI_ratio.csv'));
        ei = textscan(fid,'%s');
        fclose(fid);
        ei = ei{1}{1};
        if(str2num(ei)==1 || str2num(ei)==5)
            ei = ['0' ei];
        end

        try
            load(fullfile(folder0,sprintf('synaptic_input_EI%s.mat',ei)));
        catch
            disp(F(i).name)
            continue;
        end
        [y,x] = findpeaks(voltage,'MinPeakHeight',0);
        N(i) = length(x);
        EI(i) = str2num(ei);
        unitaryAP = zeros(2001,3,length(x));
        for j = 1:length(x)
            idcs = max(min(x(j)-1e3:x(j)+1e3,length(voltage)),1);
            y = dipoles(idcs,:);
            y(idcs==1,:) = 0;
            y(idcs==length(voltage),:) = 0;
            unitaryAP(:,:,j) = y;
        end
        unitaryAP = unitaryAP-nanmedian(unitaryAP);
        savedUnitaryAP(:,:,i) = nanmedian(unitaryAP,3);
    end
    % updateEIratio(folder,files,N,EI);
end

% Change EI ratio to get spiking frequency between 0 and 100 Hz
function updateEIratio(folder,files,N,EI)

    idcs = find(N < 5 | N>100)
    for i = 1:length(idcs)
        j = idcs(i);
        if(N(j) <= 5)
            newEI = EI(j)/2;
        else
            newEI = EI(j)*2;
        end

        fprintf('%s, %.2f, EI: %f -> %f\n',files{j},N(j),EI(j),newEI)

        csvwrite(fullfile(folder,files{j},'EI_ratio.csv'),newEI);
    end

    fid = fopen('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\changedEI.txt','w');
    fprintf(fid,'%s\n',files{idcs});
    fclose(fid)
end

% Compute unitary spectrum
function psd_unit = getUnitarySpectrum(savedUnitaryAP)
    [sa,X] = network_simulation_beluga.getHeadModel;
    fs = 16e3;
    idcs = sa.cortex2K.in_from_cortex75K;
    m = size(savedUnitaryAP,3);
    N = length(idcs);
    n = size(savedUnitaryAP,1);
    psd_unit = zeros(1e3,m);
    freq = fs/n:fs/n:fs/2;
    h = waitbar(0);
    for i = 1:N
        update_waitbar(h,i,N);
        eeg = network_simulation_beluga.getEEG(savedUnitaryAP,sa,idcs(i));
        xdft = fft(eeg);
        xdft = xdft(2:n/2+1,:);
        psdx = (1/(fs*n)) * abs(xdft).^2;
        psdx(1:end-1,:) = 2*psdx(1:end-1,:);
        psd_unit = psd_unit+psdx/N;
    end
end

function [Sii,Sij] = compute_cross_spectra(savedUnitaryAP)
    m = size(savedUnitaryAP,3);
    N = m*(m-1)/2;
    Sij = zeros(3,3,100,N);
    Sii = zeros(3,3,100,m);

    fs = 16e3;
    L = 4e3;
    f0 = fs/L:fs/L:fs/2;
    f = 10.^linspace(0,log10(3e3),100);

    h = waitbar(0);
    count = 1;
    for i=1:m
        for x = 1:3
            for y = 1:3
                R_eeg = xcorr(savedUnitaryAP(:,x,i),savedUnitaryAP(:,y,i));
                R_eeg = R_eeg(2:end);
                S12 = fft(R_eeg)/fs;
                S12 = 2*S12(2:2001).*exp(2*pi*sqrt(-1)*f0(:)*1999/fs);
                Sii(x,y,:,i) = interp1(f0,real(S12),f,'linear','extrap');
            end
        end
        for j = i+1:m
            for x = 1:3
                for y = 1:3
                    R_eeg = xcorr(savedUnitaryAP(:,x,i),savedUnitaryAP(:,y,j));
                    R_eeg = R_eeg(2:end);
                    S12 = fft(R_eeg)/fs;
                    S12 = 2*S12(2:2001).*exp(2*pi*sqrt(-1)*f0(:)*1999/fs);
                    Sij(x,y,:,count) = interp1(f0,real(S12),f,'linear','extrap');
                end
            end
            count = count + 1;
        end
        update_waitbar(h,count,N);
    end
end