warning('off','signal:findpeaks:largeMinPeakHeight');
% [sa,X] = network_simulation_beluga.getHeadModel;

folder = 'E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\neuron_models';

% Calculate spike triggered average
%{
[mtype,ei_type,layer,morph] = getMtypes(folder);
[savedUnitaryAP,N,EI,files] = getUnitaryAP(folder);
save('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\unitaryAP_all.mat','savedUnitaryAP','mtype','ei_type','layer','morph');
%}

% Compute unitary spectrum
%{
load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\unitaryAP_all.mat')
psd_unit = getUnitarySpectrum(savedUnitaryAP)
save('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\unitary_AP_PSD.mat','freq','psd_unit','mtype');
%}


% Bootstrap average
load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\mtype_abundance.mat','mtype_abundance');
[J,ID] = findgroups(mtype);

for i = 1:length(ID)
    abundance(i) = mtype_abundance(ID{i},:).Abundance;
end

mtype_abundance.count = splitapply(@(x) length(x),pEst,J(:)')';
scaledAbundance = mtype_abundance.Abundance./mtype_abundance.count;
prop = scaledAbundance(J);
F = cumsum(prop)/sum(prop);
[~,idcs] = unique(F);
m = length(J);
R = rand(m,1e3);
for i = 1:1e3
    waitbar(i/1e3)
    bs_sample = interp1(F(idcs),idcs,R(:,i),'next','extrap');
    bootstrap_Est(:,i) = nanmean(psd_unit(:,bs_sample),2);
end


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
    savedUnitaryAP = nan(2001,3,M);
    N = nan(M,1);

    fs = 16e3; % Hz

    h = waitbar(0);
    files = cell(M,1);
    tic
    for i = 1:M
        waitbar(i/M,h,round(toc/i*(M-i)/60,2))
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

    idcs = find(N==0 | N>100)
    for i = 1:length(idcs)
        j = idcs(i);
        if(N(j)==0)
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
    psd_unit = zeros(8001,m);
    N = length(idcs);
    n = fs;
    freq = fs/n:fs/n:fs/2;
    for i = 1:N
        waitbar(i/N)
        eeg = network_simulation_beluga.getEEG(savedUnitaryAP,sa,idcs(i));
        paddedAP = [zeros(7e3,m);eeg;zeros(7e3-1,m)];
        xdft = fft(paddedAP);
        xdft = xdft(2:n/2+1,:);
        psdx = (1/(fs*n)) * abs(xdft).^2;
        psdx(1:end-1,:) = 2*psdx(1:end-1,:);

        psd_unit = psd_unit+psdx/N;
    end
end