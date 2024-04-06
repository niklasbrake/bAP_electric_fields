warning('off','signal:findpeaks:largeMinPeakHeight');

load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\unitaryAPNew.mat')
allUnitaryAP = savedUnitaryAP;
folder = 'E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\neuron_models';
F = dir(folder);
F = F(3:end);
cellIDs = {F(:).name};
getUAP = @(id) allUnitaryAP(2:end,3,find(strcmp(id,cellIDs)));

folder = 'E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\sampled_mtypes';
F = dir(folder);
F = F(3:end);


EI_vec = {'01','1.5','2.1','3.1','4.5','6.6','9.7','14.1','20.6','30'};
iVec = 1:length(F);
clrs = flipud(jet(10));

% M = length(F);
M = length(iVec);
h = waitbar(0);
% psd_active = zeros(4096,length(EI_vec));
% psd_passive = zeros(4096,length(EI_vec));
N = zeros(length(EI_vec),1);
B2 = zeros(length(EI_vec),M);
apCount = zeros(length(EI_vec),M);

fs = 16e3; % Hz
h = waitbar(0);
for ii = 1:M
    update_waitbar(h,ii,M);
    i = iVec(ii);
    folder0 = fullfile(folder,F(i).name,'matlab_recordings');

    for k = 1:length(EI_vec)
        load(fullfile(folder0,sprintf('synaptic_input_EI%s_passive.mat',EI_vec{k})))
        [f0,~,psd] = eegfft(time*1e-3,detrend(dipoles(:,3)),0.5,0.4);
        psd_passive(:,k,ii) = mean(psd,2);

        load(fullfile(folder0,sprintf('synaptic_input_EI%s.mat',EI_vec{k})))
        [f0,~,psd] = eegfft(time*1e-3,detrend(dipoles(:,3)),0.5,0.4);
        psd_active(:,k,ii) = mean(psd,2);

        [y,x] = findpeaks(voltage,'MinPeakHeight',0);
        N(k) = length(x);
        % unitaryAP = zeros(1,2001);
        % for j = 1:length(x)
        %     idcs = max(min(x(j)-1e3:x(j)+1e3,length(dipoles)),1);
        %     y = dipoles(idcs,3);
        %     y(idcs==1) = 0;
        %     y(idcs==length(dipoles)) = 0;
        %     unitaryAP(j,:) = y;
        % end
        % unitaryAP = unitaryAP-median(unitaryAP,2);
        % meanAP = [meanAP;unitaryAP];
    end

    meanAP = getUAP(F(i).name);

    % meanAP = meanAP';
    % meanAP(isnan(meanAP(:))) = 0;
    % meanAP = mean(meanAP,2);

    n = length(meanAP);
    xdft = fft(meanAP);
    xdft = xdft(1:n/2+1);
    psdx = (1/(fs*n)) * abs(xdft).^2;
    psdx(2:end-1) = 2*psdx(2:end-1);
    freq = (0:fs/n:fs/2)';

    psd_unit = interp1(freq,psdx,f0,'linear','extrap');
    psd_unit = psd_unit*length(meanAP)/fs/1; % per second
    % psd_unit = psdx*length(meanAP)/fs/1; % per second

    % psd_passive = interp1(f0,psd_passive,freq,'linear','extrap');
    % psd_active = interp1(f0,psd_active,freq,'linear','extrap');

    for k = 1:10
        X0 = psd_passive(:,k,ii);
        X1 = psd_unit;
        Y = psd_active(:,k,ii);
        idcs = find(and(and(Y>1e-19,X1>1e-19),f0<1e3));
        B2(k,ii) = fminbnd(@(B1) neglnlike(B1,X0(idcs),X1(idcs),Y(idcs)),0,N(k)*2);
        apCount(k,ii) = N(k);

        Yhat = X0+B2(k,ii)*X1;
        psd_Y(:,10*(ii-1)+k) = Y;
        psd_Yhat(:,10*(ii-1)+k) = Yhat;

        Y = log(Y);
        Yhat = log(Yhat);
        R(k,ii) = 1 - sum((Y-Yhat).^2)./sum((Y-mean(Y)).^2);
    end

end
delete(h)

firingFrequency = apCount(:)/range(time)*1e3;

FT = polyfit(firingFrequency(:),R(:),1);
t = linspace(0.1,100,20);

% edges = linspace(0,max(firingFrequency),8);
edges = [0,1,2,4,10,20,40,80,Inf]
% edges = [0,quantile(firingFrequency(firingFrequency>0),linspace(0,1,5))];
bins = discretize(firingFrequency,edges);
str_edges = {'0-1','1-2','2-4','4-10','10-20','20-40','40-80','>80'};

clrs = clrsPT.sequential(length(edges)+2);
clrs = clrs(3:end,:);

figureNB;
subplot(1,3,1);
    plot(firingFrequency,B2(:),'.k')
    line([0.1,100],[0.1,100],'color','r','LineWidth',1)
    xlabel('Firing frequency (Hz)')
    ylabel('\beta_1')
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlim([1,100])
    xticks([0.1,1,10,100])
    gcaformat;
subplot(1,3,2);
    plot(firingFrequency,R(:),'.k')
    hold on;
    plot(t,polyval(FT,t),'color','r','LineWidth',2);
    ylim([0,1])
    xlabel('Firing frequency (Hz)')
    ylabel('R^2')
    set(gca,'xscale','log')
    xlim([1,100])
    xticks([0.1,1,10,100])
    gcaformat;
subplot(1,3,3);
    plot(splitapply(@(x)mean(x,2),psd_Y,bins'));
    set(gca,'ColorOrder',clrs)
    xlim([1,3e3]);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlabel('Frequency (Hz)');
    ylabel('PSD');
    gcaformat;
    colormap(clrs)
    C = colorbar;
    C.Ticks = (1/16):1/8:(1-1/16);


split_Y = splitapply(@(x)mean(x,2),psd_Y,bins');
split_Yhat = splitapply(@(x)mean(x,2),psd_Yhat,bins');


blue = [17,82,185]/255;
red = clrsPT.qualitative_CM.red;
figureNB(9,8);
subplot(2,2,1);
    plot(firingFrequency,B2(:),'.k')
    line([0.1,100],[0.1,100],'color',red,'LineWidth',1)
    xlabel('Firing frequency (Hz)')
    ylabel('\beta_1')
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlim([1,100])
    xticks([0.1,1,10,100])
    gcaformat;
subplot(2,2,2);
    plot(firingFrequency,R(:),'.k')
    hold on;
    plot(t,polyval(FT,t),'color',red,'LineWidth',2);
    ylim([0,1])
    xlabel('Firing frequency (Hz)')
    ylabel('R^2')
    set(gca,'xscale','log')
    xlim([1,100])
    xticks([0.1,1,10,100])
    gcaformat;
for i = 1:8
    subplot(4,4,8+i)
    plot(f0,split_Y(:,i),'color','k','LineWidth',2);
    hold on;
    plot(f0,split_Yhat(:,i),'color',blue,'LineWidth',1);
    xlim([1,3e3]);
    ylim([1e-4,1]);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    if(i>4)
        xlabel('Frequency (Hz)');
    else
        xticklabels({});
    end
    if(i==1 | i == 5)
        ylabel('PSD');
    else
        yticklabels({});
    end
    title([str_edges{i} ' Hz']);
    gcaformat;
end