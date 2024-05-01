warning('off','signal:findpeaks:largeMinPeakHeight');

matObj = matfile('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\unitary_AP_EI_ratio\EI_ratio.mat');
cellIDs = who(matObj);

load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\unitaryAP.mat')
files = cellfun(@(x)strrep(x,'-','_'),files,'UniformOutput',false);
allUnitaryAP = savedUnitaryAP;
getUAP = @(id) allUnitaryAP(2:end,3,find(strcmp(id,files)));


EI_vec = {'01','1.5','2.1','3.1','4.5','6.6','9.7','14.1','20.6','30'};

M = length(cellIDs);
N = zeros(length(EI_vec),1);
B2 = zeros(length(EI_vec),M);
apCount = zeros(length(EI_vec),M);


fs = 16e3; % Hz
h = waitbar(0);
for ii = 1:M
    update_waitbar(h,ii,M);
    meanAP = getUAP(cellIDs{ii});
    cell = matObj.(cellIDs{ii});

    for k = 1:length(EI_vec)
        [f0,~,psd] = eegfft(cell.passive.time*1e-3,detrend(cell.passive.dipoles(:,3,k)),0.5,0.4);
        psd_passive(:,k,ii) = mean(psd,2);

        [f0,~,psd] = eegfft(cell.active.time*1e-3,detrend(cell.active.dipoles(:,3,k)),0.5,0.4);
        psd_active(:,k,ii) = mean(psd,2);

        [y,x] = findpeaks(cell.active.voltage(:,k),'MinPeakHeight',0);
        N(k) = length(x);
    end

    n = length(meanAP);
    xdft = fft(meanAP);
    xdft = xdft(1:n/2+1);
    psdx = (1/(fs*n)) * abs(xdft).^2;
    psdx(2:end-1) = 2*psdx(2:end-1);
    freq = (0:fs/n:fs/2)';

    psd_unit = interp1(freq,psdx,f0,'linear','extrap');
    psd_unit = psd_unit*length(meanAP)/fs/1; % per second

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

firingFrequency = apCount(:)/range(cell.active.time)*1e3;

FT = polyfit(firingFrequency(:),R(:),1);
t = linspace(0.1,100,20);

blue = [17,82,185]/255;
red = clrsPT.qualitative_CM.red;
figureNB(9,4.4);
subplot(1,2,1);
    plot(firingFrequency,B2(:),'.k')
    line([0.1,100],[0.1,100],'color',red,'LineWidth',1)
    xlabel('Firing rate (Hz)')
    ylabel('Scaling factor (\beta)')
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlim([1,100])
    xticks([0.1,1,10,100])
    xticklabels([0.1,1,10,100])
    gcaformat;
subplot(1,2,2);
    plot(firingFrequency,R(:),'.k')
    hold on;
    plot(t,polyval(FT,t),'color',red,'LineWidth',2);
    ylim([0,1])
    xlabel('Firing rate (Hz)')
    ylabel('R^2')
    set(gca,'xscale','log')
    xlim([1,100])
    xticks([0.1,1,10,100])
    xticklabels([0.1,1,10,100])
gcaformat(gcf,true,8);

edges = [0,20,40,80,Inf]
bins = discretize(firingFrequency,edges);
str_edges = {'0-20','20-40','40-80','>80'};

split_Y = splitapply(@(x)mean(x,2),psd_Y,bins');
split_Yhat = splitapply(@(x)mean(x,2),psd_Yhat,bins');



figureNB(5,4);
axes('Position',[0.16, 0.62, 0.34, 0.32])
 plot(f0,split_Y(:,1),'color','k','LineWidth',2);
    hold on;
    plot(f0,split_Yhat(:,1),'color',blue,'LineWidth',1);
    xlim([1,3e3]);
    ylim([1e-4,1]);
    xticks([1,10,100,1000]);
    xticklabels({});
    set(gca,'xscale','log')
    set(gca,'yscale','log')

axes('Position',[0.6, 0.62, 0.34, 0.32])
 plot(f0,split_Y(:,2),'color','k','LineWidth',2);
    hold on;
    plot(f0,split_Yhat(:,2),'color',blue,'LineWidth',1);
    xlim([1,3e3]);
    ylim([1e-4,1]);
    xticks([1,10,100,1000]);
    xticklabels({});
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    yticklabels({});

axes('Position',[0.16, 0.21, 0.34, 0.32])
 plot(f0,split_Y(:,3),'color','k','LineWidth',2);
    hold on;
    plot(f0,split_Yhat(:,3),'color',blue,'LineWidth',1);
    xlim([1,3e3]);
    ylim([1e-4,1]);
    xticks([1,10,100,1000]);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlabel('Frequency (Hz)');

axes('Position',[0.6, 0.21, 0.34, 0.32])
    plot(f0,split_Y(:,4),'color','k','LineWidth',2);
    hold on;
    plot(f0,split_Yhat(:,4),'color',blue,'LineWidth',1);
    xlim([1,3e3]);
    ylim([1e-4,1]);
    xticks([1,10,100,1000]);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlabel('Frequency (Hz)');
    yticklabels({});

gcaformat(gcf,true,8);