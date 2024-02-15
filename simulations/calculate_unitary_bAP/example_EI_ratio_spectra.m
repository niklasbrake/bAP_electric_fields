% function example_EI_ratio_spectra
% [sa,X] = network_simulation_beluga.getHeadModel;

folder = 'E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\sampled_mtypes';
F = dir(folder);
F = F(3:end);

EI_vec = {'01','1.5','2.1','3.1','4.5','6.6','9.7','14.1','20.6','30'};
iCtx = 16e3;

% i = 46;
% i=50;
% 21    35    59
i=59;
% i = 48;
folder0 = fullfile(folder,F(i).name,'matlab_recordings');
clrs = flipud(jet(10));

N = zeros(length(EI_vec),1);
psd_active = [];
psd_passive = [];
figureNB;
meanAP = [];
for k = 1:length(EI_vec)
    load(fullfile(folder0,sprintf('synaptic_input_EI%s_passive.mat',EI_vec{k})))
    eeg = network_simulation_beluga.getEEG(dipoles,sa,iCtx);
    [f,~,psd] = eegfft(time*1e-3,detrend(eeg),0.5,0.4);
    psd_passive(:,k) = mean(psd,2);
    % Y = detrend(eeg(32001:end),'constant');
    % [psd_passive(:,k),f] = pmtm(Y,2,[],16e3);

    load(fullfile(folder0,sprintf('synaptic_input_EI%s.mat',EI_vec{k})))
    eeg = network_simulation_beluga.getEEG(dipoles,sa,iCtx);
    % Y = detrend(eeg(32001:end),'constant');
    % [psd_active(:,k),f] = pmtm(Y,2,[],16e3);
    % [freq,~,psd] = eegfft(time*1e-3,eeg,0.5,0.4);
    [f,~,psd] = eegfft(time*1e-3,detrend(eeg),0.5,0.4);
    % psd_active(:,k) = interp1(freq,mean(psd,2),f,'linear','extrap');
    psd_active(:,k) = mean(psd,2);
    [y,x] = findpeaks(voltage,'MinPeakHeight',0);
    x = x(x>0);
    N(k) = sum(x>0);
    if(N(k)==0)
        continue;
    end
    [~,I] = min(N(N>0));
    unitaryAP = zeros(1,2001);
    for j = 1:length(x)
        idcs = max(min(x(j)-1e3:x(j)+1e3,length(eeg)),1);
        y = eeg(idcs);
        y(idcs==1) = 0;
        y(idcs==length(eeg)) = 0;
        unitaryAP(j,:) = y;
    end
    unitaryAP = unitaryAP-median(unitaryAP,2);
    meanAP = [meanAP;unitaryAP];
end

firingFrequency = N(:)/range(time)*1e3;

m = size(meanAP,1);

paddedAP = [zeros(m,7e3),meanAP,zeros(m,7e3-1)]';
paddedAP(isnan(paddedAP(:))) = 0;
paddedAP = median(paddedAP,2);

n = length(paddedAP);
xdft = fft(paddedAP);
xdft = xdft(1:n/2+1);
psdx = (1/(fs*n)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:fs/n:fs/2;

psd_unit = interp1(freq,psdx,f,'linear','extrap');
psd_unit = psd_unit*length(paddedAP)/fs/1; % per second

k=2;
figureNB;
subplot(2,2,1);
    plot((-1000:1000)/16,meanAP','color',[0.5,0.5,0.5]*1.3)
    hold on;
    plot((-1000:1000)/16,nanmean(meanAP),'b','LineWidth',1);
    xlim([-15,15])
    ylabel(['Unitary AP (' char(956) 'V)'])
    xlabel('Time (ms)')
subplot(2,2,2);
    plot(f,psd_unit,'color',[0.5,0.5,0.5]*1.3);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    hold on;
    plot(f,mean(psd_unit,2),'color','b','LineWidth',2)
    xlim([1,5e3])
    ylim([1e-22,1e-16])
    xticks([1,10,100,1000])
    xticklabels([1,10,100,1000])
    xlabel('Frequency (Hz)')
    ylabel(['Power (' char(956) 'V^2/Hz)'])
subplot(2,2,3);
    plot(f,psd_active(:,k),'color','k')
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    hold on;
    xlim([1,5e3])
    % ylim([1e-22,1e-16])
    xticks([1,10,100,1000])
    xticklabels([1,10,100,1000])
    xlabel('Frequency (Hz)')
    ylabel(['Power (' char(956) 'V^2/Hz)'])
subplot(2,2,4);
    plot(f,psd_active(:,k),'color','k')
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    hold on;
    plot(f,psd_passive(:,k),'color',[0.6,0.6,0.6])
    plot(f,firingFrequency(k)*psd_unit,'color','b','LineWidth',2)
    xlim([1,5e3])
    % ylim([1e-22,1e-16])
    xticks([1,10,100,1000])
    xticklabels([1,10,100,1000])
    xlabel('Frequency (Hz)')
    ylabel(['Power (' char(956) 'V^2/Hz)'])


B2 = zeros(10,1);
X0 = mean(psd_passive,2);
figureNB(24,9);
for k = 1:10
    if(N(k)>0)
        % X1 = psd(:,k);
        X1 = psd_unit;
        Y = psd_active(:,k);
        idcs = find(and(and(Y>1e-19,X1>1e-19),f<3e3));
        B2(k) = fminbnd(@(B1) neglnlike(B1,X0(idcs),X1(idcs),Y(idcs)),0,N(k)*2);
    else
        B2(k) = 0;
    end
    % B2(k) = N(k)/2;
    subplot(2,5,k)
    plot(psd_active(:,k),'color','k')
    hold on;
    plot(psd_passive(:,k),'color',[0.6,0.6,0.6])
    plot(B2(k)*X1,'color','b','LineWidth',2)
    set(gca,'xscale','log')
    xlim([1,2e3])
    % ylim([1e-20,1e-14]);
    set(gca,'yscale','log')
    title(sprintf('%d Hz',round(firingFrequency(k))))
    xticks([1,10,100,1000])
    xticklabels([1,10,100,1000])
    xlabel('Frequency (Hz)')
    ylabel(['Power (' char(956) 'V^2/Hz)'])
    gcaformat;
end

m = ceil(max([firingFrequency;B2]));
figureNB;
    plot(firingFrequency,B2,'xk','MarkerSize',20,'LineWidth',5)
    line([0,m],[0,m],'color','r','LineWidth',2)
    xlabel('Firing frequency (Hz)')
    ylabel('\beta_1')
    gcaformat(gca,false)

%{
figureNB;
for k = 1:10
    Y = psd_active(:,k);
    B(:,k) = regress(Y,X);
    B2(k) = mean((Y-X1)./X1);
    subplot(2,5,k)
    plot(Y,'color','k')
    hold on;
    % plot(X1,'color',[0.6,0.6,0.6])
    plot(X1+N(k)/2*psd_unit(:,k),'color','b','LineWidth',2)
    set(gca,'xscale','log')
    xlim([1,2e3])
    ylim([1e-20,1e-14]);
    set(gca,'yscale','log')
    title(sprintf('%d Hz',round(N(k)/2)))
    xticks([1,10,100,1000])
    xticklabels([1,10,100,1000])
    xlabel('Frequency (Hz)')
    ylabel(['Power (' char(956) 'V^2/Hz)'])
    gcaformat;
end
%}


% figureNB;
%     for i = 1:size(psd_active,2)
%         plot(f,psd_active(:,i),'color',clrs(i,:))
%         hold on
%         set(gca,'xscale','log')
%         set(gca,'yscale','log')
%     end
%     colormap(jet(10))
%     colorbar;
%     C = get(gca,'colorbar')
%     C.Ticks = 1/22:1/10:(1-1/22)
%     C.Label.String = 'E:I ratio'
%     C.TickLabels = {'1','1.5','2.1','3.1','4.5','6.6','9.7','14.1','20.6','30'};
%     xlim([1,1e4])
%     xlabel('Frequency (Hz)')
%     ylabel(['Power (' char(956) 'V^2/Hz)'])
%     gcaformat(gcf)