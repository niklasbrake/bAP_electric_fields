warning('off','signal:findpeaks:largeMinPeakHeight');
% [sa,X] = network_simulation_beluga.getHeadModel;

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
savedUnitaryAP = nan(2001,M);
N = zeros(length(EI_vec),1);
B2 = zeros(length(EI_vec),M);
apCount = zeros(length(EI_vec),M);

fs = 16e3; % Hz

for ii = 1:M
    i = iVec(ii);
    waitbar(ii/M)
    folder0 = fullfile(folder,F(i).name,'matlab_recordings');

    meanAP = [];
    for k = 1:length(EI_vec)
        load(fullfile(folder0,sprintf('synaptic_input_EI%s_passive.mat',EI_vec{k})))
        eeg = network_simulation_beluga.getEEG(dipoles,sa,1e3);
        [f0,~,psd] = eegfft(time*1e-3,detrend(eeg),0.5,0.4);
        psd_passive(:,k) = mean(psd,2);

        load(fullfile(folder0,sprintf('synaptic_input_EI%s.mat',EI_vec{k})))
        eeg = network_simulation_beluga.getEEG(dipoles,sa,1e3);
        [f0,~,psd] = eegfft(time*1e-3,detrend(eeg),0.5,0.4);
        psd_active(:,k) = mean(psd,2);

        [y,x] = findpeaks(voltage,'MinPeakHeight',0);
        N(k) = length(x);
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

    savedUnitaryAP(:,ii) = nanmean(meanAP);
    m = size(meanAP,1);

    paddedAP = [zeros(m,7e3),meanAP,zeros(m,7e3-1)]';
    paddedAP(isnan(paddedAP(:))) = 0;
    paddedAP = mean(paddedAP,2);

    n = length(paddedAP);
    xdft = fft(paddedAP);
    xdft = xdft(1:n/2+1);
    psdx = (1/(fs*n)) * abs(xdft).^2;
    psdx(2:end-1) = 2*psdx(2:end-1);
    freq = 0:fs/n:fs/2;

    psd_unit = interp1(freq,psdx,f0,'linear','extrap');
    psd_unit = psd_unit*length(paddedAP)/fs/1; % per second

    for k = 1:10
        X0 = psd_passive(:,k);
        X1 = psd_unit;
        Y = psd_active(:,k);
        idcs = find(and(and(Y>1e-19,X1>1e-19),f0<1e3));
        B2(k,ii) = fminbnd(@(B1) neglnlike(B1,X0(idcs),X1(idcs),Y(idcs)),0,N(k)*2);
        apCount(k,ii) = N(k);
    end
end
delete(h)

firingFrequency = apCount(:)/range(time)*1e3;
figureNB;
plot(firingFrequency,B2(:),'.k')
line([0.1,70],[0.1,70],'color','r','LineWidth',1)
xlabel('Firing frequency (Hz)')
ylabel('\beta_1')
set(gca,'xscale','log')
set(gca,'yscale','log')
return;

Y = savedUnitaryAP./max(abs(savedUnitaryAP));
% Y = Y(:,idcs);
Y = Y-median(Y);
Y = Y/max(Y(:));
[coeff,score,~,~,explained] = pca(Y');


% S.mu = [-3,-2;10,2];
% S.Sigma = cat(3,[0.2,0;0,0.2],[1,-1;-1,5]);
% GMModel = fitgmdist(score(:,1:2),2,'Start',S);
GMModel = fitgmdist(score(:,1:2),2);
gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(GMModel,[x0 y0]),x,y);
J = GMModel.cluster(score(:,1:2));

figureNB(9,9);
subplot(2,1,1)
    % gscatter(score(:,1),score(:,2),layer(idcs)');
    plot(score(:,1),score(:,2),'.k','MarkerSize',10);
    hold on;
    g = gca;
    fcontour(gmPDF,[[-4,4],[-4,4]],'LineWidth',1,'MeshDensity',100)
    % fcontour(gmPDF,[[-4,4],[-4,4]],'LineWidth',1,'MeshDensity',100,'LevelList',[0.05,0.1,0.3,0.4,1])
    % plot(score(:,1),score(:,2),'.k','MarkerSize',10);
    gscatter(score(:,1),score(:,2),J)
    legend off
    xlabel('PC 1')
    ylabel('PC 2')
    set(gca,'DataAspectRatio',100-explained(1:3))
    % xlim([-2,2])
    % ylim([-1,1])
    % set(gca,'DataAspectRatio',[1,1,1])

subplot(2,2,3)
    plot((-1000:1000)/16,mean(Y(:,J==1),2),'LineWidth',1.5);
    hold on;
    plot((-1000:1000)/16,mean(Y(:,J==2),2),'LineWidth',1.5);
    % plot((-1000:1000)/16,mean(Y(:,J==3),2),'LineWidth',1.5);
    xlim([-5,5]);
    ylabel(['Unitary AP (norm.)'])
    xlabel('Time (ms)')


y1 = mean(Y(:,J==1),2);
y1 = y1-median(y1);
[psd1,f] = pmtm(y1,2,[],fs);

y2 = mean(Y(:,J==2),2);
y2 = y2-median(y2);
[psd_active,f] = pmtm(y2,2,[],fs);

% y3 = mean(Y(:,J==3),2);
% y3 = y3-median(y3);
% [psd3,f] = pmtm(y3,2,[],fs);



subplot(2,2,4)
    plot(f,psd1,'LineWidth',2);
    hold on;
    plot(f,psd_active,'LineWidth',2)
    % plot(f,psd3,'LineWidth',2)
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlabel('Frequency (Hz)')
    ylabel('Power (a.u.)')
    yticks([])

gcaformat(gcf);

Y = savedUnitaryAP-median(savedUnitaryAP);
Y = [zeros(1e4,M);Y;zeros(1e4,M)];

[psd,f] = pmtm(Y,2,[],fs);
figureNB;
subplot(1,2,1);
    plotwitherror((-1000:1000)/16,Y(1e4:1e4+2000,:),'Q','color','k','LineWidth',1);
    xlabel('Time (ms)')
    ylabel(['Unitary AP (' char(956) 'V)'])
    ylim([-2e-6,2e-6])
    xlim([-2,10])
    xticks([0:5:10])
    gcaformat
subplot(1,2,2);
    plotwitherror(f(2:end),psd(2:end,:),'Q','color','k','LineWidth',1);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlabel('Frequency (Hz)')
    ylabel(['Power (' char(956) 'V^2/Hz)'])
    gcaformat