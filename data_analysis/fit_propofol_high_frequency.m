
folder = fullfile(dataFolder,'simulations','subcritical_network','PSD');
F = dir(folder); F = F(3:end);
for i = 1:length(F)
    load(fullfile(folder,F(i).name));
    psd(:,:,i) = P;
end

y = mean(psd(:,:,end)*1.66e15,2);
Psyn = interp1(f,y,freq,'spline','extrap');

FT = polyfit(log(f),log(y),5);
Psyn(freq>f(end)) = exp(polyval(FT,log(freq(freq>f(end)))));

load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\MC_cortex_sampling\Pxy.mat')
f2 = 10.^linspace(0,log10(3e3),100);
Pxy(Pxy<0) = 0;
Pap = interp1(f2,Pxy,freq,'linear','extrap');
Pap(freq<50) = 0;

low_noise = (30*1e-3).^2;

[full_model,AP_model] = fittingmodel('eq6');
Psyn = 10.^AP_model(freq,[40e-3,6e-3,-Inf,3.5]);

figureNB;
    plot(freq,Psyn+low_noise,'-r','LineWidth',1);
    hold on;
    plot(freq,Psyn*0+low_noise,'--k');
    plot(freq,Psyn,'--r');
    plot(freq,Psyn+5e10*Pap+low_noise);
    plot(freq,Psyn+1e12*Pap+low_noise);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlim([1,1e3])
    ylim([1e-5,1e2])
    xlabel('Frequency (Hz)');
    ylabel(['PSD (' char(956) 'V^2/Hz)'])




load('E:\Research_Projects\004_Propofol\data\experiments\scalp_EEG\raw\time_series_all_channels.mat')
Cz = squeeze(TimeDomainAligned(:,2,:));
Cz_CAF = squeeze(TimeDomainAligned(:,2,:)-nanmean(TimeDomainAligned(:,1:7,:),2));
for i = 1:14
[datafreq,time,psd_Cz(:,:,i)] = eegfft(Time,Cz(:,1),2,0);
[datafreq,time,psd_Cz_CAF(:,:,i)] = eegfft(Time,Cz_CAF(:,1),2,0);
end

pre = squeeze(nanmedian(psd_Cz(:,time<-200,:),2));
post = squeeze(nanmedian(psd_Cz(:,time>0,:),2));
post_CAF = squeeze(nanmedian(psd_Cz_CAF(:,time>0,:),2));
pre_CAF = squeeze(nanmedian(psd_Cz_CAF(:,time>-200,:),2));

noise_floor = 2.8e-4;
figureNB;
subplot(1,2,1);
    plotwitherror(datafreq,pre_CAF,'M','LineWidth',1,'color','g');
    plot(freq,10.^AP_model(freq,[15e-3,3e-3,-Inf,2.4])+noise_floor,'k');
    plot(freq,freq*0+noise_floor,'--k');
    xlim([1,250])
    set(gca,'xscale','log')
    set(gca,'yscale','log')

    plot(freq,10.^AP_model(freq,[15e-3,4e-3,-Inf,2.4])+noise_floor+10.^AP_model(freq,[2.5e-3,0.5e-3,-Inf,2.15]),'b','LineWidth',1);
subplot(1,2,2);
    plotwitherror(datafreq,post_CAF,'M','LineWidth',1,'color','g');
    plot(freq,10.^AP_model(freq,[35e-3,7e-3,-Inf,2.5])+noise_floor,'k');
    hold on;
    plot(freq,freq*0+noise_floor,'--k');
    xlim([1,250])
    set(gca,'xscale','log')
    set(gca,'yscale','log')

