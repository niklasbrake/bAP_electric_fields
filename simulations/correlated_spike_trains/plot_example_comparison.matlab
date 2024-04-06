
load('E:\Research_Projects\004_Propofol\data\experiments\scalp_EEG\raw\time_series_all_channels.mat')
Cz = squeeze(TimeDomainAligned(:,2,:));
Cz_CAF = squeeze(TimeDomainAligned(:,2,:)-nanmean(TimeDomainAligned(:,1:7,:),2));
for i = 1:14
[freq,time,psd_Cz(:,:,i)] = eegfft(Time,Cz(:,1),2,0);
[freq,time,psd_Cz_CAF(:,:,i)] = eegfft(Time,Cz_CAF(:,1),2,0);
end

pre = squeeze(nanmedian(psd_Cz(:,time<-200,:),2));
post = squeeze(nanmedian(psd_Cz(:,time>0,:),2));
post_CAF = squeeze(nanmedian(psd_Cz_CAF(:,time>0,:),2));
pre_CAF = squeeze(nanmedian(psd_Cz_CAF(:,time>-200,:),2));

figureNB;
    plotwitherror(freq,post_CAF,'M','LineWidth',1,'color','r');
    plotwitherror(freq,pre_CAF,'M','LineWidth',1,'color','g');
    low_noise = (8/1e3)^2;
    plot(freq,freq*0+low_noise,'--k');

    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlim([1,512])