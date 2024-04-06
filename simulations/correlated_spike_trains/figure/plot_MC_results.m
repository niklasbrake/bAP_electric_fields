extracted_Scheer2006;
close all;

load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\MC_cortex_sampling\Pxx_brute.mat')
load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\MC_cortex_sampling\Pxy.mat');
count_D = count_D';
fs =16e3;
L = 2e3;
f = fs/L:fs/L:fs/2;
mu = Pxy_D./count_D;
CI = 1.96*sqrt(SSExy_D./(count_D-1))./sqrt(count_D);

%{
figureNB;
    plot(dValues,count_D,'k','Linewidth',1)
    xlabel('Pairwise distance (mm)');
    ylabel('n samples');
    title(sprintf('Total samples: %d',sum(count_D)));
    gcaformat;


x = [f,flip(f)];
figureNB(30,12);
    y = [mu-CI;flip(mu+CI)];
    for i = 1:50:length(dValues)
        plot3(f,f*0+dValues(i),mu(:,i),'k','Linewidth',1)
        hold on;
        fill3(x,x*0+dValues(i),y(:,i),'k','FaceAlpha',0.2,'EdgeColor','none');
        set(gca,'xscale','log')
        drawnow;
    end
    ylim([0,dValues(end)])
zlim([-1e-17,1e-17])
zlabel(['PSD (' char(956) 'V^2/Hz)'])
ylabel('Distance (mm)')
xlabel('Frequency (Hz)')
view([-70,50]);
ylim([0,40]);


x = [f,flip(f)];
figureNB(30,12);
    y = [mu-CI;flip(mu+CI)];
    for i = 1:50:length(dValues)
        plot3(f,f*0+dValues(i),mu(:,i).*exp(-dValues(i).^2/6),'k','Linewidth',1)
        hold on;
        fill3(x,x*0+dValues(i),y(:,i).*exp(-dValues(i).^2/6),'k','FaceAlpha',0.2,'EdgeColor','none');
        set(gca,'xscale','log')
        drawnow;
    end
    ylim([0,dValues(end)])
zlim([-1e-17,1e-17])
zlabel(['PSD (' char(956) 'V^2/Hz)'])
ylabel('Distance (mm)')
xlabel('Frequency (Hz)')
view([-70,50]);
ylim([0,40]);
%}


x = [f,flip(f)];
y = [mu-CI;flip(mu+CI)];
y(isinf(y(:))) = nan;


load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\cortical_area\pairwise_distance.mat')
dV = mean(diff(dValues));
DV1 = [dValues,dValues(end)+mean(diff(dValues))];
dN = interp1(0.5*(rValues(1:end-1)+rValues(2:end)),diff(total_area)./diff(rValues'),dValues,'linear','extrap'); 
dN = mean(dN,2)'; % d(area) / d(radius) (r)
dN = 2e5*dN; % d(neurons) / d(radius) (r)

N2 = 16e9;

Rxy = nansum(mu.*exp(-dValues.^2/6).*dN*dV,2)/N2;
Rxy_CI = nansum(y.*exp(-dValues.^2/6).*dN*dV,2)/N2;

figureNB;
subplot(1,2,1);
    plot(freq,p_avg,'k','Linewidth',1);
    xlim([10,1e3])
    set(gca,'xscale','log')
    xlabel('Frequency (Hz)')
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    gcaformat;

subplot(1,2,2);
    plot(f,Rxy,'k','Linewidth',1);
    hold on;
    fill(x,Rxy_CI,'k','FaceAlpha',0.2,'EdgeColor','none');;
    set(gca,'xscale','log')
    xlim([10,1e3])
    set(gca,'xscale','log')
    xlabel('Frequency (Hz)')
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    gcaformat;

low_noise = (16/1e3)^2;;
lam = 8.5;

figureNB(8,7);
    plot(f0,S,'r','Linewidth',1);
    hold on;

    %N Noise threshold (doi.org/10.1088/0967-3334/27/2/002)
    low_noise = (8*1e-3).^2;
    plot([1,3e3],low_noise*[1,1],'k','LineWidth',1,'LineStyle','--')

    % plot(f,lam*N2*p_avg,'color',[0.6,0.6,0.6],'LineWidth',1)
    % plot(f,lam*N2*p_avg + 0.2*lam*N2*(N2-1)*Rxy,'color',[0,0,0],'LineWidth',1)

    dS = [0,10.^linspace(-3,-1.5,6)];
    for i = 1:length(dS)
        sig = dS(i);
        B = exp(-(2*pi*f(:)*sig).^2/2);
        if(i==1)
            plot(f,lam*N2*p_avg + 0.2*lam*N2*(N2-1)*B.*Rxy,'color',[0,0,0],'LineWidth',1)
        else
            plot(f,lam*N2*p_avg + 0.2*lam*N2*(N2-1)*B.*Rxy,'color',[0.6,0.6,0.6],'LineWidth',1)
        end
    end

    sig = sqrt(2)*8e-3;
    B = exp(-(2*pi*f(:)*sig).^2/2);
    plot(f,lam*N2*p_avg + 0.2*lam*N2*(N2-1)*B.*Rxy,'color','b','LineWidth',1)
    set(gca,'xscale','log')
    set(gca,'yscale','log');
    xlabel('Frequency (Hz)');
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    xlim([30,3e3])
    xticks([30,300,3e3])
    xticklabels([30,300,3e3])
    gcaformat;





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


figureNB(8,7);
    plotwitherror(freq,post_CAF,'M','LineWidth',1,'color','r');
    plotwitherror(freq,pre_CAF,'M','LineWidth',1,'color','g');
    low_noise = (16/1e3)^2;
    plot(freq,freq*0+low_noise,'--k');
    plot(f,lam*N2*p_avg + 0.2*lam*N2*(N2-1)*B.*Rxy,'color','k','LineWidth',1)

    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlim([30,300])
    xticks([30,100,300])
    xticklabels([30,100,300]);

    set(gca,'xscale','log')
    set(gca,'yscale','log');
    xlabel('Frequency (Hz)');
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    gcaformat;