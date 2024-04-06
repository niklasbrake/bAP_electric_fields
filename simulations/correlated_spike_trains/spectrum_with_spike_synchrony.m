load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\unitaryAP_all.mat')
[sa,X] = network_simulation_beluga.getHeadModel;

N = 1e3;
sp = synthetic_spikes(0.2,1e-3,N);

s0 = 1e-3;%*sqrt(2);
sp2 = zeros(size(sp));
sp3 = zeros(size(sp));
for i = 1:N
    idcs = find(sp(:,i));
    idcs = round(min(max(idcs+16e3*s0*randn(size(idcs)),1),size(sp2,1)));
    sp2(idcs,i) = 1;

    idcs = randi(size(sp,1),length(idcs),1);
    sp3(idcs,i) = 1;
end

bs_sample = sample_blue_neurons(N);
uAP = network_simulation_beluga.getEEG(savedUnitaryAP,sa,43e3);
for i = 1:N
    sp_ap(:,i) = filter(uAP(:,bs_sample(i)),1,sp(:,i));
    sp2_ap(:,i) = filter(uAP(:,bs_sample(i)),1,sp2(:,i));
    sp3_ap(:,i) = filter(uAP(:,bs_sample(i)),1,sp3(:,i));
end

[psd1,f] = pmtm(sum(sp_ap,2),2,[],16e3);
psd1 = smooth(psd1,100);

psd2 = pmtm(sum(sp2_ap,2),2,[],16e3);
psd2 = smooth(psd2,100);

psd3 = pmtm(sum(sp3_ap,2),2,[],16e3);
psd3 = smooth(psd3,100);

[psd0,f0] = pmtm(sum(uAP(:,bs_sample),2),2,[],16e3);
psd0 = smooth(psd0,100);

psd00 = sum(pmtm(uAP(:,bs_sample),2,[],16e3),2);
psd00 = smooth(psd00,100);


figureNB;
    plot(f,psd,'k','LineWidth',1)
    hold on;
    plot(f,psd2,'r','LineWidth',1)
    plot(f,psd3,'b','LineWidth',1)
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlim([1,3e3])
    xlabel('Frequency (Hz)')
    ylabel(['Spectral density (' char(956) 'V^2/Hz)'])
    ylim([1e-16,1e-11])
    gcaformat

    yyaxis right;
    plot(f0,psd0,'g','LineWidth',1)
    hold on;
    plot(f0,psd00,'m','LineWidth',1)
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlim([1,3e3])
    xlabel('Frequency (Hz)')
    ylabel(['Spectral density (' char(956) 'V^2/Hz)'])
    ylim([1e-15,1e-10])
    % ylim([1e-16,1e-12])