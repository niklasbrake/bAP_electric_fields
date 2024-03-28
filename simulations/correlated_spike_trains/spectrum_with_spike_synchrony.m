load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\unitaryAP_all.mat')
[sa,X] = network_simulation_beluga.getHeadModel;

sp = synthetic_spikes(0.2,1e-3,1e3);
sp2 = zeros(size(sp));

s0 = 8*sqrt(2);
for i = 1:1e3
    idcs = find(sp(:,i));
    idcs = round(min(max(idcs+16e3*s0*randn(size(idcs)),1),size(sp2,1)));
    sp2(idcs,i) = 1;
end

bs_sample = sample_blue_neurons(1e3);
uAP = network_simulation_beluga.getEEG(savedUnitaryAP,sa,43e3);
for i = 1:length(R)
    sp_ap(:,i) = filter(uAP(:,bs_sample(i)),1,sp(:,i));
    sp2_ap(:,i) = filter(uAP(:,bs_sample(i)),1,sp2(:,i));
end

[psd,f] = pmtm(sum(sp_ap,2),2,[],16e3);
psd = smooth(psd,100);

psd2 = pmtm(sum(sp2_ap,2),2,[],16e3);
psd2 = smooth(psd2,100);


figureNB;
    plot(f,psd,'k','LineWidth',1)
    hold on;
    plot(f,psd2,'r','LineWidth',1)
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlim([1,3e3])
    xlabel('Frequency (Hz)')
    ylabel(['Spectral density (' char(956) 'V^2/Hz)'])
    ylim([1e-16,1e-12])
    gcaformat