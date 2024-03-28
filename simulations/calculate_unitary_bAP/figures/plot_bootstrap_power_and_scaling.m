load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\unitaryAP_all.mat')
load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\mtype_abundance.mat','mtype_abundance');
load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\unitary_AP_PSD.mat','freq','psd_unit','mtype');
[J,ID] = findgroups(mtype);

for i = 1:length(ID)
    abundance(i) = mtype_abundance(ID{i},:).Abundance;
end

mtype_abundance.count = splitapply(@(x) sum(x),ones(size(J)),J)';
scaledAbundance = mtype_abundance.Abundance./mtype_abundance.count;
prop = scaledAbundance(J);
F = cumsum(prop)/sum(prop);
[~,idcs] = unique(F);
m = length(J);
R = rand(m,1e3);
bootstrap_Est = zeros(size(psd_unit,1),1e3);
h = waitbar(0);
for i = 1:1e3
    update_waitbar(h,i,1e3);
    bs_sample = interp1(F(idcs),idcs,R(:,i),'next','extrap');
    bootstrap_Est(:,i) = nanmean(psd_unit(:,bs_sample),2);
end

calculate_dendrite_asymmetry;

C = mtype_abundance(mtype,:).Abundance;



figureNB(14.5,4.5);
subplot(1,3,1);
    scatter(vecnorm(asym_idx(~ei_type,:),2,2),pEst(~ei_type),1+C(~ei_type),zeros(sum(~ei_type),1),'filled','MarkerFaceAlpha',0.3);
    hold on;
    scatter(vecnorm(asym_idx(ei_type,:),2,2),pEst(ei_type),1+C(ei_type),ones(sum(ei_type),1),'filled','MarkerFaceAlpha',0.7);
    colormap([0,0,1;1,0,0]);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    t = log10(get(gca,'xlim'));
    t = linspace(t(1),t(end),1e3);
    FT = fitlm((vecnorm(asym_idx(:,:),2,2)),(pEst(:)'),'intercept',false,'RobustOpts',true);
    FT.Rsquared.Ordinary
    plot(10.^t,FT.predict(10.^t(:)),'-k')
    xlabel('Dendrite asymmetry index')
    ylabel(['Unitary AP power (' char(956) 'V^2)'])
    ylim([3e-16,1e-12])
    gcaformat
subplot(1,3,2)
    plotwitherror(freq,bootstrap_Est,'SE','color','k','LineWidth',1)
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    ylabel(['Unit apEEG power (' char(956) 'V^2/Hz)'])
    xlabel('Frequency (Hz)')
    xlim([1,2e3])
    xticks([1,10,100,1000])
    xticklabels([1,10,100,1000])
    gcaformat;
subplot(1,3,3)
    plot(firingFrequency,B2(:),'.k')
    line([0.1,70],[0.1,70],'color','r','LineWidth',1)
    xlim([0.1,100]);
    ylim([0.1,100]);
    yticks([0.1,1,10,100]);
    yticklabels([0.1,1,10,100])
    xticks([0.1,1,10,100]);
    xticklabels([0.1,1,10,100])
    xlabel('Firing frequency (Hz)')
    ylabel('apEEG scaling factor')
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    gcaformat;