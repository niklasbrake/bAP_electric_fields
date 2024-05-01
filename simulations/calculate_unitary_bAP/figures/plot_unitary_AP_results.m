load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\unitaryAP.mat','mtype','ei_type')
load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\mtype_abundance.mat','mtype_abundance');
load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\unitarySpectrum.mat')
% m=1000;
% idcs = sample_blue_neurons(1035*m);
% idcs = reshape(idcs,[m,1035]);

% Rxx = zeros(size(psd,1),1);
% for i = 1:m
%     Rxx = Rxx+mean(psd(:,idcs(i,:)),2)/m;
% end

calculate_dendrite_asymmetry;
ai = vecnorm(asym_idx,2,2);

red = clrsPT.qualitative_CM.red;
blue = [17,82,185]/255;

[f0,S] = import_Scheer2006;

low_noise = (8e-3).^2;
lam = 8.5;
R = 0.2;
N2 = 16e9;

pEst = sum(psd)*mean(diff(freq));
C = mtype_abundance(mtype,:).Abundance;
figureNB(13.2,4.5);
axes('Position',[0.1, 0.22, 0.23, 0.68])
    scatter(vecnorm(asym_idx(~ei_type,:),2,2),pEst(~ei_type),1+C(~ei_type),zeros(sum(~ei_type),1),'filled','MarkerFaceAlpha',0.3);
    hold on;
    scatter(vecnorm(asym_idx(ei_type,:),2,2),pEst(ei_type),1+C(ei_type),ones(sum(ei_type),1),'filled','MarkerFaceAlpha',0.7);
    colormap([0.5,0.5,0.8;0.8,0.5,0.5]);
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
    xticklabels([100,1000,10000])
    gcaformat(gca,true,8);
axes('Position',[0.425, 0.22, 0.23, 0.68])
    plot(freq,Rxx,'color',blue,'Linewidth',1);
    xlim([10,1e3])
    set(gca,'xscale','log')
    xlabel('Frequency (Hz)')
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    xlim([1,1e3])
    xticks([1,10,100,1000]);
    xticklabels([1,10,100,1000]);
    gcaformat;
axes('Position',[0.75, 0.22, 0.23, 0.68])
    plot(f0,S,'color',0*[0.6,0.6,0.6],'Linewidth',1);
    hold on;

    %N Noise threshold (doi.org/10.1088/0967-3334/27/2/002)
    plot([1,3e3],low_noise*[1,1],'color',red,'LineWidth',1.5,'LineStyle','-')

    plot(freq,0.1*N2*Rxx,'color',blue,'LineWidth',1)
    plot(freq,1*N2*Rxx,'color',blue,'LineWidth',1)
    plot(freq,10*N2*Rxx,'color',blue,'LineWidth',1)
    plot(freq,100*N2*Rxx,'color',blue,'LineWidth',1)

    set(gca,'xscale','log')
    set(gca,'yscale','log');
    xlabel('Frequency (Hz)');
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    ylim([1e-8,1e1])
    yticks([1e-8,1e-4,1e0])
    xlim([1,3e3])
    xticks([1,10,100,1000]);
    xticklabels([1,10,100,1000]);
gcaformat(gcf,true,8)