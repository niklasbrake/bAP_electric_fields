load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\unitaryAPNew.mat')
load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\mtype_abundance.mat','mtype_abundance');

calculate_dendrite_asymmetry
V = squeeze(max(vecnorm(savedUnitaryAP,2,2)));
% X = [N(:),abs(asym_idx(:,2)./N),V(:)];
% I = (10*(X(:,2)-80)+(X(:,1)-1e3) < X(:,3));
% ai = vecnorm(asym_idx,2,2);
% I = ai<1e5;
% I([21,24,41])=0;
ai = vecnorm(asym_idx,2,2);

red = clrsPT.qualitative_CM.red;
blue = [17,82,185]/255;

load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\MC_cortex_sampling\Pxx_brute.mat')
Rxx = zeros(4e3,1);
Rxx(1:2:end) = p_avg;
Rxx(2:2:end) = p_avg;
load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\MC_cortex_sampling\Pxy.mat');
count_D = count_D';
fs =16e3;
L = 8e3;
f = fs/L:fs/L:fs/2;
mu = Pxy_D./count_D;
CI = 1.96*sqrt(SSExy_D./(count_D-1))./sqrt(count_D);

x = [f,flip(f)];
y = [mu-CI;flip(mu+CI)];
y(isinf(y(:))) = nan;


load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\cortical_area\pairwise_distance.mat')
dV = mean(diff(dValues));
DV1 = [dValues,dValues(end)+mean(diff(dValues))];
dN = interp1(0.5*(rValues(1:end-1)+rValues(2:end)),diff(total_area)./diff(rValues'),dValues,'linear','extrap'); 
dN = mean(dN,2)'; % d(area) / d(radius) (r)
dN = dN/mean(total_area(end,:))*16e9; % d(neurons) / d(radius) (r)

N2 = 16e9;

Rxy = nansum(mu.*exp(-dValues.^2/6).*dN*dV,2)/N2;
Rxy_CI = nansum(y.*exp(-dValues.^2/6).*dN*dV,2)/N2;


[f0,S] = import_Scheer2006;

low_noise = (8e-3).^2;
lam = 8.5;
R = 0.2;

pEst = sum(psd)*mean(diff(freq));
C = mtype_abundance(mtype,:).Abundance;
figureNB(13.2,4.5);
axes('Position',[0.1, 0.22, 0.23, 0.68])
    scatter(vecnorm(asym_idx(~ei_type,:),2,2),pEst(~ei_type),1+C(~ei_type),zeros(sum(~ei_type),1),'filled','MarkerFaceAlpha',0.3);
    hold on;
    scatter(vecnorm(asym_idx(ei_type,:),2,2),pEst(ei_type),1+C(ei_type),ones(sum(ei_type),1),'filled','MarkerFaceAlpha',0.7);
    colormap([blue;red]);
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
    plot(f(2:2:end),Rxx(2:2:end),'k','Linewidth',1);
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


    clrs = clrsPT.sequential(7);
    clrs = 0.6+0*clrs(4:7,:);
    plot(f,0.1*N2*Rxx,'color',clrs(1,:),'LineWidth',1)
    plot(f,1*N2*Rxx,'color',clrs(2,:),'LineWidth',1)
    plot(f,10*N2*Rxx,'color',clrs(3,:),'LineWidth',1)
    plot(f,100*N2*Rxx,'color',clrs(4,:),'LineWidth',1)

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