load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\unitaryAP_all.mat')
load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\mtype_abundance.mat','mtype_abundance');

calculate_dendrite_asymmetry
V = squeeze(max(vecnorm(savedUnitaryAP,2,2)));
% X = [N(:),abs(asym_idx(:,2)./N),V(:)];
% I = (10*(X(:,2)-80)+(X(:,1)-1e3) < X(:,3));
% ai = vecnorm(asym_idx,2,2);
% I = ai<1e5;
% I([21,24,41])=0;
ai = vecnorm(asym_idx,2,2);

C = mtype_abundance(mtype,:).Abundance;
figureNB(9,7.5);
    scatter(vecnorm(asym_idx(~ei_type,:),2,2),pEst(~ei_type),1+2*C(~ei_type),zeros(sum(~ei_type),1),'filled','MarkerFaceAlpha',0.3);
    hold on;
    scatter(vecnorm(asym_idx(ei_type,:),2,2),pEst(ei_type),1+2*C(ei_type),ones(sum(ei_type),1),'filled','MarkerFaceAlpha',0.7);
    colormap([0,0,1;1,0,0]);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    t = log10(get(gca,'xlim'));
    t = linspace(t(1),t(end),1e3);
    FT = fitlm((vecnorm(asym_idx(:,:),2,2)),(pEst(:)'),'intercept',false,'RobustOpts',true)
    plot(10.^t,FT.predict(10.^t(:)),'-k')
    xlabel('Dendrite asymmetry index')
    ylabel(['Unitary AP power (' char(956) 'V^2)'])
    ylim([3e-16,1e-12])
    gcaformat

t = linspace(1,7e5,1e3)';
figureNB;
    scatter(ai,V,20,I,'filled');
    hold on;
    FT = fitlm(ai(I==1),V(I==1));
    FT.Rsquared.Ordinary
    plot(t,FT.predict(t),'k')
    FT = fitlm(ai(I==0),V(I==0));
    FT.Rsquared.Ordinary
    plot(t,FT.predict(t),'k')
    ylim([0,400]);
    xlabel('Dendrite asymmetry index')
    ylabel(['Unitary AP dipole magnitude (nA' char(956) 'm)'])


[l,ID] = findgroups(layer);;
clrs = 5-findgroups(layer);
abun = 2+2*mtype_abundance(mtype,:).Abundance;
figureNB(20,8);
    scatter(ai(ei_type),V(ei_type),abun(ei_type),clrs(ei_type),'filled','Marker','o')
    hold on;
    scatter(ai(~ei_type),V(~ei_type),abun(~ei_type),clrs(~ei_type),'Marker','o','Linewidth',1.5)
    colormap(clrsPT.lines(5))
    set(gca,'CLim',[0,5]);
    C = colorbar
    C.Ticks = 5*(1/10:1/5:(1-1/10));
    C.TickLabels = flip({'L1','L2/3','L4','L5','L6'});
    xlabel('Dendrite asymmetry index')
    ylabel(['Unitary AP dipole magnitude (nA' char(956) 'm)'])





figureNB(20,8);
for i = 1:5
    subplot(2,3,i);
    jj = (l==i);
    scatter(vecnorm(asym_idx,2,2),V(:),20,[0.6,0.6,0.6],'filled','Marker','o')
    hold on;
    idcs = find(ei_type.*jj);
    scatter(vecnorm(asym_idx(idcs,:),2,2),V(idcs),20,'r','filled','Marker','o')
    hold on;
    idcs = find(~ei_type.*jj);
    scatter(vecnorm(asym_idx(idcs,:),2,2),V(idcs),20,'b','filled','Marker','o')
    xlabel('Dendrite asymmetry index')
    ylabel(['Unitary AP dipole magnitude (nA' char(956) 'm)'])
    title(ID{i})
    gcaformat;
end