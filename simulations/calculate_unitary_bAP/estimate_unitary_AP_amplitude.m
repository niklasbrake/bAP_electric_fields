% load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP\unitary_AP_analysis.mat')
load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\unitaryAP_all.mat')

load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\mtype_abundance.mat','mtype_abundance');

% folder = 'E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\sampled_mtypes';
[J,ID] = findgroups(mtype);

for i = 1:length(ID)
    abundance(i) = mtype_abundance(ID{i},:).Abundance;
end

m_pEst = splitapply(@(x) nanmean(x,2),pEst,J(:)');
s0_pEst = splitapply(@(x) -1.96*stderror(x')+nanmean(x,2),pEst,J(:)');
s1_pEst = splitapply(@(x) 1.96*stderror(x')+nanmean(x,2),pEst,J(:)');

% Bootstrap average
mtype_abundance.count = splitapply(@(x) length(x),pEst,J(:)')';
scaledAbundance = mtype_abundance.Abundance./mtype_abundance.count;
prop = scaledAbundance(J);
F = cumsum(prop)/sum(prop);
[~,idcs] = unique(F);
m = length(J);
R = rand(m,10000);
for i = 1:1e3
    bs_sample = interp1(F(idcs),idcs,R(:,i),'next','extrap');
    bootstrap_Est(i) = nanmean(pEst(bs_sample));
end


% Bootstrap average (e-cells)
idcs0 = find(ei_type);
prop = scaledAbundance(J);
prop = prop(idcs0);
F = cumsum(prop)/sum(prop);
[~,idcs] = unique(F);
m = length(J);
R = rand(m,10000);
for i = 1:1e3
    bs_sample = interp1(F(idcs),idcs0(idcs),R(:,i),'next','extrap');
    ebootstrap_Est(i) = nanmean(pEst(bs_sample));
end

% Bootstrap average (i-cells)
idcs0 = find(~ei_type);
prop = scaledAbundance(J);
prop = prop(idcs0);
F = cumsum(prop)/sum(prop);
[~,idcs] = unique(F);
m = length(J);
R = rand(m,10000);
for i = 1:1e3
    bs_sample = interp1(F(idcs),idcs0(idcs),R(:,i),'next','extrap');
    ibootstrap_Est(i) = nanmean(pEst(bs_sample));
end


fprintf('Excitatory cell uAP energy: %.1f +/- %.1f zV^2\n',10.^(log10(mean(ebootstrap_Est))+15),10.^(log10(std(ebootstrap_Est))+15))

fprintf('Inhibitory cell uAP energy: %.1f +/- %.1f zV^2\n',10.^(log10(mean(ibootstrap_Est))+15),10.^(log10(std(ibootstrap_Est))+15))

figureNB(30,15);
    bar(m_pEst,'EdgeColor','k','FaceColor','k');
    line([1:55;1:55],[s0_pEst;s1_pEst],'color','k')
    xticks(1:55)
    xticklabels(cellfun(@(x)strrep(x,'_','\_'),ID,'UniformOutput',false))
    set(gca,'yscale','log')
    ylabel(['Unitary AP power (' char(956) 'V^2)'])
    gcaformat
    set(gca,'XTickLabelRotation',45)
    ylim([1e-15,2e-13])
    yticks([1e-15,1e-14,1e-13])
    yticks([1e-15,1e-14,1e-13])


figureNB(8,13)
axes('Position',[0.16, 0.11, 0.60, 0.81])
    plot(pEst,56-J,'.k')
    set(gca,'xscale','log')
    xlabel(['Unitary AP power (' char(956) 'V^2)'])
    yticks(1:55)
    set(gca,'FontSize',5)
    yticklabels(cellfun(@(x)strrep(x,'_','\_'),flip(ID),'UniformOutput',false))
    gcaformat;
    yax = get(gca,'yaxis');
    ylim([0,56])
    yax.FontSize = 5;

axes('Position',[0.76, 0.11, 0.09, 0.81])
    scatter(ones(55,1),55-(0:54),max(2,abundance),abundance,'filled');
    ylim([0,56])
    axis off
    idcs = find(abundance>2);
    for i = idcs
        text(2,56-i,sprintf('%.1f%c',abundance(i),char(37)),'FontSize',6,'VerticalAlignment','middle');
    end
    colormap(flip(copper))


[I,ID2] = findgroups(layer);
figureNB(9,4);
for i = 1:5
    idcs = find(and(I==i,ei_type==0));
    boxplotNB(3*(i-1),log10(pEst(idcs)),'b',5);
    hold on;
    idcs = find(and(I==i,ei_type==1));
    if(i>1)
        boxplotNB(3*(i-1)+1,log10(pEst(idcs)),'r',15);
    end
end
xticks([0.5:3:12.5])
xticklabels({'Layer 1','Layer 2/3','Layer 4','Layer 5','Layer 6'})
xlim([-1,14])
ylabel(['Unitary AP power (' char(956) 'V^2)'])
yticks([-15:-13])
yticklabels({'10^{-15}','10^{-14}','10^{-13}'});
gcaformat



figureNB(5,6);
    boxplotNB(1,log10(pEst(ei_type==0)),'b',15);
    hold on;
    boxplotNB(2,log10(pEst(ei_type==1)),'r',15);
    ylabel(['Unitary AP power (' char(956) 'V^2)'])
    yticks([-15:-13])
    yticklabels({'10^{-15}','10^{-14}','10^{-13}'});
    xticks([1,2])
    xticklabels({'Inhibitory','Excitatory'})
    gcaformat



relCont = (m_pEst(:).*abundance(:))./nansum(m_pEst(:).*abundance(:));
[~,I] = sort(relCont,'descend');
ei = splitapply(@mean,ei_type,J)';
y = relCont(I);
figureNB;
    bar(y,'k');
    xticks(1:55)
    xticklabels(cellfun(@(x)strrep(x,'_','\_'),ID(I),'UniformOutput',false))
    yyaxis right;
    hold on;
    idEI = find(~ei(I));
    plot([idEI;55],[cumsum(y(idEI));sum(y(idEI))],'-b','Linewidth',1)
    idEI = find(ei(I));
    plot([idEI;55],[cumsum(y(idEI));sum(y(idEI))],'-r','Linewidth',1)
    gcaformat
    set(gca,'XTickLabelRotation',45)



for i = 1:5
    idcs = find(and(I==i,ei_type==0));
    x(2*(i-1)+1) = 3*(i-1);
    y{2*(i-1)+1} = log10(pEst(idcs))

    idcs = find(and(I==i,ei_type==1));
    x(2*(i-1)+2) = 3*(i-1)+1;
    y{2*(i-1)+2} = log10(pEst(idcs))

    clrs(2*(i-1)+1,:) = [0,0,1];
    clrs(2*(i-1)+2,:) = [1,0,0];
end
x = x([1,3:end]);
y = y([1,3:end]);
clrs = clrs([1,3:end],:);
