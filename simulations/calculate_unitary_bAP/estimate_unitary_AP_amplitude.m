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



i0 = find(strcmp(ID,'L6_UTPC'))
i1 = find(strcmp(ID,'L6_NBC'))
    temp1 = J(i0);
    temp2 = pEst(i0);
    temp3 = ID(i0);
    temp4 = abundance(i0);
    J(i0) = J(i1);
    pEst(i0) = pEst(i1);
    ID(i0) = ID(i1);
    abundance(i0) = abundance(i1);
    J(i1) = temp1;
    pEst(i1) = temp2;
    ID(i1) = temp3;
    abundance(i1) = temp4;


i0 = find(strcmp(ID,'L6_TPC_L4'));
i1 = find(strcmp(ID,'L6_NBC'))
    temp1 = J(i0);
    temp2 = pEst(i0);
    temp3 = ID(i0);
    temp4 = abundance(i0);
    J(i0) = J(i1);
    pEst(i0) = pEst(i1);
    ID(i0) = ID(i1);
    abundance(i0) = abundance(i1);
    J(i1) = temp1;
    pEst(i1) = temp2;
    ID(i1) = temp3;
    abundance(i1) = temp4;

i0 = find(strcmp(ID,'L5_TTPC1'))
i1 = find(strcmp(ID,'L5_NBC'))
    temp1 = J(i0);
    temp2 = pEst(i0);
    temp3 = ID(i0);
    temp4 = abundance(i0);
    J(i0) = J(i1);
    pEst(i0) = pEst(i1);
    ID(i0) = ID(i1);
    abundance(i0) = abundance(i1);
    J(i1) = temp1;
    pEst(i1) = temp2;
    ID(i1) = temp3;
    abundance(i1) = temp4;

figureNB(12,6);
axes('Position',[0.1, 0.22, 0.8, 0.62])
    plot(J,pEst,'.k')
    set(gca,'yscale','log')
    ylabel(['Unitary AP power (' char(956) 'V^2)'])
    % xticks(1:55)
    set(gca,'FontSize',5)
    % xticklabels(cellfun(@(x)strrep(x,'_','\_'),ID,'UniformOutput',false))
    gcaformat;
    xax = get(gca,'xaxis');
    xlim([0,56])
    set(gca,'XTickLabelRotation',90)
    xax.FontSize = 5;

axes('Position',[0.1, 0.9, 0.8, 0.07])
    scatter(1:55,ones(55,1),max(2,abundance),abundance,'filled');
    xlim([0,56])
    axis off
    idcs = find(abundance>6.5);
    for i = idcs
        text(i,2,sprintf('%.1f%c',abundance(i),char(37)),'FontSize',6,'HorizontalAlignment','center');
        % text(i,2,sprintf('%.1f',abundance(i)),'FontSize',6,'HorizontalAlignment','center');
    end
    colormap(flip(copper))




figureNB(10.4,5);
axes('Position',[0.14, 0.18, 0.8, 0.61])
    i = find(~ei_type);
    S = scatter(J(i),pEst(i),8,[0.5,0.5,0.8],'filled','MarkerEdgeColor','k');
    hold on;
    i = find(ei_type);
    S = scatter(J(i),pEst(i),10,[0.8,0.5,0.5],'filled','MarkerEdgeColor','k');
    set(gca,'yscale','log')
    ylabel(['Unitary AP power (' char(956) 'V^2)'])
    gcaformat;
    xax = get(gca,'xaxis');
    xax.MinorTickValues = 1:55;
    xax.MinorTick = 'on';
    xlim([0,56])
    xlabel('Neuron morphology type index')
    set(gca,'FontSize',8)

axes('Position',[0.14, 0.92, 0.8, 0.05])
    [~,I] = sort(abundance,'descend');
    for i = I
        scatter(i,1,2*max(2,abundance(i)),abundance(i),'filled');
        hold on;
    end
    xlim([0,56])
    axis off
    colormap(flip(copper))

axes('Position',[0.145 0.83 0.28 0.06])
    exA = [2,5,10,15];
    scatter([0,0.36,0.65,1],zeros(1,4),2*exA,exA,'filled');
    text(0.05,0,'<2%','FontSize',8,'VerticalAlignment','middle')
    text(0.41,0,'5%','FontSize',8,'VerticalAlignment','middle')
    text(0.7,0,'10%','FontSize',8,'VerticalAlignment','middle')
    text(1.06,0,'15%','FontSize',8,'VerticalAlignment','middle')
    ylim([0,1]);
    xlim([-0.1,1]);
    T = text(-0.1,1,'Abundance','FontSize',8,'HorizontalAlignment','left');
    axis off;


axes('Position',[0.88 0.83 0.02 0.06])

    S = scatter(0,1,12,[0.8,0.5,0.5],'filled','MarkerEdgeColor','k');
    text(1,1,'Ex.','FontSize',8,'VerticalAlignment','middle')

    hold on;
    S = scatter(0,0,12,[0.5,0.5,0.8],'filled','MarkerEdgeColor','k');
    text(1,0,'In.','FontSize',8,'VerticalAlignment','middle')
    axis off;
    ylim([0,1]);