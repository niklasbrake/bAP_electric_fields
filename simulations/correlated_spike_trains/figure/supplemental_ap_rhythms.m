[f0,S] = import_paralytic_data;
S = fillgaps(nanmean(S(:,1:8),2),5);

R = 0.2;
lam = 1;
dF = 5*2.^[0:5];
sig = 60e-3;
t = linspace(-20,20,1e3);
low_noise = 1e-3;

clrs = clrsPT.sequential(length(dF)+3);
clrs = clrs(4:end,:);
blue = [17,82,185]/255;

[full_model,AP_model] = fittingmodel('eq6');
Psyn = 10.^AP_model(f,[20e-3,4e-3,-Inf,3.6])+10.^AP_model(f,[3e-3,1e-3,-Inf,3.3])+1e-3;

lam = 1;
A = 10.^linspace(-3,2,500);
P0 = [];

for i = 1:length(A)
    P0(i,:) = 10*log10((lam*N2*Rxx+A(i)*N2*(N2-1)*Rxy)./Psyn(:));
end
[full_model,AP_model] = fittingmodel('eq6');
Psyn = 10.^AP_model(0:1e3,[20e-3,4e-3,-Inf,3.6])+10.^AP_model(0:1e3,[3e-3,1e-3,-Inf,3])+1e-3;
[XX,YY] = meshgrid(f,A);


figureNB(13.2,3.75);
axes('Position',[0.05, 0.26, 0.15, 0.57])
    t = linspace(-60,60,1e3);
    plot(t,R*lam/sqrt(2*pi*sig.^2)*exp(-1e-6*t.^2./2/sig.^2).*cos(2*pi*t*1e-3*40),'color',blue,'LineWidth',2)
    ylim([-20,20]);
    axis off;
    gcaformat;
axes('Position',[0.315, 0.26, 0.17, 0.57])
    plot(f0,S,'color',[0.6,0.6,0.6],'Linewidth',1);
    hold on;
    plot([1,3e3],low_noise*[1,1],'r','LineWidth',0.5,'LineStyle','-')
    plot(0:1e3,Psyn,'k','LineWidth',1)
    F0 = 80;
    B = exp(-2*(pi*(f(:)-F0)*sig).^2);
    plot(f,lam*N2*Rxx+R*lam*N2*(N2-1)*Rxy,'color',blue*0.4+0.6,'LineWidth',1,'LineStyle',':')
    plot(f,lam*N2*Rxx + R*lam*N2*(N2-1)*B.*Rxy,'color',blue,'LineWidth',1)
    set(gca,'xscale','log')
    set(gca,'yscale','log');
    xlabel('Frequency (Hz)');
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    gcaformat;
    set(gca,'FontSize',8)
    xticks([1,10,100,1e3])
    xticklabels([1,10,100,1e3])
    xlim([1,1e3])
    ylim([1e-4,1e2])
axes('Position',[0.6, 0.26, 0.25, 0.57])
    surf(XX,YY,0*P0-1,P0,'LineStyle','none')
    view([0,90]);
    set(gca,'xscale','log')
    set(gca,'yscale','log');
    set(gca,'CLim',[-10,20])
    CM = flip(bone(1e3)); CM = CM(1:800,:);
    colormap(CM);
    grid off
    xlabel('Oscillation peak frequency (Hz)')
    ylabel('\lambdaR_{max}')
    ylim([1e-4,10])
    yticks([1e-4,1e-2,1e0]);
    xlim([1,1e3])
    xticks([1,10,100,1000])
    xticklabels([1,10,100,1000])
    C = colorbar;
    C.Label.String = 'Rel. apEEG power (dB)';
    C.Position(1) = 0.89;
    gcaformat(gca,true,8);

    line([1,1],[1e-4,1],'color','k','lineWidth',0.75);
    line([1,1e3],[1e-4,1e-4],'color','k','lineWidth',0.75);