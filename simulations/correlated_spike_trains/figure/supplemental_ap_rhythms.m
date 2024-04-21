Rxx = Pxx/N;
Rxy = Pxy/N;

CI = 1.96*sqrt(SSE_xy/(N-1))/sqrt(N);

x = [f,flip(f)];
y = [Rxy-CI;flip(Rxy+CI)];
y(isinf(y(:))) = nan;


lam = 8.5;
figureNB(9,4.4);
    plot(f,Rxy,'k','Linewidth',1);
    hold on;
    fill(x,y,'k','FaceAlpha',0.2,'EdgeColor','none');;
    set(gca,'xscale','log')
    xlim([1,1e3])
    xticks([1,10,100,1000]);
    xticklabels([1,10,100,1000]);

    set(gca,'xscale','log')
    xlabel('Frequency (Hz)')
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    gcaformat;



[f0,S] = import_paralytic_data;

R = 0.2;
lam = 5;

dF = 5*2.^[0:3];
sig = 8e-3*sqrt(2);
% sig = 1e2;

clrs = clrsPT.sequential(length(dF)+3);
clrs = clrs(4:end,:);
blue = [17,82,185]/255;

figureNB(11,9);

subplot(2,2,2);
t = linspace(-60,60,1e3);
plot(t,0.04*exp(-1e-6*t.^2./2/sig.^2).*cos(2*pi*t*1e-3*40),'color',blue,'LineWidth',1)
xlabel('Lag (ms)')
ylabel('Correlation')
gcaformat(gca,true,8)

subplot(2,2,3);
i = 4;
    plot(f0,S(:,1:8),'k','Linewidth',1);
    hold on;

    plot([1,3e3],low_noise*[1,1],'r','LineWidth',1.5,'LineStyle','-')

    F0 = dF(i);
    B = exp(-2*(pi*(f(:)-F0)*sig).^2);
    plot(f,lam*N2*Rxx+R*lam*N2*(N2-1)*Rxy,'color',[0.6,0.6,0.6],'LineWidth',1,'LineStyle',':')
    plot(f,lam*N2*Rxx+0*R*lam*N2*(N2-1)*Rxy,'color',[0.6,0.6,0.6],'LineWidth',1,'LineStyle',':')

    plot(f,lam*N2*Rxx + R*lam*N2*(N2-1)*B.*Rxy,'color',blue,'LineWidth',1)

    set(gca,'yscale','log');
    xlabel('Frequency (Hz)');
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    ylim([1e-4,1e2])
    yticks([1e-4,1e-2,1e0,1e2])
    xlim([0,100]);

subplot(2,2,4);
for i = 1:4
    F0 = dF(i);
    B = exp(-2*(pi*(f(:)-F0)*sig).^2);
    mdl = interp1(f,lam*N2*Rxx + R*lam*N2*(N2-1)*B.*Rxy,f0(:));
    plot(f0,mdl./nanmean(S(:,1:8),2),'color',clrs(i,:),'LineWidth',1)
    hold on;
end


xlabel('Frequency (Hz)');
ylabel('Fraction apEEG')
xlim([0,100])
C = colorbar;
colormap(clrs(1:4,:));
C.Label.String = 'Synchronization frequency (Hz)';
C.Ticks = [1/9:1/4:1-1/9];
C.TickLabels = dF(1:4);

gcaformat(gcf,true,8);

