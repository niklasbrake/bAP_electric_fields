load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\unitarySpectrum.mat','Rxx')
temp = zeros(8e3,1);
temp(1:2:end) = Rxx;
temp(2:2:end) = Rxx;
Rxx = temp;

load('C:\Users\brake\Desktop\all_chains.mat');
fs =16e3;
L = 16e3;
f = fs/L:fs/L:fs/2;
mu = Pxy_D./count_D;
CI = 1.96*sqrt(SSExy_D./(count_D-1))./sqrt(count_D);

x = [f,flip(f)];
y = [mu-CI;flip(mu+CI)];
y(isinf(y(:))) = nan;


load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\cortical_area\MC_data.mat','dValues');
load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\cortical_area\pairwise_distance.mat')
N2 = 16e9;
dV = mean(diff(dValues));
DV1 = [dValues,dValues(end)+mean(diff(dValues))];
dN = interp1(0.5*(rValues(1:end-1)+rValues(2:end)),diff(total_area)./diff(rValues'),dValues,'linear','extrap');
dN = mean(dN,2)'; % d(area) / d(radius) (r)
dN = dN/mean(total_area(end,:))*16e9; % d(neurons) / d(radius) (r)

sigx2 = 3;
A = exp(-dValues.^2/(2*sigx2)).*dN*dV/N2;

Rxy = smooth(nansum(mu.*A,2));
Rxy_CI = nansum(y.*A.^2,2);


[f0,S] = import_Scheer2006;

low_noise = (8e-3).^2;
lam = 1.75;%8.5;
R = 0.2;

red = clrsPT.qualitative_CM.red;

%{
figureNB(8.5,4)
axes('Position',[0.09, 0.19, 0.35, 0.71])
    plot(f(2:2:end),Rxx(2:2:end),'k','Linewidth',1);
    xlim([10,1e3])
    set(gca,'xscale','log')
    xlabel('Frequency (Hz)')
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    xlim([1,1e3])
    xticks([1,10,100,1000]);
    xticklabels([1,10,100,1000]);
    gcaformat;
axes('Position',[0.61, 0.19, 0.35, 0.71])
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

lam = 1;%8.5;
figureNB(9,4.4);
axes('Position',[0.1, 0.21, 0.35, 0.70])
    plot(f,Rxy,'k','Linewidth',1);
    hold on;
    fill(x,Rxy_CI,'k','FaceAlpha',0.2,'EdgeColor','none');;
    set(gca,'xscale','log')
    xlim([1,1e3])
    xticks([1,10,100,1000]);
    xticklabels([1,10,100,1000]);

    set(gca,'xscale','log')
    xlabel('Frequency (Hz)')
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    gcaformat;
axes('Position',[0.60, 0.21, 0.35, 0.70])
    plot(f0,S,'color',0*[0.6,0.6,0.6],'Linewidth',1);
    hold on;

    %N Noise threshold (doi.org/10.1088/0967-3334/27/2/002)
    plot([1,3e3],low_noise*[1,1],'color',red,'LineWidth',1.5,'LineStyle','-')

    plot(f,lam*N2*Rxx,'color',clrs(1,:),'LineWidth',1)

    dR = [0.001,0.01,0.1,1];
    for i = 1:length(dR)
        plot(f,lam*N2*Rxx+dR(i)*lam*N2*(N2-1)*Rxy,'color',clrs(i,:),'LineWidth',1)
        hold on
    end
    set(gca,'xscale','log')
    set(gca,'yscale','log');
    xlabel('Frequency (Hz)');
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    xlim([1,3e3])
    ylim([1e-6,1e1]);
    xticks([1,10,100,1000]);
    xticklabels([1,10,100,1000]);

gcaformat(gcf,true,8);
%}


lam = 1;%8.5;


figureNB(6.9,5.7);
    plot(f0,S,'color',0*[0.6,0.6,0.6],'Linewidth',1);
    hold on;

    %N Noise threshold (doi.org/10.1088/0967-3334/27/2/002)
    plot([1,3e3],low_noise*[1,1],'color',red,'LineWidth',1.5,'LineStyle','-')

    % plot(f,lam*N2*Rxx,'color',[0.6,0.6,0.6],'LineWidth',1)
    % plot(f,lam*N2*Rxx + 0.2*lam*N2*(N2-1)*Rxy,'color',[0,0,0],'LineWidth',1)

    blue = [17,82,185]/255;
    dS = [0,1e-3*2.^(0:5),Inf];
    dS = [0,2,5,10,25,50,Inf]*1e-3;
    for i = 1:length(dS)
        sig = dS(i);
        B = exp(-2*(pi*f(:)*sig).^2);
        plot(f,lam*N2*Rxx + R*lam*N2*(N2-1)*B.*Rxy,'color',blue,'LineWidth',1)
    end


    % sig = 50e-3;
    % B = exp(-2*(pi*f(:)*sig).^2);
    % plot(f,lam*N2*Rxx+ 0.1*lam*N2*(N2-1)*B.*Rxy,'color',blue,'LineWidth',1)

    set(gca,'xscale','log')
    set(gca,'yscale','log');
    xlabel('Frequency (Hz)');
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    xlim([1,3e3])
    ylim([1e-7,1e1]);
    xticks([1,10,100,1000]);
    xticklabels([1,10,100,1000]);
    yticks([1e-6,1e-4,1e-2,1e0])
    gcaformat(gca,true,8);




% low_noise = 2*low_noise;
low_noise = 1e-4;
PN = @(A,lam,sig) lam*N2*Rxx + A*lam*N2*(N2-1)*exp(-(2*pi*f(:)*sig).^2).*Rxy;

lam = 10.^linspace(-1,2,50);
sig = 10.^linspace(-3,-1,100)*1e3;
A = 10.^linspace(-3,0,5);
M = zeros(length(lam),length(sig));
[XX,YY] = meshgrid(sig,lam);

w = 0.125;
red = clrsPT.qualitative_CM.red;
figureNB(13.2,2.8);
for i = 1:length(A)
    axes('Position',[0.07+(0.91-w-0.07)*(i-1)/4,0.3,0.13,0.6]);
    for j = 1:lengsig = 10e-3;
B = exp(-2*(pi*(f(:)-130)*sig).^2);
G = lam*N2*Rxx + 0.1*N2*(N2-1)*B.*Rxy;


Gsyn = exp(-(f-50).^2/10^2);

figureNB(9,6);
Psyn = 10.^AP_model(f,[20e-3,4e-3,-Inf,3.6])+10.^AP_model(f,[3e-3,1e-3,-Inf,3]);
Y1 = (1+3*A+Gsyn+1./f.^2).*Psyn+0*G+additive_noise;
Psyn2 = 10.^AP_model(f,[35e-3,4e-3,-Inf,4])+0.1*10.^AP_model(f,[3e-3,1e-3,-Inf,3]);
Y2 = (1+3*A+Gsyn+1./f.^2).*Psyn2+0*G+0.5*additive_noise;
subplot(2,2,1);
    plot(f,Psyn,'color','r');
    hold on;
    plot(f,additive_noise+f*0,'color','r');
    yticks([1e-3,1e-1,1e1])
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    plot(f,Y1,'color',[0.6,0.6,0.6],'LineWidth',1);
    ylim([1e-3,1e1]);
    set(gca,'yscale','log')
    xlim([0,200])
    xlabel('Frequency (Hz)');
subplot(2,2,2);
    plot(f,Psyn2,'color','r');
    hold on;
    plot(f,0.5*additive_noise+f*0,'color','r');
    yticks([1e-3,1e-1,1e1])
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    plot(f,Y1,'color',[0.6,0.6,0.6],'LineWidth',1);
    plot(f,Y2,'k','LineWidth',1);
    % set(gca,'xscale','log')
    ylim([1e-3,1e1]);
    set(gca,'yscale','log')
    xlim([0,200])
    xlabel('Frequency (Hz)');

Z1 = Y1;
Z2 = Y2;
subplot(2,3,4);
    fill([f;f(end);.1],[10*log10(Z2./Z1);0;0],'k','LineWidth',1,'FaceAlpha',0.4);
    xlim([1,200])
    xlabel('Frequency (Hz)');
    set(gca,'xscale','log');
    xticks([1,10,100])
    xticklabels([1,10,100])

Z1 = Y1./(Psyn+additive_noise);
Z2 = Y2./(Psyn2+0.5*additive_noise);
subplot(2,3,5);
    fill([f;f(end);.1],[10*log10(Z2./Z1);0;0],'k','LineWidth',1,'FaceAlpha',0.4);
    % hold on;
    % plot(f,Z2,'color','k','LineWidth',1);
    xlim([1,200])
    xlabel('Frequency (Hz)');
    set(gca,'xscale','log');
    xticks([1,10,100])
    xticklabels([1,10,100])

Z1 = (Y1-additive_noise)./Psyn;
Z2 = (Y2-0.5*additive_noise)./Psyn2;
subplot(2,3,6);
    fill([f;f(end);.1],[10*log10(Z2./Z1);0;0],'k','LineWidth',1,'FaceAlpha',0.4);
    % hold on;
    % plot(f,Z2,'color','k','LineWidth',1);
    xlim([1,200])
    xlabel('Frequency (Hz)');
    set(gca,'xscale','log');
    xticks([1,10,100])
    xticklabels([1,10,100])


sig = 10e-3;
B = exp(-2*(pi*(f(:)-130)*sig).^2);
G = lam*N2*Rxx + 0.1*N2*(N2-1)*B.*Rxy;


Gsyn = exp(-(f-50).^2/10^2);

figureNB(9,6);
Psyn = 10.^AP_model(f,[20e-3,4e-3,-Inf,3.6])+10.^AP_model(f,[3e-3,1e-3,-Inf,3]);
Y1 = (1+3*A+Gsyn+1./f.^2).*Psyn+2*G+additive_noise;
Psyn2 = 10.^AP_model(f,[35e-3,4e-3,-Inf,4])+0.1*10.^AP_model(f,[3e-3,1e-3,-Inf,3]);
Y2 = (1+3*A+Gsyn+1./f.^2).*Psyn2+2*G+0.5*additive_noise;
subplot(2,2,1);
    plot(f,Psyn,'color','r');
    hold on;
    plot(f,additive_noise+f*0,'color','r');
    yticks([1e-3,1e-1,1e1])
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    plot(f,Y1,'color',[0.6,0.6,0.6],'LineWidth',1);
    ylim([1e-3,1e1]);
    set(gca,'yscale','log')
    xlim([0,200])
    xlabel('Frequency (Hz)');
subplot(2,2,2);
    plot(f,Psyn2,'color','r');
    hold on;
    plot(f,0.5*additive_noise+f*0,'color','r');
    yticks([1e-3,1e-1,1e1])
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    plot(f,Y1,'color',[0.6,0.6,0.6],'LineWidth',1);
    plot(f,Y2,'k','LineWidth',1);
    % set(gca,'xscale','log')
    ylim([1e-3,1e1]);
    set(gca,'yscale','log')
    xlim([0,200])
    xlabel('Frequency (Hz)');

Z1 = Y1;
Z2 = Y2;
subplot(2,3,4);
    fill([f;f(end);.1],[10*log10(Z2./Z1);0;0],'k','LineWidth',1,'FaceAlpha',0.4);
    xlim([1,200])
    xlabel('Frequency (Hz)');
    set(gca,'xscale','log');
    xticks([1,10,100])
    xticklabels([1,10,100])

Z1 = Y1./(Psyn+additive_noise);
Z2 = Y2./(Psyn2+0.5*additive_noise);
subplot(2,3,5);
    fill([f;f(end);.1],[10*log10(Z2./Z1);0;0],'k','LineWidth',1,'FaceAlpha',0.4);
    % hold on;
    % plot(f,Z2,'color','k','LineWidth',1);
    xlim([1,200])
    xlabel('Frequency (Hz)');
    set(gca,'xscale','log');
    xticks([1,10,100])
    xticklabels([1,10,100])

Z1 = (Y1-additive_noise)./Psyn;
Z2 = (Y2-0.5*additive_noise)./Psyn2;
subplot(2,3,6);
    fill([f;f(end);.1],[10*log10(Z2./Z1);0;0],'k','LineWidth',1,'FaceAlpha',0.4);
    % hold on;
    % plot(f,Z2,'color','k','LineWidth',1);
    xlim([1,200])
    xlabel('Frequency (Hz)');
    set(gca,'xscale','log');
    xticks([1,10,100])
    xticklabels([1,10,100])
th(lam)
        for k = 1:length(sig)
            P = PN(A(i),lam(j),1e-3*sig(k));
            [M(j,k),I(j,k)] = max(P(f>30));
        end
    end
    % imagesc(sig*1e3,lam,log10(M));
    surf(XX,YY,0*M,log10(M),'LineStyle','none')
    view([0,90]);
    hold on;
    C = contour(XX,YY,log10(M),log10([low_noise,low_noise]),'color',red,'linewidth',1);
    % set(gca,'CLim',[log10(low_noise),2])
    if(i==4)
        fill([10,90,90,10],[0.11,0.11,4,4],'k','LineStyle','--','FaceColor','none');
    end
    set(gca,'CLim',[-4,2]);
    ylim(10.^[-1,2])
    xlim(10.^[0,2])
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    line(1e3*[1e-3,0.1],[100,100],[0,0]-1e6,'color','k','LineWidth',0.75)
    line(1e3*[0.1,0.1],[0.1,100],[0,0]-1e6,'color','k','LineWidth',0.75)

    xticks([1,10,100])
    xticklabels([1,10,100])
    yticks([0.1,1,10,100])
    yticklabels([0.1,1,10,100])
    gcaformat
    xlabel('Jitter (ms)','FontSize',8)
    if(i==1)
        ylabel('Firing rate (Hz)','FontSize',8)
    end
    set(gca,'FontSize',8);
    colormap(flip(bone))
    drawnow;
end
C = colorbar;
C.Position = [0.92,0.3,0.01,0.6];
C.Ticks = [-4:2:2];
C.Label.Position = [5,-0.8,0];
C.TickLabels = {'10^{-4}','10^{-2}','10^{0}','10^{2}'};
C.Label.String = ['Max PSD (' char(956) 'V^2/Hz)'];






low_noise = 2*low_noise;
PN = @(A,lam,sig) lam*N2*Rxx + A*lam*N2*(N2-1)*exp(-(2*pi*f(:)*sig).^2).*Rxy;

lam = 10.^linspace(-1,2,50);
sig = 10.^linspace(-3,-1,100)*1e3;
A = 10.^linspace(-3,0,5);
M = zeros(length(lam),length(sig));
[XX,YY] = meshgrid(sig,lam);

[~,I] = unique(f0);
S0 = interp1(f0(I),S(I),f(f<=30));

w = 0.125;
red = clrsPT.qualitative_CM.red;
figureNB(13.2,2.8);
for i = 1:length(A)
    axes('Position',[0.07+(0.91-w-0.07)*(i-1)/4,0.3,0.13,0.6]);
    for j = 1:length(lam)
        for k = 1:length(sig)
            P = PN(A(i),lam(j),1e-3*sig(k));
            M(j,k) = max(10*log10(P(f<=30)./S0(:)));
        end
    end
    % imagesc(sig*1e3,lam,log10(M));
    surf(XX,YY,0*M,M,'LineStyle','none')
    view([0,90]);
    hold on;
    [C,h] = contour(XX,YY,M,[-20,-10,0],'color','k','linewidth',1);
    if(i==4)
        fill([10,90,90,10],[0.11,0.11,4,4],'k','LineStyle','--','FaceColor','none');
    end
    set(gca,'CLim',[-20,20])
    ylim(10.^[-1,2])
    xlim(10.^[0,2])
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    line(1e3*[1e-3,0.1],[100,100],[0,0]-1e6,'color','k','LineWidth',0.75)
    line(1e3*[0.1,0.1],[0.1,100],[0,0]-1e6,'color','k','LineWidth',0.75)

    xticks([1,10,100])
    xticklabels([1,10,100])
    yticks([0.1,1,10,100])
    yticklabels([0.1,1,10,100])
    gcaformat
    xlabel('Jitter (ms)','FontSize',8)
    if(i==1)
        ylabel('Firing rate (Hz)','FontSize',8)
    end
    set(gca,'FontSize',8);
    colormap(flip(bone))
    drawnow;
end
C = colorbar;
C.Position = [0.92,0.3,0.01,0.6];
% C.Ticks = [-4:2:2];
% C.Label.Position = [5,-0.8,0];
% C.TickLabels = {'10^{-4}','10^{-2}','10^{0}','10^{2}'};
C.Label.String = ['Rel. power (dB)'];

