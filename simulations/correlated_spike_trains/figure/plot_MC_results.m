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

sigx2 = 5;

Rxy = nansum(mu.*exp(-dValues.^2/(2*sigx2)).*dN*dV,2)/N2;
Rxy_CI = nansum(y.*exp(-dValues.^2/(2*sigx2)).*dN*dV,2)/N2;


[f0,S] = import_Scheer2006;

low_noise = (8e-3).^2;
lam = 1;%8.5;
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


figureNB(6.5,5.5);
    plot(f0,S,'color',0*[0.6,0.6,0.6],'Linewidth',1);
    hold on;

    %N Noise threshold (doi.org/10.1088/0967-3334/27/2/002)
    plot([1,3e3],low_noise*[1,1],'color',red,'LineWidth',1.5,'LineStyle','-')

    % plot(f,lam*N2*Rxx,'color',[0.6,0.6,0.6],'LineWidth',1)
    % plot(f,lam*N2*Rxx + 0.2*lam*N2*(N2-1)*Rxy,'color',[0,0,0],'LineWidth',1)

    dS = [0,1e-3*2.^(0:5),Inf];
    for i = 1:length(dS)
        sig = dS(i);
        B = exp(-2*(pi*f(:)*sig).^2);
        if(i==1 || i==length(dS))
            plot(f,lam*N2*Rxx + R*lam*N2*(N2-1)*B.*Rxy,'color',0.6+[0,0,0],'LineWidth',1)
        else
            plot(f,lam*N2*Rxx + R*lam*N2*(N2-1)*B.*Rxy,'color',[0.6,0.6,0.6],'LineWidth',1)
        end
    end

    blue = [17,82,185]/255;


    sig = 8e-3*sqrt(2);
    B = exp(-2*(pi*f(:)*sig).^2);
    plot(f,lam*N2*Rxx+ 0.2*lam*N2*(N2-1)*B.*Rxy,'color',blue,'LineWidth',1)

    % sig = 50e-3;
    % B = exp(-2*(pi*f(:)*sig).^2);
    % plot(f,lam*N2*Rxx+ 0.1*lam*N2*(N2-1)*B.*Rxy,'color',blue,'LineWidth',1)

    set(gca,'xscale','log')
    set(gca,'yscale','log');
    xlabel('Frequency (Hz)');
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    xlim([1,3e3])
    ylim([1e-6,1e1]);
    xticks([1,10,100,1000]);
    xticklabels([1,10,100,1000]);
    gcaformat(gca,true,8);





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
    for j = 1:length(lam)
        for k = 1:length(sig)
            P = PN(A(i),lam(j),1e-3*sig(k));
            [M(j,k),I(j,k)] = max(P(f>40));
        end
    end
    % imagesc(sig*1e3,lam,log10(M));
    surf(XX,YY,0*M,log10(M),'LineStyle','none')
    view([0,90]);
    hold on;
    C = contour(XX,YY,log10(M),log10([low_noise,low_noise]),'color',red,'linewidth',1);
    % set(gca,'CLim',[log10(low_noise),2])
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



