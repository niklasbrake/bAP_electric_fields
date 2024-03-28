
[sa,X] = network_simulation_beluga.getHeadModel;
load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\unitaryAP_all.mat')
uAP = network_simulation_beluga.getEEG(savedUnitaryAP,sa,43e3);

dt = 1/16e3;
tmax = 1;
T = tmax/dt+1;
t = 0:dt:tmax;
N = 100;

m1 = 20;%floor(tmax*100);

X1 = zeros(T,N);
X2 = zeros(T,N);
Y1 = zeros(T,N);
Y2 = zeros(T,N);
t0 = (tmax-250e-3)*rand(m1,1)+125e-3;

A = 1;
S = 1e-4;


I = sample_blue_neurons(N);
y = uAP(:,I);

% Caluclate cross power spectrum
fs = 16e3;
R_eeg = [];
L = 4e3;
f0 = fs/L:fs/L:fs/2;
sxy = zeros(length(f0),N*(N-1)/2);
sxx = zeros(length(f0),N);
count = 1;
for i=1:N
    R_eeg = xcorr(y(:,i),y(:,i));
    R_eeg = R_eeg(2:end);
    S11 = fft(R_eeg)/fs;
    S11 = 2*S11(2:2001).*exp(2*pi*sqrt(-1)*f0(:)*1999/fs);
    sxx(:,i) = real(S11);
    for j = i+1:N
        R_eeg = xcorr(y(:,i),y(:,j));
        R_eeg = R_eeg(2:end);
        S12 = fft(R_eeg)/fs;
        S12 = 2*S12(2:2001).*exp(2*pi*sqrt(-1)*f0(:)*1999/fs);

        R_eeg = xcorr(y(:,j),y(:,i));
        R_eeg = R_eeg(2:end);
        S21 = fft(R_eeg)/fs;
        S21 = 2*S21(2:2001).*exp(2*pi*sqrt(-1)*f0(:)*1999/fs);

        sxy(:,count) = S12+S21;
        count = count + 1;
    end
end
sxy = A*round(A*m1)*sxy;
sxx = round(A*m1)*sxx;
%%%%

for i = 1:N
    % m2 = poissrnd(A*m1);
    m2 = round(A*m1);
    s1 = t0(randperm(m1,m2))+S*randn(m2,1);
    s2 = tmax*rand(m2,1);
    
    it = interp1(t,1:T,s1,'nearest','extrap');
    [I,J] = findgroups(it);
    X1(J,i) = splitapply(@length,I,I);

    it = interp1(t,1:T,s2,'nearest','extrap');
    [I,J] = findgroups(it);
    X2(J,i) = splitapply(@length,I,I);
    
    Y1(:,i) = filter(y(:,i),1,X1(:,i));
    Y2(:,i) = filter(y(:,i),1,X2(:,i));
end

Y1 = Y1(1001:end,:);
Y2 = Y2(1001:end,:);
[psdN,f] = pmtm(sum(Y1,2),2,[],16e3);
[psdi,f] = pmtm(Y1,2,[],16e3);
[psd2,f] = pmtm(sum(Y2,2),2,[],16e3);


psd0 = sum(psdi,2);
sxy0 = interp1(f0,sum(sxy,2),f)/size(Y1,1);
sxx0 = interp1(f0,sum(sxx,2),f)/size(Y1,1);

% B = N*(N-1)*mean(cos(2*pi*f(:)*S.*randn(1,N)),2);
B = exp(-(2*pi*f*S*sqrt(2)).^2/2);

fig = figureNB(13.5,5.5);
subplot(1,2,1);
    plot(f,psd0,'color',[0.6,0.6,0.6]);
    hold on;
    plot(f,sxx0,'color','k','LineWidth',1,'LineStyle','--');
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlim([5,2e3])
    ylim([1e-17,1e-11])
    yticks([1e-16,1e-14,1e-12])
    xticks([10,100,1000])
    xticklabels([10,100,1000])
    xlabel('Frequency (Hz)');
    ylabel('Power');
    grid on
    set(gca,'MinorGridLineStyle','-')
    set(gca,'MinorGridAlpha',0.1)
    set(gca,'GridAlpha',0.3);
    set(gca,'LineWidth',0.75);
    set(gca,'FontSize',10);
    L = legend('Simulation','Theory');
    L.ItemTokenSize = [20,10]; L.FontSize = 8;
    set(gca,'FontName','LM Roman 10')
subplot(1,2,2);
    plot(f,psdN,'color',[0.6,0.6,0.6]);
    hold on;
    plot(f,sxx0+B.*sxy0,'color','k','LineWidth',1,'LineStyle','--')
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlim([5,2e3])
    ylim([1e-17,1e-11])
    yticks([1e-16,1e-14,1e-12])
    xticks([10,100,1000])
    xticklabels([10,100,1000])
    xlabel('Frequency (Hz)');
    ylabel('Power');
    grid on
    set(gca,'MinorGridLineStyle','-')
    set(gca,'MinorGridAlpha',0.1)
    set(gca,'GridAlpha',0.3);
    set(gca,'LineWidth',0.75);
    set(gca,'FontSize',10);
    set(gca,'FontName','LM Roman 10')

set(fig,'PaperPosition',[0,0,fig.Position(3:4)],'PaperUnits','centimeters','PaperSize',fig.Position(3:4),'Renderer','Painters');
print('E:/Research_Projects/005_Aperiodic_EEG/unitary_APs/manuscript/figures/workspace/figure4/theory_cross_spectrum.png','-dpng','-r600');