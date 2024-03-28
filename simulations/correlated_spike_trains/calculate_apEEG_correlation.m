load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\spike_synchrony\cross_correlation_uapEEG.mat','R_eeg');

dt = 1/16e3;

A = linspace(0,0.5,200);
S = linspace(1e-3,20e-3,200);
t = dt*(-1e3:1e3)';
baseline_rate = 10*dt; % 10 Hz
spike_xcorr = @(a,s) a*exp(-t.^2/(2*2*s^2))/sqrt(2*pi*2*s^2);

R = mean(R_eeg,2);

for i = 1:length(A)
    for j = 1:length(S)
        a = A(i);
        s = S(j);
        expected_corr(i,j) = sum(R.*spike_xcorr(a,s)*dt);
    end
end

s0 = sqrt(2)*8e-3;
t = 1/16*(-1e3:1e3);

figureNB(9,5);
axes('Position',[0.11, 0.20, 0.30, 0.72])
    plotwitherror(t,R_eeg,'CI','color','k')
    xlim([-10,10])
    ylim([-0.1,0.28])
    xlabel('Lag (ms)')
    ylabel('uAP correlation')
axes('Position',[0.55, 0.20, 0.32, 0.73])
    imagesc(S*1e3,A,log10(abs(expected_corr)))
    hold on
    plot(s0*1e3,0.2,'xr','MarkerSize',10,'LineWidth',2)
    text(10,0.22,sprintf('%.1d',sum(R.*spike_xcorr(0.2,s0))),'color','k')
    axis xy;
    set(gca,'Clim',[-4,-1]);
    ylabel('R_{noise}')
    xlabel('Timescale (ms)')
    colormap(flip(bone))
    C = colorbar;
    xlabel('Jitter (ms)')
    ylabel('Spike time correlation')
    C.Label.String = 'apEEG correlation'
    C.TickLabels = cellfun(@(x) ['10^{' x '}'],C.TickLabels,'UniformOutput',false);
