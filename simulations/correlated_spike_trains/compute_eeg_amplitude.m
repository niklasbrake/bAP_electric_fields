% Get baseline EEG spectrum from propofol cohort
dataFolder = 'E:\Research_Projects\004_Propofol\manuscript\Nature Communications\final_submission\manuscript_source_data';
data = load(fullfile(dataFolder,'EEG_data','electrode2_Cz.mat'));
data.baseline = squeeze(nanmedian(data.psd(:,data.time<-1,:),2));

% Get average correlation due to cortex geometry, assuming spatial decay in correlations
load(fullfile(dataFolder,'cortex_anatomy','anatomy_cortical_pairwise_distance_distribution.mat'));
signed_area = A;
N = 16e9;
dMids = 0.5*(rValues(2:end)+rValues(1:end-1));
nrnCount = mean(diff(signed_area),2)*200000;
nrnCount(end) = N-sum(nrnCount(1:end-1));
scaling = @(rho_max,sigma) N+N*(N-1)*rho_max*sum(exp(-dMids.^2/sigma).*nrnCount)/sum(nrnCount');

% Compute maximal correlation in apEEG from spike synchrony parameters
load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\spike_synchrony\cross_correlation_uapEEG.mat','R_eeg');
R = mean(R_eeg,2);
dt = 1/16e3;
t = dt*(-1e3:1e3)';
rho_max = @(a,s) sum(R.*a.*exp(-t.^2/(2*2*s^2))/sqrt(2*pi*2*s^2)*dt);

% Get simulated passive spectrum
model = load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\unitary_AP_PSD.mat');
ap_unit = mean(model.bootstrap_Est,2);
% P0 = interp1(model.freq,ap_unit,data.freq);

pars.noise_correlation = 0*0.2;
pars.jitter_window = 100e-3; % 8 ms
pars.spatial_scale = 3.5; % 5 mm
pars.firing_rate = 8.5; % Hz

scaled_uAP = @(P) P.firing_rate*ap_unit*scaling(rho_max(P.noise_correlation,P.jitter_window),P.spatial_scale);


figureNB;
    pars.noise_correlation = 0;
    plot(model.freq,scaled_uAP(pars),'color','k','LineWidth',1);
    hold on;
    pars.noise_correlation = 0.2;
    pars.jitter_window = 8e-3*sqrt(2);
    plot(model.freq,scaled_uAP(pars),'color','b','LineWidth',1)

    pars.noise_correlation = 0.02;
    pars.jitter_window = sqrt(10e-3);
    plot(model.freq,scaled_uAP(pars),'color','r','LineWidth',1)

    %N Noise threshold (doi.org/10.1088/0967-3334/27/2/002)
    low_noise = (8/1e3)^2;;
    plot([1,2e3],low_noise*[1,1],'k','LineWidth',1,'LineStyle','--')


    set(gca,'xscale','log')
    set(gca,'yscale','log');
    xlim([1,3e3])
    xticks([1,10,100,1e3,3e3])
    xticklabels({1,10,100,'1k','3k'});
    ytk = ([1,10,100,1e3,1e4]/1e3).^2;
    yticks(ytk)
    ylim([ytk(1),ytk(end)]);
    % yticklabels({1,10,100,'1k','10k'});

    xlabel('Frequency (Hz)');
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    gcaformat;

return;
figureNB(10,5.8);
axes('Units','centimeters','Position',[1.38,0.94,8.21,4.52]);
    pars.noise_correlation = 0;
    plot(model.freq,1e3*sqrt(scaled_uAP(pars)),'color','k','LineWidth',1);
    hold on;
    pars.noise_correlation = 0.05;
    pars.jitter_window = 10e-3;
    plot(model.freq,1e3*sqrt(scaled_uAP(pars)),'color','r','LineWidth',1)
    pars.noise_correlation = 0.2;
    pars.jitter_window = 8e-3;
    plot(model.freq,1e3*sqrt(scaled_uAP(pars)),'color','b','LineWidth',1)

    %N Noise threshold (doi.org/10.1088/0967-3334/27/2/002)
    low_noise = 8;
    plot([1,3e3],low_noise*[1,1],'k','LineWidth',1,'LineStyle','--')


    set(gca,'xscale','log')
    set(gca,'yscale','log');
    xlim([1,3e3])
    xticks([1,10,100,1e3,3e3])
    xticklabels({1,10,100,'1k','3k'});
    ytk = [1,10,100,1e3,1e4];
    yticks(ytk)
    ylim([ytk(1),ytk(end)]);
    yticklabels({1,10,100,'1k','10k'});

    xlabel('Frequency (Hz)');
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    grid on;
    set(gca,'MinorGridLineStyle','-')
    set(gca,'MinorGridLineWidth',0.2)
    set(gca,'LineWidth',0.2)
    set(gca,'TickDir','in')

figureNB;
    plotwitherror(freq,post,'M','LineWidth',1,'color',[0.6,0.6,0.6]);
    pars.noise_correlation = 0;
    plot(model.freq,scaled_uAP(pars),'color','k','LineWidth',1)
    pars.noise_correlation = 0.05;
    pars.jitter_window = 10e-3;
    plot(model.freq,scaled_uAP(pars),'color','r','LineWidth',1)
    pars.noise_correlation = 0.2;
    pars.jitter_window = 100e-3;
    plot(model.freq,scaled_uAP(pars),'color','b','LineWidth',1)

    set(gca,'yscale','log');
    % xlim([1,150])
    % xticks([1,50,100,150])
    xlabel('Frequency (Hz)');
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    gcaformat;


scaled_uAP1 = @(P) P.firing_rate*psd*scaling(rho_max(P.noise_correlation,P.jitter_window),P.spatial_scale);
scaled_uAP2 = @(P) P.firing_rate*psd2*scaling(rho_max(P.noise_correlation,P.jitter_window),P.spatial_scale);

figureNB(10,5.8);
axes('Units','centimeters','Position',[1.38,0.94,8.21,4.52]);
    pars.noise_correlation = 0;
    plot(f,1e3*sqrt(scaled_uAP(pars)),'color','k','LineWidth',1);
    hold on;
    pars.noise_correlation = 0.2;
    pars.jitter_window = 8e-3;
    plot(f,1e3*sqrt(scaled_uAP1(pars)),'color','r','LineWidth',1)
    pars.noise_correlation = 0.2;
    pars.jitter_window = 8e-3;
    plot(f,1e3*sqrt(scaled_uAP2(pars)),'color','b','LineWidth',1)

    %N Noise threshold (doi.org/10.1088/0967-3334/27/2/002)
    low_noise = 8;
    plot([1,3e3],low_noise*[1,1],'k','LineWidth',1,'LineStyle','--')


    set(gca,'xscale','log')
    set(gca,'yscale','log');
    xlim([1,3e3])
    xticks([1,10,100,1e3,3e3])
    xticklabels({1,10,100,'1k','3k'});
    ytk = [1,10,100,1e3,1e4];
    yticks(ytk)
    ylim([ytk(1),ytk(end)]);
    yticklabels({1,10,100,'1k','10k'});

    xlabel('Frequency (Hz)');
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    grid on;
    set(gca,'MinorGridLineStyle','-')
    set(gca,'MinorGridLineWidth',0.2)
    set(gca,'LineWidth',0.2)
    set(gca,'TickDir','in')
    grid off
    xticks([])
    xlabel('')
    ylabel('')
    yticks([])
    set(gcf,'color','none')
    set(gca,'color','none')