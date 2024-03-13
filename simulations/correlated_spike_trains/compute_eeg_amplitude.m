% Get average correlation due to cortex geometry, assuming spatial decay in correlations
dataFolder = 'E:\Research_Projects\004_Propofol\manuscript\Nature Communications\final_submission\manuscript_source_data';
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
t = dt*(-2e3:2e3)';
rho_max = @(a,s) sum(R.*a.*exp(-t.^2/(2*2*s^2))/sqrt(2*pi*2*s^2)*dt);

% Get simulated passive spectrum
model = load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\unitary_AP_PSD.mat');
ap_unit = mean(model.bootstrap_Est,2);
P0 = interp1(model.freq,ap_unit,data.freq);

pars.noise_correlation = 0*0.2;
pars.jitter_window = 100e-3; % 8 ms
pars.spatial_scale = 3.5; % 5 mm
pars.firing_rate = 8.5; % Hz

scaled_uAP = @(P) P.firing_rate*P0*scaling(rho_max(P.noise_correlation,P.jitter_window),P.spatial_scale);

% Get baseline EEG spectrum from propofol cohort
data = load(fullfile(dataFolder,'EEG_data','electrode2_Cz.mat'));
data.baseline = squeeze(nanmedian(data.psd(:,data.time>0,:),2));

figureNB;
    plotwitherror(data.freq,data.baseline,'M','LineWidth',1,'color',[0.6,0.6,0.6]);
    pars.noise_correlation = 0;
    plot(data.freq,scaled_uAP(pars),'color','k','LineWidth',1)
    pars.noise_correlation = 0.05;
    pars.jitter_window = 10e-3;
    plot(data.freq,scaled_uAP(pars),'color','r','LineWidth',1)
    pars.noise_correlation = 0.2;
    pars.jitter_window = 8e-3;
    plot(data.freq,scaled_uAP(pars),'color','b','LineWidth',1) 

    set(gca,'yscale','log');
    xlim([1,150])
    xticks([1,50,100,150])
    xlabel('Frequency (Hz)');
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    gcaformat;