[full_model,AP_model] = fittingmodel('eq6');

[f0,S] = import_paralytic_data;

low_noise = 1e-3;


figureNB(13.2,5.8);
axes('Position',[0.08, 0.14, 0.38, 0.8])
    plot(f0,S(:,9:16),'color','k');
    hold on;
    plot([1,3e3],low_noise*[1,1],'k','LineWidth',1,'LineStyle','--')
    Psyn = 10.^AP_model(f0,[3e-3,1e-3,-Inf,4.6]); 
    plot(f0,Psyn,'color',blue,'LineStyle','--')
    Psyn = 10.^AP_model(f0,[20e-3,4e-3,-Inf,3.6]); 
    plot(f0,Psyn,'color',blue,'LineStyle','--')
    Psyn = 10.^AP_model(f0,[20e-3,4e-3,-Inf,3.6])+10.^AP_model(f0,[3e-3,1e-3,-Inf,4.6]); 
    plot(f0,Psyn,'color',blue,'LineWidth',1)
    set(gca,'yscale','log')
    xlabel('Frequency (Hz)')
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    title('Unparalyzed','FontWeight','normal')
    xlim([0,100]);
    gcaformat;
    set(gca,'FontSize',8)

    % sig = 8e-3*sqrt(2);
    % B = exp(-(2*pi*freq(:)*sig).^2/2);
    % plot(freq,lam*N2*Rxx2+ R*lam*N2*(N2-1)*B.*Rxy2,'color',blue,'LineWidth',1)
axes('Position',[0.58, 0.14, 0.38, 0.8])
    plot(f0,S(:,1:8),'color','k');
    hold on;
    plot([1,3e3],low_noise*[1,1],'k','LineWidth',1,'LineStyle','--')
    Psyn = 10.^AP_model(f0,[3e-3,1e-3,-Inf,3.3]); 
    plot(f0,Psyn,'color',blue,'LineStyle','--')
    Psyn = 10.^AP_model(f0,[20e-3,4e-3,-Inf,3.6]); 
    plot(f0,Psyn,'color',blue,'LineStyle','--')
    Psyn = 10.^AP_model(f0,[20e-3,4e-3,-Inf,3.6])+10.^AP_model(f0,[3e-3,1e-3,-Inf,3.3]); 
    plot(f0,Psyn,'color',blue,'LineWidth',1)
    set(gca,'yscale','log')
    xlabel('Frequency (Hz)')
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    title('20 mg cisatracurium','FontWeight','normal')
    xlim([0,100]);
    gcaformat;
    set(gca,'FontSize',8)

    % sig = 8e-3*sqrt(2);
    % B = exp(-(2*pi*freq(:)*sig).^2/2);
    % plot(freq,lam*N2*Rxx2+ R*lam*N2*(N2-1)*B.*Rxy2,'color',blue,'LineWidth',1)
