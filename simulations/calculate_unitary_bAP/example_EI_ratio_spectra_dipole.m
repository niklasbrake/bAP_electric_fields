
EI_vec = {'01','1.5','2.1','3.1','4.5','6.6','9.7','14.1','20.6','30'};
EI = cellfun(@(x) str2num(x),EI_vec);
iCtx = 16e3;

i=50;
% i=60;
folder0 = 'E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\unitary_AP_EI_ratio\L6_BPC_cADpyr231_3\matlab_recordings';
clrs = flipud(jet(10));

meanQx = [];
meanQy = [];
meanQz = [];
N = zeros(length(EI_vec),1);
firingFrequency = zeros(size(N));
fs = 16e3;
psd_active = [];
psd_passive = [];

for k = 1:length(EI_vec)
    load(fullfile(folder0,sprintf('synaptic_input_EI%s_passive.mat',EI_vec{k})))
    for j = 1:3
        [f,~,psd] = eegfft(time*1e-3,detrend(dipoles(:,j)),0.5,0.4);
        psd_passive(:,j,k) = mean(psd,2);
    end

    dipoles_passive = dipoles;
    if(k==1)
        passiveQ = dipoles;
    end

    load(fullfile(folder0,sprintf('synaptic_input_EI%s.mat',EI_vec{k})))
    for j = 1:3
        [f,~,psd] = eegfft(time*1e-3,detrend(dipoles(:,j)),0.5,0.4);
        psd_active(:,j,k) = mean(psd,2);
    end

    if(k==1)
        fullQ = dipoles;
        tSave = time;
    end

    dipoles = dipoles-dipoles_passive;

    [y,x] = findpeaks(voltage,'MinPeakHeight',0);
    N(k) = length(x);
    firingFrequency(k) = N(k)/range(time)*1e3;
    if(N(k)<2)
        continue;
    end
    Qx = zeros(1,2001);
    Qy = zeros(1,2001);
    Qz = zeros(1,2001);
    for j = 1:length(x)
        idcs = max(min(x(j)-1e3:x(j)+1e3,length(dipoles)),1);
        y = dipoles(idcs,1);
        y(idcs==1) = 0;
        y(idcs==length(dipoles)) = 0;
        Qx(j,:) = y;

        y = dipoles(idcs,2);
        y(idcs==1) = 0;
        y(idcs==length(dipoles)) = 0;
        Qy(j,:) = y;

        y = dipoles(idcs,3);
        y(idcs==1) = 0;
        y(idcs==length(dipoles)) = 0;
        Qz(j,:) = y;
    end
    Qx = Qx-median(Qx,2);
    Qy = Qy-median(Qy,2);
    Qz = Qz-median(Qz,2);
    meanQx = [meanQx;Qx];
    meanQy = [meanQy;Qy];
    meanQz = [meanQz;Qz];
end


m = size(meanQx,1);

paddedAP = [zeros(m,7e3),meanQx,zeros(m,7e3-1)]';
paddedAP = meanQx';
paddedAP(isnan(paddedAP(:))) = 0;
paddedAP = detrend(median(paddedAP,2),'constant');

n = length(paddedAP);
xdft = fft(paddedAP);
xdft = xdft(1:n/2+1);
psdx = abs(xdft/fs).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:fs/n:fs/2;
psd_Qx = interp1(freq,psdx,f,'linear','extrap');

m = size(meanQy,1);

paddedAP = [zeros(m,7e3),meanQy,zeros(m,7e3-1)]';
paddedAP = meanQy';
paddedAP(isnan(paddedAP(:))) = 0;
paddedAP = detrend(median(paddedAP,2),'constant');

n = length(paddedAP);
xdft = fft(paddedAP);
xdft = xdft(1:n/2+1);
psdx = abs(xdft/fs).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:fs/n:fs/2;
psd_Qy = interp1(freq,psdx,f,'linear','extrap');


m = size(meanQz,1);

paddedAP = [zeros(m,7e3),meanQz,zeros(m,7e3-1)]';
paddedAP = meanQz';
paddedAP(isnan(paddedAP(:))) = 0;
paddedAP = detrend(median(paddedAP,2),'constant');

n = length(paddedAP);
xdft = fft(paddedAP);
xdft = xdft(1:n/2+1);
psdx = abs(xdft/fs).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:fs/n:fs/2;
psd_Qz = interp1(freq,psdx,f,'linear','extrap');

psd = psd_Qx; psd(:,2) = psd_Qy; psd(:,3) = psd_Qz;
meanQ{1} = meanQx; meanQ{2} = meanQy; meanQ{3} = meanQz;



B2 = zeros(10,1);
X0 = mean(psd_passive(:,3,:),3);
figureNB(12,5);
for k = 1:6
    if(N(k)>0)
        % X1 = psd(:,k);
        X1 = psd(:,3);
        Y = psd_active(:,3,k);
        idcs = find(and(and(Y>1e-19,X1>1e-19),f<3e3));
        B2(k) = fminbnd(@(B1) neglnlike(B1,X0(idcs),X1(idcs),Y(idcs)),0,N(k)*2);
    else
        B2(k) = 0;
    end
end

m = ceil(max([firingFrequency;B2]));


blue = [17,82,185]/255;
red = clrsPT.qualitative_CM.red;

k=1;
figureNB(13.2,8.5);
for i = 1:3
    axes('Position',[0.27 0.83-(0.16+0.02)*(i-1) 0.24 0.15-0.02])
        plot(tSave,fullQ(:,i)-median(fullQ(:,i)),'color','k')
        hold on;
        plot(tSave,passiveQ(:,i)-median(passiveQ(:,i)),'color',[0.6,0.6,0.6],'LineWidth',1);
        axis off;
        xlim([0.2,2]*1e3);
        yl = get(gca,'ylim');
        ylim([yl(2)-270,yl(2)]);
    axes('Position',[0.56 0.83-0.16*(i-1) 0.11 0.13-0.02])
        plot(psd_passive(:,i,k),'color',[0.6,0.6,0.6],'LineWidth',1)
        hold on;
        plot(psd_active(:,i,k),'color','k','LineWidth',1)
        set(gca,'xscale','log')
        ylabel('log PSD')
        yticks([]);
        xlim([1,1e3])
        ylim([1e-8,1e1]);
        % yticks([1e-25,1e-20,1e-15])
        set(gca,'yscale','log')
        xticks([1,10,100,1000])
        if(i<3)
            xticklabels({});
            xlabel('');
        else
            % xticklabels([1,10,100,1000])
            xlabel('Frequency (Hz)')
        end
        % ylabel(['Power (' char(956) 'V^2/Hz)'])
        gcaformat;
        set(gca,'FontSize',8);
end

for i = 1:3
    axes('Position',[0.72 0.83-(0.16+0.02)*(i-1) 0.1 0.15-0.02])
        plot((-1000:1000)/16,meanQ{i}','color',blue*0.5+0.5)
        hold on;
        plot((-1000:1000)/16,nanmedian(meanQ{i}),'color',blue,'LineWidth',1);
        xlim([-10,15])
        ylabel(['Unitary AP (' char(956) 'V)'])
        xlabel('Time (ms)')
        axis off;
    axes('Position',[0.87 0.83-0.16*(i-1) 0.11 0.13-0.02])
        plot(f,psd(:,i),'color',blue,'LineWidth',1);
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        ylabel('log PSD')
        yticks([]);
        xlim([1,1e3])
        ylim([1e-6,1]);
        xticks([1,10,100,1000])
        if(i<3)
            xticklabels({});
            xlabel('');
        else
            % xticklabels([1,10,100,1000])
            xlabel('Frequency (Hz)')
        end
        % ylabel(['Power (' char(956) 'V^2/Hz)'])
        gcaformat;
        set(gca,'FontSize',8);
end

axes('Position',[0.06,0.16,0.12,0.2]);
    plot(1./EI,firingFrequency,'.-k','MarkerSize',10,'LineWidth',1)
    hold on;
    plot(1./EI,0*firingFrequency,'.-','MarkerSize',10,'LineWidth',1,'color',[0.6,0.6,0.6])
    set(gca,'xscale','log')
    xticks([1/30,1/10,1/3,1])
    xticklabels({'1:30','1:10','1:3','1:1'})
    xlabel('EI ratio (\lambda_E:\lambda_I)')
    ylabel('Firing rate (Hz)')
    gcaformat
    set(gca,'FontSize',8);
    xlim([1/10,1])
idcs = [5,4,1];
for k = 1:3
    axes('Position',[0.3+(k-1)*0.17,0.16,0.12,0.2]);
    plot(psd_passive(:,3,idcs(k)),'color',[0.6,0.6,0.6],'LineWidth',1)
    hold on;
    plot(psd_active(:,3,idcs(k)),'color','k','LineWidth',2)
    plot(B2(idcs(k))*X1,'color',blue,'LineWidth',1)
    set(gca,'xscale','log')
    xlim([1,1e3])
    yticks([1e-5,1e0])
    ylim([1e-8,1e1]);
    % ylim([1e-25,1e-14]);
    % yticks([1e-25,1e-20,1e-15])
    set(gca,'yscale','log')
    % title(sprintf('%.1f Hz',firingFrequency(k)))
    xticks([1,10,100,1000])
    % xticklabels([1,10,100,1000])
    xlabel('Frequency (Hz)')
    if(k==1)
        ylabel(['PSD (nA^2' char(956) 'm^2)'])
    end
    gcaformat;
    set(gca,'FontSize',8);
end

axes('Position',[0.86,0.16,0.12,0.2]);
    plot(firingFrequency,B2,'.k','MarkerSize',10,'LineWidth',1)
    line([0,10],[0,10],'color','k','LineWidth',1)
    xlabel('Firing rate (Hz)')
    % ylabel('\beta_1')
    ylabel('Scaling factor (\beta)')
    gcaformat(gca)
    set(gca,'FontSize',8);

gcaformat(gcf,true,8)