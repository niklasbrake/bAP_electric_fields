[sa,X] = network_simulation_beluga.getHeadModel;

folder = 'E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\other_code\representative_models';
F = dir(folder);
F = F(3:end);

for i = 1:length(F)
    items = split(F(i).name,'_');
    layer{i} = items{1};
    morph{i} = items{2};
    if(strcmp(items{3},'L1') || strcmp(items{3},'L4'))
        morph{i} = [items{2} '_' items{3}];
    end
    mtype{i} = [layer{i} '_' morph{i}];
end
unique(layer)
unique(morph)
e_mtypes = {'PC','SS','SP','STPC','UTPC','TTPC1','TTPC2','TPC','BPC','IPC'};
% i_mtypes = setdiff(unique(morph),e_mtypes);
ei_type = ismember(morph,e_mtypes);

EI_vec = {'01','1.5','2.1','3.1','4.5','6.6','9.7','14.1','20.6','30'};

h = waitbar(0);
for i = 1:100
    waitbar(i/100,h)
    folder0 = fullfile(folder,F(i).name,'matlab_recordings');

    N = zeros(length(EI_vec),1);
    psd2 = zeros(16385,length(EI_vec));
    count2 = 1;
    saveAP = [];
    for k = 1:length(EI_vec)
        load(fullfile(folder0,sprintf('synaptic_input_EI%s.mat',EI_vec{k})))
        eeg = network_simulation_beluga.getEEG(dipoles,sa,1e3);
        Y = detrend(eeg(16001:end));
        [y,x] = findpeaks(voltage,'MinPeakHeight',0);
        [psd2(:,k),f0] = pmtm(Y,2,[],16e3);
        N(k) = sum(x>16e3);
        [~,I] = min(N(N>0));
        if(and(N(k)>0,N(k)<30))
            unitaryAP = zeros(1,2001);
            count = 1;
            for j = 1:length(x)
                if(x(j)<16e3)
                    continue;
                end
                idcs = max(min(x(j)-1e3:x(j)+1e3,length(eeg)),1);
                % y = dipoles(idcs,3);
                y = eeg(idcs(:));
                y(idcs==1,:) = 0;
                y(idcs==length(eeg),:) = 0;
                if(j>1)
                    idcs0 = find(idcs<=x(j-1)+100);
                    y(idcs0,:) = nan;
                end
                if(j<length(x))
                    idcs0 = find(idcs>=x(j+1)-100);
                    y(idcs0,:) = nan;
                end
                unitaryAP(count,:) = y(:)';
                count = count + 1;
            end
            saveAP(count2,:) = nanmedian(unitaryAP,1);
            count2 = count2+1;
        end
    end
    if(count2==1)
        savedUnitaryAP(:,i) = nan;
    else
        savedUnitaryAP(:,i) = sum(saveAP.*N(and(N>0,N<30)),1)/sum(N(and(N>0,N<30)));
    end
end
delete(h);

idcs = ~isnan(savedUnitaryAP(1,:));

figureNB(3.2,3.2)
    histogram(log10(pEst2),'BinWidth',0.2,'EdgeColor','none','FaceColor',[0.1,0.1,0.1])
    xlim([-18,-12])
    xticks([-18:2:-12])
    xticklabels({'10^{-18}','10^{-16}','10^{-14}','10^{-12}'});
    xlabel({'Unitary AP EEG response',['average power (' char(956) 'V^2)']})
    yticks([]);




Y = savedUnitaryAP./max(abs(savedUnitaryAP));
% Y = Y(:,idcs);
Y = Y-median(Y);
Y = Y/max(Y(:));
[coeff,score,~,~,explained] = pca(Y');


% S.mu = [-3,-2;10,2];
% S.Sigma = cat(3,[0.2,0;0,0.2],[1,-1;-1,5]);
% GMModel = fitgmdist(score(:,1:2),2,'Start',S);
GMModel = fitgmdist(score(:,1:2),2);
gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(GMModel,[x0 y0]),x,y);
J = GMModel.cluster(score(:,1:2));

figureNB(9,9);
subplot(2,1,1)
    % gscatter(score(:,1),score(:,2),layer(idcs)');
    plot(score(:,1),score(:,2),'.k','MarkerSize',10);
    hold on;
    g = gca;
    fcontour(gmPDF,[[-4,4],[-4,4]],'LineWidth',1,'MeshDensity',100)
    % fcontour(gmPDF,[[-4,4],[-4,4]],'LineWidth',1,'MeshDensity',100,'LevelList',[0.05,0.1,0.3,0.4,1])
    % plot(score(:,1),score(:,2),'.k','MarkerSize',10);
    gscatter(score(:,1),score(:,2),J)
    legend off
    xlabel('PC 1')
    ylabel('PC 2')
    set(gca,'DataAspectRatio',100-explained(1:3))
    % xlim([-2,2])
    % ylim([-1,1])
    % set(gca,'DataAspectRatio',[1,1,1])

subplot(2,2,3)
    plot((-1000:1000)/16,mean(Y(:,J==1),2),'LineWidth',1.5);
    hold on;
    plot((-1000:1000)/16,mean(Y(:,J==2),2),'LineWidth',1.5);
    % plot((-1000:1000)/16,mean(Y(:,J==3),2),'LineWidth',1.5);
    xlim([-5,5]);
    ylabel(['Unitary AP (norm.)'])
    xlabel('Time (ms)')


y1 = mean(Y(:,J==1),2);
y1 = y1-median(y1);
[psd1,f] = pmtm(y1,2,[],16e3);

y2 = mean(Y(:,J==2),2);
y2 = y2-median(y2);
[psd2,f] = pmtm(y2,2,[],16e3);

% y3 = mean(Y(:,J==3),2);
% y3 = y3-median(y3);
% [psd3,f] = pmtm(y3,2,[],16e3);



subplot(2,2,4)
    plot(f,psd1,'LineWidth',2);
    hold on;
    plot(f,psd2,'LineWidth',2)
    % plot(f,psd3,'LineWidth',2)
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlabel('Frequency (Hz)')
    ylabel('Power (a.u.)')
    yticks([])

gcaformat(gcf);

Y = savedUnitaryAP(:,idcs)-median(savedUnitaryAP(:,idcs));
Y = [zeros(1e4,99);Y;zeros(1e4,99)];

[psd,f] = pmtm(Y,2,[],16e3);
figureNB;
subplot(1,2,1);
    plotwitherror((-1000:1000)/16,Y(1e4:1e4+2000,:),'Q','color','k','LineWidth',1);
    xlabel('Time (ms)')
    ylabel(['Unitary AP (' char(956) 'V)'])
    ylim([-2e-6,2e-6])
    xlim([-2,10])
    xticks([0:5:10])
    gcaformat
subplot(1,2,2);
    plotwitherror(f(2:end),psd(2:end,:),'Q','color','k','LineWidth',1);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlabel('Frequency (Hz)')
    ylabel(['Power (' char(956) 'V^2/Hz)'])
    gcaformat




[I,J] = findgroups(morph);
grouped = splitapply(@(x) mean(x,2),savedUnitaryAP,I(:)');
figureNB;
for i=  1:length(J)
    subplot(5,5,i)
    % plot((-1000:1000)/16,grouped(:,i),'color','k','LineWidth',1.5)

    plot((-1000:1000)/16,savedUnitaryAP(:,I==i));
    title(J{i})
    xlim([-10,10])
    ylim([-5,5]*1e-6)
    axis off;
end