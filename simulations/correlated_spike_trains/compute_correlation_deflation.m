
load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\unitaryAP_all.mat');
[sa,X] = network_simulation_beluga.getHeadModel;
uAP = network_simulation_beluga.getEEG(savedUnitaryAP,sa,43e3);


r2 = @(i,j) xcorr(uAP(:,i),uAP(:,j))/(std(uAP(:,i))*std(uAP(:,j)))/size(uAP,1);

m = size(uAP,2);
t = m*(m-1)/2;
count = 0;
R2 = zeros(m,m);
lag = zeros(m,m);

R_x = zeros(m,m);
R_y = zeros(m,m);
R_z = zeros(m,m);

for i = 1:m
    for j = i+1:m
        count = count+1;
        waitbar(count/t);
        R_x(i,j) = corr(savedUnitaryAP(:,1,i),savedUnitaryAP(:,1,j));
        R_y(i,j) = corr(savedUnitaryAP(:,2,i),savedUnitaryAP(:,2,j));
        R_z(i,j) = corr(savedUnitaryAP(:,3,i),savedUnitaryAP(:,3,j));
        % [~,lag(i,j)] = max(abs(temp));
        % R2(i,j) = temp(lag(i,j));
    end
end

lag = lag+eye(m)+lag';
R2 = R2+eye(m)+R2';


R_x = R_x + R_x' + eye(m);
R_y = R_y + R_y' + eye(m);
R_z = R_z + R_z' + eye(m);

figureNB;
subplot(1,3,1);
    imagesc(R_x);
    set(gca,'CLim',[-1,1]);
    axis square
subplot(1,3,2);
    imagesc(R_y);
    set(gca,'CLim',[-1,1]);
    axis square
subplot(1,3,3);
    imagesc(R_z);
    set(gca,'CLim',[-1,1]);
    axis square
    colorbar

colormap(clrsPT.diverging(1e3));


for i = 1:1e3
    waitbar(i/1e3)
    bs_sample = interp1(F(idcs),idcs,R(:,i),'next','extrap');
    bootstrap_Est(:,i) = nanmean(psd_unit(:,bs_sample),2);
end

for i = 1:1e3
    waitbar(i/1e3)
    bs_sample = interp1(F(idcs),idcs,R(:,i),'next','extrap');
    bootstrap_Est(:,i) = nanmean(psd_unit(:,bs_sample),2);
end

