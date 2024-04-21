load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\unitaryAPNew.mat')
% load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\unitaryAP_all.mat');

n = 2000;
lag = 2e3;
T = 2*lag+1;

R_xx = zeros(T,n); R_yx = zeros(T,n); R_zx = zeros(T,n);
R_xy = zeros(T,n); R_yy = zeros(T,n); R_zy = zeros(T,n);
R_xz = zeros(T,n); R_yz = zeros(T,n); R_zz = zeros(T,n);

R11 = zeros(3,3,n);
R22 = zeros(3,3,n);

bs_sample = reshape(sample_blue_neurons(2*n),2,[]);
h = waitbar(0);
for count = 1:n
    update_waitbar(h,count,n);
    i = bs_sample(1,count);
    j = bs_sample(2,count);
    R_xx(:,count) = xcorr(savedUnitaryAP(:,1,i),savedUnitaryAP(:,1,j),lag,'unbiased');
    R_xy(:,count) = xcorr(savedUnitaryAP(:,1,i),savedUnitaryAP(:,2,j),lag,'unbiased');
    R_xz(:,count) = xcorr(savedUnitaryAP(:,1,i),savedUnitaryAP(:,3,j),lag,'unbiased');

    R_yx(:,count) = xcorr(savedUnitaryAP(:,2,i),savedUnitaryAP(:,1,j),lag,'unbiased');
    R_yy(:,count) = xcorr(savedUnitaryAP(:,2,i),savedUnitaryAP(:,2,j),lag,'unbiased');
    R_yz(:,count) = xcorr(savedUnitaryAP(:,2,i),savedUnitaryAP(:,3,j),lag,'unbiased');

    R_zx(:,count) = xcorr(savedUnitaryAP(:,3,i),savedUnitaryAP(:,1,j),lag,'unbiased');
    R_zy(:,count) = xcorr(savedUnitaryAP(:,3,i),savedUnitaryAP(:,2,j),lag,'unbiased');
    R_zz(:,count) = xcorr(savedUnitaryAP(:,3,i),savedUnitaryAP(:,3,j),lag,'unbiased');

    R11(:,:,count) = cov(savedUnitaryAP(:,:,i))*2e3;
    R22(:,:,count) = cov(savedUnitaryAP(:,:,j))*2e3;
end


R = zeros(3,3,n,T);
R(1,1,:,:) = R_xx';
R(1,2,:,:) = R_xy';
R(1,3,:,:) = R_xz';
R(2,1,:,:) = R_yx';
R(2,2,:,:) = R_yy';
R(2,3,:,:) = R_yz';
R(3,1,:,:) = R_zx';
R(3,2,:,:) = R_zy';
R(3,3,:,:) = R_zz';
R = permute(R,[4,3,1,2]);
t = (-lag:lag)/16;

figureNB(9,9);
axes('Position',[0,0,1/3,1/3]);
    plotwitherror(t,1./(abs(t)<10).*R(:,:,1,1),'CI','LineWidth',1,'color','k');
    xlabel('Lag (ms)');
    ylabel('Correlation')
    axis off;
    xlim([-12,12]);
    ylim([-2,6]);
    % title('Q_{xx}')
axes('Position',[0,1/3,1/3,1/3]);
    plotwitherror(t,1./(abs(t)<10).*R(:,:,1,2),'CI','LineWidth',1,'color','k');
    xlabel('Lag (ms)');
    ylabel('Correlation')
    axis off;
    xlim([-12,12]);
    ylim([-2,6]);
    % title('Q_{xy}')
axes('Position',[0,2/3,1/3,1/3]);
    plotwitherror(t,1./(abs(t)<10).*R(:,:,1,3),'CI','LineWidth',1,'color','k');
    xlabel('Lag (ms)');
    ylabel('Correlation')
    axis off;
    xlim([-12,12]);
    ylim([-2,6]);
    % title('Q_{xz}')
axes('Position',[1/3,0,1/3,1/3]);
    plotwitherror(t,1./(abs(t)<10).*R(:,:,2,1),'CI','LineWidth',1,'color','k');
    xlabel('Lag (ms)');
    ylabel('Correlation')
    axis off;
    xlim([-12,12]);
    ylim([-2,6]);
    % title('Q_{yx}')
axes('Position',[1/3,1/3,1/3,1/3]);
    plotwitherror(t,1./(abs(t)<10).*R(:,:,2,2),'CI','LineWidth',1,'color','k');
    xlabel('Lag (ms)');
    ylabel('Correlation')
    axis off;
    xlim([-12,12]);
    ylim([-2,6]);
    % title('Q_{yy}')
axes('Position',[1/3,2/3,1/3,1/3]);
    plotwitherror(t,1./(abs(t)<10).*R(:,:,2,3),'CI','LineWidth',1,'color','k');
    xlabel('Lag (ms)');
    ylabel('Correlation')
    axis off;
    xlim([-12,12]);
    ylim([-2,6]);
    % title('Q_{yz}')
axes('Position',[2/3,0,1/3,1/3]);
    plotwitherror(t,1./(abs(t)<10).*R(:,:,3,1),'CI','LineWidth',1,'color','k');
    xlabel('Lag (ms)');
    ylabel('Correlation')
    axis off;
    xlim([-12,12]);
    ylim([-2,6]);
    % title('Q_{zx}')
axes('Position',[2/3,1/3,1/3,1/3]);
    plotwitherror(t,1./(abs(t)<10).*R(:,:,3,2),'CI','LineWidth',1,'color','k');
    xlabel('Lag (ms)');
    ylabel('Correlation')
    axis off;
    xlim([-12,12]);
    ylim([-2,6]);
    % title('Q_{zy}')
axes('Position',[2/3,2/3,1/3,1/3]);
    plotwitherror(t,1./(abs(t)<10).*R(:,:,3,3),'CI','LineWidth',1,'color','k');
    xlabel('Lag (ms)');
    ylabel('Correlation')
    axis off;
    xlim([-12,12]);
    ylim([-2,6]);
    % title('Q_{zz}')
gcaformat(gcf);
