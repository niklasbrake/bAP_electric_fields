load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\unitaryAPNew.mat')
% load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\unitaryAP_all.mat');

n = 2000;
lag = 1e3;
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

R = R;
for count = 1:n
    R2(:,:,count,:) = R(:,:,count,:)./sqrt(diag(R11(:,:,count)).*diag(R22(:,:,count))');
end
R2 = permute(R2,[4,3,1,2]);

% Compute lead field in normal vector cooordinates
[sa,X] = network_simulation_beluga.getHeadModel;
M = size(X.vertices,1);
Lxyz = zeros(M,3);
for idx = 1:M
    L0 = squeeze(sa.cortex75K.V_fem(49,idx,:))'; % czIDX = 49
    vz = sa.cortex75K.normals(idx,:);
    [vx,vy] = getOrthBasis(vz);
    A = [vx(:),vy(:),vz(:)];
    Lxyz(idx,:) = 1e-6*L0*A;
end
E = permute(Lxyz,[2,3,1]).*permute(Lxyz,[3,2,1]);

E = E(:,:,sa.cortex10K.in_from_cortex75K);
M = size(E,3);

h = waitbar(0);
R_eeg = zeros(T,n);
for i = 1:M
    update_waitbar(h,i,M);
    Exy = squeeze(sum(sum(E(:,:,i).*R)));
    Exx = squeeze(sum(sum(E(:,:,i).*R11)));
    Eyy = squeeze(sum(sum(E(:,:,i).*R22)));
    R_eeg = R_eeg + 1/M*(Exy./sqrt(Exx.*Eyy))';
end


figureNB(3,3);
    plotwitherror((-lag:lag)/16,R_eeg,'CI','LineWidth',1,'color','k');
    xlabel('Lag (ms)');
    ylabel('Correlation')
    xlim([-10,10]);
    ylim([-2,6]);
    title('uapEEG')

R = permute(R,[4,3,1,2]);


figureNB(9,9);
subplot(3,3,1);
    plotwitherror((-lag:lag)/16,R(:,:,1,1),'CI','LineWidth',1,'color','k');
    xlabel('Lag (ms)');
    ylabel('Correlation')
    axis off;
    xlim([-10,10]);
    ylim([-2,6]);
    % title('Q_{xx}')
subplot(3,3,2);
    plotwitherror((-lag:lag)/16,R(:,:,1,2),'CI','LineWidth',1,'color','k');
    xlabel('Lag (ms)');
    ylabel('Correlation')
    axis off;
    xlim([-10,10]);
    ylim([-2,6]);
    % title('Q_{xy}')
subplot(3,3,3);
    plotwitherror((-lag:lag)/16,R(:,:,1,3),'CI','LineWidth',1,'color','k');
    xlabel('Lag (ms)');
    ylabel('Correlation')
    axis off;
    xlim([-10,10]);
    ylim([-2,6]);
    % title('Q_{xz}')
subplot(3,3,4);
    plotwitherror((-lag:lag)/16,R(:,:,2,1),'CI','LineWidth',1,'color','k');
    xlabel('Lag (ms)');
    ylabel('Correlation')
    axis off;
    xlim([-10,10]);
    ylim([-2,6]);
    % title('Q_{yx}')
subplot(3,3,5);
    plotwitherror((-lag:lag)/16,R(:,:,2,2),'CI','LineWidth',1,'color','k');
    xlabel('Lag (ms)');
    ylabel('Correlation')
    axis off;
    xlim([-10,10]);
    ylim([-2,6]);
    % title('Q_{yy}')
subplot(3,3,6);
    plotwitherror((-lag:lag)/16,R(:,:,2,3),'CI','LineWidth',1,'color','k');
    xlabel('Lag (ms)');
    ylabel('Correlation')
    axis off;
    xlim([-10,10]);
    ylim([-2,6]);
    % title('Q_{yz}')
subplot(3,3,7);
    plotwitherror((-lag:lag)/16,R(:,:,3,1),'CI','LineWidth',1,'color','k');
    xlabel('Lag (ms)');
    ylabel('Correlation')
    axis off;
    xlim([-10,10]);
    ylim([-2,6]);
    % title('Q_{zx}')
subplot(3,3,8);
    plotwitherror((-lag:lag)/16,R(:,:,3,2),'CI','LineWidth',1,'color','k');
    xlabel('Lag (ms)');
    ylabel('Correlation')
    axis off;
    xlim([-10,10]);
    ylim([-2,6]);
    % title('Q_{zy}')
subplot(3,3,9);
    plotwitherror((-lag:lag)/16,R(:,:,3,3),'CI','LineWidth',1,'color','k');
    xlabel('Lag (ms)');
    ylabel('Correlation')
    axis off;
    xlim([-10,10]);
    ylim([-2,6]);
    % title('Q_{zz}')
gcaformat(gcf);
