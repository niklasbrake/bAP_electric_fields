% [sa,X] = network_simulation_beluga.getHeadModel;
dt = 1/16e3;
tmax = 1;
N = tmax/dt+1;
t = 0:dt:tmax;

m1 = poissrnd(tmax*100);

t0 = tmax*rand(m1,1);

A = 0.2;
S = 1e-3;

m2 = floor(A*m1);
s1 = t0(randperm(m1,m2));%+S*randn(m2,1);
s2 = t0(randperm(m1,m2));%+S*randn(m2,1);

[common,I1,I2] = intersect(s1,s2);
s1 = s1+S*randn(m2,1);
s2 = s2+S*randn(m2,1);

figureNB(8,6);
axes('Position',[0.05,0.58,0.6,0.4]);
    R = raster([s1*0;s2*0+1],[s1;s2],gcf);
    R.Color = [0.6,0.6,0.6];
    hold on;
    R = raster(common*0,s1(I1),gcf);
    R.Color = 'k';
    R = raster(common*0+1,s2(I2),gcf);
    R.Color = 'k';
    xlim([0,tmax]);
    ylim([-0.2,2.1]);
    axis off;
    plot(common(2)+1e-3*[-30,30,30,-30,-30],[-0.2,-0.2,2.1,2.1,-0.2],':k','LineWidth',1)

axes('Position',[0.7,0.58,0.25,0.4]);
c = linspace(-5,5,1e3);
    fill([c,flip(c)],[0*c,0.9*exp(-c.^2./(2*(1e3*S)^2))],[0.8,0.8,0.8],'EdgeColor','none');
    hold on;
    % plot(c,exp(-c.^2./(2*8^2)),'color','r','LineWidth',1);
    line([0,0],[0,0.9],'color',0*[0.6,0.6,0.6],'LineWidth',1)
    del = 1e3*(s1(I1(2))-common(2));
    line([0,0]+del,[0,0.9],'color','k','LineWidth',1)

    fill([c,flip(c)],[1+0*c,1+0.9*exp(-c.^2./(2*(1e3*S)^2))],[0.8,0.8,0.8],'EdgeColor','none');
    hold on;
    % plot(c,1.2+exp(-c.^2./(2*8^2)),'color','r','LineWidth',1);
    line([0,0],1+[0,0.9],'color',0*[0.6,0.6,0.6],'LineWidth',1)
    del = 1e3*(s2(I2(2))-common(2));
    line([0,0]+del,1+[0,0.9],'color','k','LineWidth',1)

    plot([-5,5,5,-5,-5],[-0.2,-0.2,2.1,2.1,-0.2],':k','LineWidth',1)
    axis off;


load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\unitaryAP_all.mat')

X = zeros(N,2);
it = interp1(t,1:N,s1,'nearest','extrap');
[I,J] = findgroups(it);
X(J,1) = splitapply(@length,I,I);

it = interp1(t,1:N,s2,'nearest','extrap');
[I,J] = findgroups(it);
X(J,2) = splitapply(@length,I,I);

I = sample_blue_neurons(2);
% I = [535,312];


uAP = network_simulation_beluga.getEEG(savedUnitaryAP,sa,43e3);
y1 = uAP(:,I(1));
y2 = uAP(:,I(2));

Y(:,1) = filter(y1,1,X(:,1));
Y(:,2) = filter(y2,1,X(:,2));

t0 = 0.5*(s1(I1(2))+s2(I2(2)));


Z = Y+filter(exp(-(0:1/16:50)/2),1,5e-8*randn(size(Y)));
M1 = quantile(Z(:),0.001)
M2 = quantile(Z(:),0.999)
Z = (Z-M1)./(M2-M1);

axes('Position',[0.05,0.03,0.6,0.4]);
    hold on;
    plot(t-1/16,0.9*Z(:,1),'color','k');
    plot(t-1/16,1+0.9*Z(:,2),'color','r')
    F = fill(t0+[-30,30,30,-30,-30]*1e-3,[-0.2,-0.2,2.1,2.1,-0.2],[0.6,0.6,0.6],'FaceAlpha',0.2,'EdgeColor','k','LineStyle',':','LineWidth',1);
    F.FaceColor = 'none';
    axis off;
    ylim([-0.2,2.1]);
    xlim([0,1]);


axes('Position',[0.7,0.03,0.25,0.4]);
    F = fill(t0+[-5,5,5,-5,-5]*1e-3,[-0.2,-0.2,2.1,2.1,-0.2],[0.6,0.6,0.6],'FaceAlpha',0.2,'EdgeColor','k','LineStyle',':','LineWidth',1);
    F.FaceColor = 'none';
    hold on;
    plot(t-1/16,0.9*Z(:,1),'LineWidth',1,'color','k');
    plot(t-1/16,1+0.9*Z(:,2),'LineWidth',1,'color','r')

    % for i = 1:length(common)
    %     line([1,1]*s1(I1(i)),[60,45],'color','r','LineWidth',1)
    %     line([1,1]*s2(I2(i)),[75,60],'color','k','LineWidth',1)
    % end
    xlim(t0+[-5e-3,5e-3]);
    ylim([-0.2,2.1]);
    axis off;




% figureNB(8,3.5);
% axes('Position',[0.05,0.2,0.6,0.7]);
%     set(gca,'ColorOrder',[0,0,0;1,0,0]);
%     hold on;
%     plot(t-1/16,Y+filter(exp(-(0:1/16:50)/2),1,randn(size(Y))))
%     fill(t0+[-30,30,30,-30,-30]*1e-3,[-60,-60,85,85,-60],[0.6,0.6,0.6],'FaceAlpha',0.2,'EdgeColor','none');
%     axis off;
%     ylim([-60,100]);
%     xlim([0,1]);


% axes('Position',[0.7,0.2,0.25,0.7]);
%     fill(t0+[-5,5,5,-5,-5]*1e-3,[-60,-60,85,85,-60],[0.6,0.6,0.6],'FaceAlpha',0.2,'EdgeColor','none');
%     set(gca,'ColorOrder',[0,0,0;1,0,0]);
%     hold on;
%     plot(t-1/16,Y+filter(exp(-(0:1/16:50)/2),1,randn(size(Y))),'LineWidth',1)

%     for i = 1:length(common)
%         line([1,1]*s1(I1(i)),[60,45],'color','r','LineWidth',1)
%         line([1,1]*s2(I2(i)),[75,60],'color','k','LineWidth',1)
%     end
%     xlim(t0+[-5e-3,5e-3]);
%     axis off;
%     ylim([-60,100]);
