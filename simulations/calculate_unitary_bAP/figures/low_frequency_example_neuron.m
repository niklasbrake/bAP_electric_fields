data1=load('C:\Users\brake\Desktop\L6_UTPC_cADpyr231_5\matlab_recordings\synaptic_input_EI12.55.mat')
data2=load('C:\Users\brake\Desktop\L6_UTPC_cADpyr231_5\matlab_recordings\synaptic_input_EI12.55_passive.mat')
[y,x] = findpeaks(data1.voltage,data1.time,'MinPeakHeight',0);


% del = data1.dipoles(2:end,2,1)-data2.dipoles(2:end,2,1);
% del(del>=17.5) = 17.5;
% del(del<=-17.5) = -17.5;
% data1.dipoles(abs(data1.dipoles)>17.5) = nan;
% data2.dipoles(abs(data2.dipoles)>17.5) = nan;

red = clrsPT.qualitative_CM.red;
figureNB(13.2,10)
subplot(2,1,1);
    plot(data1.time(2:end),data1.dipoles(2:end,2,1),'color','k')
    hold on;
    plot(data1.time(2:end),data2.dipoles(2:end,2,1),'color',red)
    ylabel(['q_y (nA' char(956) 'm)'])
    % ylim([-20,20])

subplot(2,1,2);
    plot(data1.time(2:end),del,'k')
    line([0,1e4],[0,0],'color',red)
    ylim([-10,10])
    ylabel(['q_y active - passive (nA' char(956) 'm)'])


% for i = 1:length(x)
%     line(x(i)+[-50,50],0+[17,18],'color','k','LineWidth',1)
%     line(x(i)+[-50,50],1+[17,18],'color','k','LineWidth',1)

%     line(x(i)+[-50,50],0-[18,17],'color','k','LineWidth',1)
%     line(x(i)+[-50,50],-1-[18,17],'color','k','LineWidth',1)
% end

gcaformat(gcf,true,8)