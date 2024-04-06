y1 = uAP(:,1);
y2 = uAP(:,2);

fs = 16e3;
N = length(y1);
f = fs/N:fs/N:fs/2;

SN = fft(xcorr(y1+y2,y1+y2,1e3))/(fs*N);
S1 = fft(xcorr(y1,y1,1e3))/(fs*N);
S2 = fft(xcorr(y2,y2,1e3))/(fs*N);
S12 = fft(xcorr(y1,y2,1e3))/(fs*N);
SN = 2*abs(SN(2:1001));
S1 = 2*abs(S1(2:1001));
S2 = 2*abs(S2(2:1001));
S12 = 2*S12(2:1001);

S12 = real(S12.*exp(2*pi*sqrt(-1)*f(:)*1001/fs));
% S12 = 0*real(S12);
% S12 = 0.5*(S1+S2);

figureNB;
    plot(f,SN,'b')
    hold on;
    plot(f,S1+S2+2*S12,'c')
    set(gca,'xscale','log')
% set(gca,'yscale','log')



y1 = uAP(:,1);
y2 = [uAP(1700:2001,2);uAP(1:1699,2)];
y2 = [uAP(1960:2001,2);uAP(1:1959,2)];
y4 = [y2(1960:2001);y2(1:1959)];


fs = 16e3;
N = length(y1);
f = fs/N:fs/N:fs/2;

SN = fft(xcorr(y1+y2,y1+y2,1e3))/(fs*N);
S1 = fft(xcorr(y1,y1,1e3))/(fs*N);
S2 = fft(xcorr(y2,y2,1e3))/(fs*N);
S12 = fft(xcorr(y1,y2,1e3))/(fs*N);
SN = 2*abs(SN(2:1001));
S1 = 2*abs(S1(2:1001));
S2 = 2*abs(S2(2:1001));
S12 = 2*S12(2:1001);

S12 = real(S12.*exp(2*pi*sqrt(-1)*f(:)*1001/fs));

% figureNB;
    plot(f,SN,'r')
    hold on;
    plot(f,S1+S2+2*S12,'g')
    set(gca,'xscale','log')
% set(gca,'yscale','log')