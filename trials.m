FILE = 'eric.wav';
x = linspace(0,2*pi,2*pi*1000);
fc = 1;
y = cos(2*pi*fc*x);
yf = abs(fftshift(fft(y)));
f = linspace(-fc/2,fc/2,2*pi*1000);


plot(f,yf);