%-------------------------------------------------clear--------------------------------------------------------------------
clc;
close all;
clear all;
%-------------------------------------------------read file--------------------------------------------------------------------
FILE = 'eric.wav';
fc = 100000;
[yt, fs]= audioread(FILE);
fs_res = 5*fc;
f_filter = 4000;
%-------------------------------------------------plotting--------------------------------------------------------------------
xt = linspace(0,length(yt)/fs,length(yt));
%plot_in_time(yt,fs,'Message in time domain');
plot_in_frequency(real(fftshift(fft(yt))),fs, 'message in frequency');
%-------------------------------------------------message in frequency domain--------------------------------------------------------------------
yf = fftshift(fft(yt));
xf = linspace(-fs/2,fs/2,length(yf));
%-------------------------------------------------filter--------------------------------------------------------------------
filter = generate_filter(length(yf),fs,f_filter);
yf_filtered = filter.*yf;
%-------------------------------------------------back to time domain--------------------------------------------------------------------
yt_filtered= real(ifft(ifftshift(yf_filtered)));
plot_in_time(yt_filtered, fs, 'Signal in time domain after filter' )
%-------------------------------------------------sound--------------------------------------------------------------------
% sound(yt_filtered,fs);
%-------------------------------------------------plotting--------------------------------------------------------------------
xt = linspace(0,length(yt_filtered)/fs, length(yt_filtered));
%-------------------------------------------------resmple--------------------------------------------------------------------
yt_resampled = resample(yt_filtered,fs_res,fs);
%-------------------------------------------------carrier signal--------------------------------------------------------------------
t = linspace(0,length(yt_resampled)/fs_res, length(yt_resampled)); %(x2-x1)/(n-1) = 1/5*fc, linspace(x1,x2,n) lama radawan yeegy neb2a nshofha
carrier = cos(2*pi*fc*t).';
carriersin = sin(2*pi*fc*t).';
yt_carrier=carrier;
yf_carrier=fftshift(fft(yt_carrier));
xt_carrier=t;
xf_carrier=linspace(-fc/2,fc/2,length(yt_carrier));
%------------------------------------------------FM Modulation------------------------------------------------------------------------
%------------------------------------------------NBFM------------------------------------------------------------------------
A = max(abs(yt));
kf = pi; %de ely btefre2 fel amplitude bta3 elcarrier
beta = (kf*A)/(2*pi*fs_res);
m_int = kf.*cumsum(yt_resampled).'; % Integrating Msg
St = A.*cos(2*pi*fc*t + m_int);
%plot_in_time(St, fs_res, 'modulated in time');
plot_in_frequency(abs(fftshift(fft(St))),fs_res,'modulated signal in frequency');
%St = A.*cos(2*pi*fc*t)-m_int;%.*sin(2*pi*fc*t);
fourier = fftshift(fft(St));
%------------------------------------------------Demodulation------------------------------------------------------------------------
%------------------------------convert FM to AM----------------------------
%AM = A.*(2*pi*fc+kf*yt_resampled).*sin(2*pi*fc*t+m_int);
AM = diff(St);
AM = [0 AM];
%plot_in_frequency(abs(fftshift(fft(AM))), fs, 'lel');
[yt_demod, yf_demod] = env_demod(AM,fs_res,fs,0,0);

%plot_in_time(yt_demod,fs,'Envelope Detector');
%plot_in_frequency(yf_demod,fs,'demodulated signal in frequency');

sound(yt_demod,fs);