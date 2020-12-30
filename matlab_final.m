%-------------------------------------------------clear--------------------------------------------------------------------
clc;
close all;
clear all;
%-------------------------------------------------read file--------------------------------------------------------------------
FILE = 'eric.wav';
fc = 100000;
fs_res = 5*fc;
modulation_index = 0.5;
[yt, fs]= audioread(FILE);
f_filter = 4000;
%-------------------------------------------------plotting--------------------------------------------------------------------
t = linspace(0,length(yt)/fs,length(yt));
% plot_signal(t,yt,'Signal in time domain','time','Value');
%-------------------------------------------------message in frequency domain--------------------------------------------------------------------
yf = fftshift(fft(yt));
f = linspace(-fs/2,fs/2,length(yf));
%-------------------------------------------------plotting--------------------------------------------------------------------
% plot_signal(f,real(yf),'Signal in frequency domain','Frequency','Value');
%-------------------------------------------------filter--------------------------------------------------------------------
filter = generate_filter(length(yf),f,f_filter);
yf_filtered = filter.*yf;
%-------------------------------------------------plotting--------------------------------------------------------------------
% plot_signal(f,real(filter),'Filter','Frequency','Value');
% plot_signal(f,real(yf_filtered),'Filtered message in frequency domain','Frequency', 'Value');
%-------------------------------------------------back to time domain--------------------------------------------------------------------
yt_filtered= real(ifft(ifftshift(yf_filtered)));
%-------------------------------------------------sound--------------------------------------------------------------------
% sound(yt_filtered,fs);
%-------------------------------------------------plotting--------------------------------------------------------------------
t = linspace(0,length(yt_filtered)/fs, length(yt_filtered));
% plot_signal(t, yt_filtered, 'Filtered message in time domain','Time','Value');
%-------------------------------------------------resmple--------------------------------------------------------------------
yt_resampled = resample(yt_filtered,fs_res,fs);
f_res = linspace(-fs_res/2,fs_res/2,length(yt_resampled));
%-------------------------------------------------carrier signal--------------------------------------------------------------------
t_res = linspace(0,length(yt_resampled)/fs_res, length(yt_resampled)); %(x2-x1)/(n-1) = 1/5*fc, linspace(x1,x2,n)
carrier = cos(2*pi*fc*t_res).';
%-------------------------------------------------DSBSC modulation--------------------------------------------------------------------
yt_dsbsc = carrier.*yt_resampled;
yf_dsbsc = fftshift(fft(yt_dsbsc));
%sound(yt_dsbsc, fs);
%-------------------------------------------------plotting--------------------------------------------------------------------
%plot_signal(t_res,yt_dsbsc,'DSBSC modulated signal in time domain', 'Time', 'Value');
%plot_signal(f_res,real(yf_dsbsc),'DSBSC modulated signal in frequency domain','Frequency','Value');
%-------------------------------------------------DSBTC modulation--------------------------------------------------------------------
A = max(abs(yt_resampled));
message_normalized = yt_resampled./A;
yt_dsbtc = A.*(1+modulation_index.*message_normalized).*carrier;
yf_dsbtc = fftshift(fft(yt_dsbtc));
%-------------------------------------------------plotting--------------------------------------------------------------------
plot_signal(t_res,yt_dsbtc,'DSBTC modulated signal in time domain', 'Time', 'Value');
plot_signal(f_res,real(yf_dsbtc),'DSBTC modulated signal in frequency domain', 'Frequency', 'Value');
%-------------------------------------------------DSBSC envelope detector demodulation--------------------------------------------------------------------
% yt_dsbsc_demod_ = abs(hilbert(yt_dsbsc));
% yt_dsbsc_demod_res_ = resample(yt_dsbsc_demod_ ,fs,fs_res);
% yf_dsbsc_demod_res_ = fftshift(fft(yt_dsbsc_dtc,real(yf_dsbtc),'DSBTC modulated signal in frequency domain','Frequency','Value'));
%-------------------------------------------------DSBSC coherent demodulation--------------------------------------------------------------------
% filter = generate_filter(length(yt_dsbsc),t_res,f_filter);
% yt_dsbsc_demod = (yt_dsbsc.*carrier).*filter;
% plot_signal(t_res,yt_dsbsc_demod,'DSBSC demodulated signal in time domain', 'Time', 'Value');
% yt_dsbsc_demod_res = resample(yt_dsbsc_demod ,fs,fs_res);
% yf_dsbsc_demod_res = fftshift(fft(yt_dsbsc_demod_res));
% xf_dsbsc_demod_res = linspace(-fs/2,fs/2,length(yf_dsbsc_demod_res));
% plot_signal(xf_dsbsc_demod_res, yf_dsbsc_demod_res, 'demodulated signal in frequency domain', 'Frequency', 'Value');
% sound(yt_dsbsc_demod_res, fs);
% sound(yt_dsbsc_demod_res_, fs);
% plot_signal(xf_dsbsc_demod_res, yf_dsbsc_demod_res_, 'DSBSC demodulated signal using envelope detector', 'Frequency', 'Value');
% %-------------------------------------------------DSBTC envelope detector demodulation--------------------------------------------------------------------
% yt_dsbtc_demod = abs(hilbert(yt_dsbtc));
% yt_dsbtc_demod_res = resample(yt_dsbtc_demod ,fs,fs_res);
% yf_dsbtc_demod_res = fftshift(fft(yt_dsbtc_demod_res));
% sound(yt_dsbtc_demod_res, fs);
% plot_signal(xf_dsbsc_demod_res, yf_dsbtc_demod_res, 'demodulated signal in frequency domain', 'Frequency', 'Value');
% %-------------------------------------------------modulation test-----------------------------------------------------------------
% [yt_modulated, t_res] = modulate(yt_resampled, fc, fs_res);
% yf_modulated = fftshift(fft(yt_modulated));
% f = linspace(-fc/2,fc/2,length(yf_modulated));
% plot_signal(t_res,yt_modulated,'','','');
% plot_signal(f,yf_modulated,'','','');
% 
% yt_dsbtc = awgn(yt_dsbtc, 50);
% %-------------------------------------------------DSBTC envelope detector demodulation--------------------------------------------------------------------
% yt_dsbtc_demod = abs(hilbert(yt_dsbtc));
% yt_dsbtc_demod_res = resample(yt_dsbtc_demod ,fs,fs_res);
% yf_dsbtc_demod_res = fftshift(fft(yt_dsbtc_demod_res));
% sound(yt_dsbtc_demod_res, fs);
% plot_signal(xf_dsbsc_demod_res, yf_dsbtc_demod_res, 'demodulated signal in frequency domain', 'Frequency', 'Value');
% %---------------------------------------------------------------------------------------------------------------
% receivedSignal = awgn(signal, 30);

