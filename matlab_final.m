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
%plot_in_time(yt,fs,'Signal in time domain');

%-------------------------------------------------message in frequency domain--------------------------------------------------------------------
yf = fftshift(fft(yt));
%-------------------------------------------------plotting--------------------------------------------------------------------
%plot_in_frequency(real(yf),fs,'Signal in frequency domain');

%-------------------------------------------------filter--------------------------------------------------------------------
filter = generate_filter(length(yf),fs,f_filter);
yf_filtered = filter.*yf;
%-------------------------------------------------plotting--------------------------------------------------------------------
%plot_in_frequency(real(filter),fs,'Filter');
%plot_in_frequency(real(yf_filtered),fs,'Filtered message in frequency domain');

%-------------------------------------------------back to time domain--------------------------------------------------------------------
yt_filtered= real(ifft(ifftshift(yf_filtered)));
%-------------------------------------------------sound--------------------------------------------------------------------
% sound(yt_filtered,fs);
%-------------------------------------------------plotting--------------------------------------------------------------------
t = linspace(0,length(yt_filtered)/fs, length(yt_filtered));
%plot_in_time(yt_filtered,fs, 'Filtered message in time domain');
%-------------------------------------------------resmple--------------------------------------------------------------------
yt_resampled = resample(yt_filtered,fs_res,fs);

%-------------------------------------------------carrier signal--------------------------------------------------------------------
t_carr = linspace(0,length(yt_resampled)/fs_res, length(yt_resampled)); %(x2-x1)/(n-1) = 1/5*fc, linspace(x1,x2,n)
carrier_t = cos(2*pi*fc*t_carr).';
carrier_f = fftshift(fft(carrier_t));




%-------------------------------------------------DSBSC modulation--------------------------------------------------------------------
yt_sc = carrier_t.*yt_resampled;
yf_sc = fftshift(fft(yt_sc));
%sound(yt_sc, fs);
%-------------------------------------------------plotting--------------------------------------------------------------------
%plot_in_time(yt_sc,fs_res,'DSBSC modulated signal in time domain');
%plot_in_frequency(real(yf_sc),fs_res,'DSBSC modulated signal in frequency domain');




%-------------------------------------------------DSBTC modulation--------------------------------------------------------------------
max_message = max(abs(yt_resampled));
A = 2*max_message;
%message_normalized = yt_resampled./A;
%dada = (A+yt_resampled)*carrier_t;
%yt_tc = A.*(1+modulation_index.*message_normalized).*carrier_t;
yt_tc = (A+yt_resampled).*carrier_t;
yf_tc = fftshift(fft(yt_tc));
%-------------------------------------------------plotting--------------------------------------------------------------------
%plot_in_time(yt_tc,fs_res,'DSBTC modulated signal in time domain');
%plot_in_frequency(real(yf_tc),fs_res,'DSBTC modulated signal in frequency domain');
%sound(yt_tc,fs);


%-------------------------------------------------DSBSC envelope detector demodulation--------------------------------------------------------------------
[yt_sc_env, yf_sc_env] = env_demod(yt_sc,fs_res,fs,0,0);
%-------------------------------------------------plotting--------------------------------------------------------------------
%plot_in_time(yt_sc_env, fs, 'DSBSC demodulated signal in time domain using envelope detector');
%plot_in_frequency(real(yf_sc_env), fs, 'DSBSC demodulated signal in frequency domain using envelope detector');
%-------------------------------------------------sound--------------------------------------------------------------------
%sound(yt_sc_env,fs); %should not be good



%-------------------------------------------------DSBTC envelope detector demodulation--------------------------------------------------------------------
[yt_tc_env, yf_tc_env] = env_demod(yt_tc,fs_res,fs,0,0);
%-------------------------------------------------plotting--------------------------------------------------------------------
%plot_in_time(yt_tc_env, fs, 'DSBTC demodulated signal in time domain using envelope detector');
%plot_in_frequency(real(yf_tc_env), fs, 'DSBTC demodulated signal in frequency domain using envelope detector');
%-------------------------------------------------sound--------------------------------------------------------------------
%sound(yt_tc_env,fs); %should be good



%-------------------------------------------------DSBSC coherent demodulation--------------------------------------------------------------------

[yt_sc_coh, yf_sc_coh] = coh_demod(yt_sc,fs_res,fs,0,0,fc,0,f_filter);
%-------------------------------------------------plotting--------------------------------------------------------------------
plot_in_time(real(yt_sc_coh), fs, 'DSBSC demodulated signal in time domain using coherent detector');
plot_in_frequency(real(yf_sc_coh), fs, 'DSBSC demodulated signal in frequency domain using coherent detector');
%-------------------------------------------------sound--------------------------------------------------------------------
sound(real(yt_sc_coh),fs); %should be good



























%-------------------------------------------------DSBSC coherent demodulation--------------------------------------------------------------------

















%-------------------------------------------------DSBTC envelope detector demodulation--------------------------------------------------------------------
% 
% yt_dsbtc_demod_res = resample(abs(hilbert(yt_tc)),fs,fs_res);
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
% yt_tc = awgn(yt_tc, 50);
% %-------------------------------------------------DSBTC envelope detector demodulation--------------------------------------------------------------------
% yt_dsbtc_demod = abs(hilbert(yt_tc));
% yt_dsbtc_demod_res = resample(yt_dsbtc_demod ,fs,fs_res);
% yf_dsbtc_demod_res = fftshift(fft(yt_dsbtc_demod_res));
% sound(yt_dsbtc_demod_res, fs);
% plot_signal(xf_dsbsc_demod_res, yf_dsbtc_demod_res, 'demodulated signal in frequency domain', 'Frequency', 'Value');
% %---------------------------------------------------------------------------------------------------------------
% receivedSignal = awgn(signal, 30);
% 
