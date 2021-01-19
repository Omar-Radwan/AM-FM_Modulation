%{
    plot_in_time
    plot_in_frequency
    generate_filter
    env_demod
    coh_demod
%}


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
max_amp = max(abs(yt));
f_filter = 4000;
%-------------------------------------------------plotting--------------------------------------------------------------------
plot_in_time(yt,fs,'Signal in time domain');
%-------------------------------------------------message in frequency domain--------------------------------------------------------------------
yf = fftshift(fft(yt));
%-------------------------------------------------plotting--------------------------------------------------------------------
plot_in_frequency(real(yf),fs,'Signal in frequency domain');
%-------------------------------------------------filter--------------------------------------------------------------------
filter = generate_filter(length(yf),fs,f_filter);
yf_filtered = filter.*yf;
%-------------------------------------------------plotting--------------------------------------------------------------------
plot_in_frequency(real(filter),fs,'Filter');
plot_in_frequency(real(yf_filtered),fs,'Filtered message in frequency domain');

%-------------------------------------------------back to time domain--------------------------------------------------------------------
yt_filtered= real(ifft(ifftshift(yf_filtered)));
%-------------------------------------------------sound--------------------------------------------------------------------
%sound(yt_filtered,fs);
%-------------------------------------------------plotting--------------------------------------------------------------------
t = linspace(0,length(yt_filtered)/fs, length(yt_filtered));
plot_in_time(yt_filtered,fs, 'Filtered signal in time domain');
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
plot_in_time(yt_sc,fs_res,'DSBSC modulated signal in time domain');
plot_in_frequency(real(yf_sc),fs_res,'DSBSC modulated signal in frequency domain');

%-------------------------------------------------DSBTC modulation--------------------------------------------------------------------
max_message = max(abs(yt_resampled));
A = 2*max_message;
yt_tc = (A+yt_resampled).*carrier_t;
yf_tc = fftshift(fft(yt_tc));
%-------------------------------------------------plotting--------------------------------------------------------------------
plot_in_time(yt_tc,fs_res,'DSBTC modulated signal in time domain');
plot_in_frequency(real(yf_tc),fs_res,'DSBTC modulated signal in frequency domain');
%sound(yt_tc,fs);


%-------------------------------------------------DSBSC envelope detector demodulation--------------------------------------------------------------------
[yt_sc_env, yf_sc_env] = env_demod(yt_sc,fs_res,fs,0,0);
%-------------------------------------------------plotting--------------------------------------------------------------------
plot_in_time(yt_sc_env, fs, 'DSBSC demodulated signal in time domain using envelope detector');
plot_in_frequency(real(yf_sc_env), fs, 'DSBSC demodulated signal in frequency domain using envelope detector');
%-------------------------------------------------sound--------------------------------------------------------------------
%sound(yt_sc_env,fs); %should not be good




%#################################################################################################################################################################
%-------------------------------------------------DSBTC no snr envelope detector demodulation--------------------------------------------------------------------
[yt_tc_env, yf_tc_env] = env_demod(yt_tc,fs_res,fs,0,0);
%-------------------------------------------------plotting--------------------------------------------------------------------
plot_in_time(yt_tc_env, fs, 'DSBTC demodulated signal no snr, in time domain using envelope detector');
plot_in_frequency(real(yf_tc_env), fs, 'DSBTC demodulated signal no snr in frequency domain using envelope detector');
%-------------------------------------------------sound--------------------------------------------------------------------
%sound(yt_tc_env,fs); %should be good


%-------------------------------------------------DSBTC 0 snr envelope detector demodulation--------------------------------------------------------------------
[yt_tc_env, yf_tc_env] = env_demod(yt_tc,fs_res,fs,1,0);
%-------------------------------------------------plotting--------------------------------------------------------------------
plot_in_time(yt_tc_env, fs, 'DSBTC demodulated signal snr=0, in time domain using envelope detector');
plot_in_frequency(real(yf_tc_env), fs, 'DSBTC demodulated signal snr=0, in frequency domain using envelope detector');
%-------------------------------------------------sound--------------------------------------------------------------------
%sound(yt_tc_env,fs); 


%-------------------------------------------------DSBTC 10 snr envelope detector demodulation--------------------------------------------------------------------
[yt_tc_env, yf_tc_env] = env_demod(yt_tc,fs_res,fs,1,10);
%-------------------------------------------------plotting--------------------------------------------------------------------
plot_in_time(yt_tc_env, fs, 'DSBTC demodulated signal snr=10, in time domain using envelope detector');
plot_in_frequency(real(yf_tc_env), fs, 'DSBTC demodulated signal snr=10, in frequency domain using envelope detector');
%-------------------------------------------------sound--------------------------------------------------------------------
%sound(yt_tc_env,fs); 

%-------------------------------------------------DSBTC 30 snr envelope detector demodulation--------------------------------------------------------------------
[yt_tc_env, yf_tc_env] = env_demod(yt_tc,fs_res,fs,1,30);
%-------------------------------------------------plotting--------------------------------------------------------------------
plot_in_time(yt_tc_env, fs, 'DSBTC demodulated signal snr=30, in time domain using envelope detector');
plot_in_frequency(real(yf_tc_env), fs, 'DSBTC demodulated signal snr=30, in frequency domain using envelope detector');
%-------------------------------------------------sound--------------------------------------------------------------------
%sound(yt_tc_env,fs); 
%#################################################################################################################################################################






































%#################################################################################################################################################################
%-------------------------------------------------DSBSC no snr coherent demodulation--------------------------------------------------------------------
[yt_sc_coh, yf_sc_coh] = coh_demod(yt_sc,fs_res,fs,0,0,fc,0,f_filter);
%-------------------------------------------------plotting--------------------------------------------------------------------
plot_in_time(real(yt_sc_coh), fs, 'DSBSC demodulated signal no snr in time domain using coherent detector');
plot_in_frequency(real(yf_sc_coh), fs, 'DSBSC demodulated signal no snr in frequency domain using coherent detector');
%-------------------------------------------------sound--------------------------------------------------------------------
% sound(real(yt_sc_coh),fs); %should be good

%-------------------------------------------------DSBSC 0 snr coherent demodulation--------------------------------------------------------------------
[yt_sc_coh, yf_sc_coh] = coh_demod(yt_sc,fs_res,fs,1,0,fc,0,f_filter);
%-------------------------------------------------plotting--------------------------------------------------------------------
plot_in_time(real(yt_sc_coh), fs, 'DSBSC demodulated signal snr=0, in time domain using coherent detector');
plot_in_frequency(real(yf_sc_coh), fs, 'DSBSC demodulated signal snr=0, in frequency domain using coherent detector');
%-------------------------------------------------sound--------------------------------------------------------------------
%sound(real(yt_sc_coh),fs); %should be good


%-------------------------------------------------DSBSC 10 snr coherent demodulation--------------------------------------------------------------------
[yt_sc_coh, yf_sc_coh] = coh_demod(yt_sc,fs_res,fs,1,10,fc,0,f_filter);
%-------------------------------------------------plotting--------------------------------------------------------------------
plot_in_time(real(yt_sc_coh), fs, 'DSBSC demodulated signal snr=10, in time domain using coherent detector');
plot_in_frequency(real(yf_sc_coh), fs, 'DSBSC demodulated signal snr=10, in frequency domain using coherent detector');
%-------------------------------------------------sound--------------------------------------------------------------------
%sound(real(yt_sc_coh),fs); 


%-------------------------------------------------DSBSC 30 snr coherent demodulation--------------------------------------------------------------------
[yt_sc_coh, yf_sc_coh] = coh_demod(yt_sc,fs_res,fs,1,30,fc,0,f_filter);
%-------------------------------------------------plotting--------------------------------------------------------------------
plot_in_time(real(yt_sc_coh), fs, 'DSBSC demodulated signal snr=30, in time domain using coherent detector');
plot_in_frequency(real(yf_sc_coh), fs, 'DSBSC demodulated signal snr=30 in frequency domain using coherent detector');
%-------------------------------------------------sound--------------------------------------------------------------------
%sound(real(yt_sc_coh),fs);
%#################################################################################################################################################################






%-------------------------------------------------DSBSC no snr coherent demodulation with fc error--------------------------------------------------------------------
[yt_sc_coh, yf_sc_coh] = coh_demod(yt_sc,fs_res,fs,0,0,fc+100,0,f_filter);
%-------------------------------------------------plotting--------------------------------------------------------------------
plot_in_time(real(yt_sc_coh), fs, 'DSBSC demodulated signal in time domain using coherent detector with fc error');
plot_in_frequency(real(yf_sc_coh), fs, 'DSBSC demodulated signal in frequency domain using coherent detector with fc error');
%-------------------------------------------------sound--------------------------------------------------------------------
% sound(real(yt_sc_coh),fs); 









%-------------------------------------------------DSBSC no snr coherent demodulation with phase shift --------------------------------------------------------------------
phaseshift = degtorad(20);
[yt_sc_coh, yf_sc_coh] = coh_demod(yt_sc,fs_res,fs,0,0,fc,phaseshift,f_filter);
%-------------------------------------------------plotting--------------------------------------------------------------------
plot_in_time(real(yt_sc_coh), fs, 'DSBSC demodulated signal in time domain using coherent detector with phase shift');
plot_in_frequency(real(yf_sc_coh), fs, 'DSBSC demodulated signal in frequency domain using coherent detector with phase shift');
%-------------------------------------------------sound--------------------------------------------------------------------
%sound(real(yt_sc_coh),fs); 












