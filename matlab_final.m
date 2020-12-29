%-------------------------------------------------clear--------------------------------------------------------------------
clc
close all 
clear all
%-------------------------------------------------read file--------------------------------------------------------------------
FILE = 'eric.wav';
f_c = 100000;
f_s_resampled = 5*f_c;
[y_t, f_s]= audioread(FILE);
%-------------------------------------------------message in frequency domain--------------------------------------------------------------------
y_f = fftshift(fft(y_t));
x_f = linspace(-f_s/2,f_s/2,length(y_f));
%-------------------------------------------------plotting--------------------------------------------------------------------
%plot_signal(x_f,real(y_f),'Message in frequency domain','Frequency','Value');
%-------------------------------------------------filter--------------------------------------------------------------------
filter = generate_filter(length(y_f),x_f,4000);
y_f_filtered = filter.*y_f;
%-------------------------------------------------plotting--------------------------------------------------------------------
%plot_signal(x_f,real(filter),'Filter','Frequency','Value');
%plot_signal(x_f,real(y_f_filtered),'Filtered message in frequency domain','Frequency', 'Value');
%-------------------------------------------------back to time domain--------------------------------------------------------------------
y_t_filtered= real(ifft(ifftshift(y_f_filtered)));
%-------------------------------------------------sound--------------------------------------------------------------------
%sound(y_t_filtered,f_s);
%-------------------------------------------------plotting--------------------------------------------------------------------
x_t = linspace(0,length(y_t_filtered)/f_s, length(y_t_filtered));
%plot_signal(x_t, y_t_filtered, 'Filtered message in time domain','Time','Value'); 
%-------------------------------------------------resmple--------------------------------------------------------------------
y_t_resampled = resample(y_t_filtered,f_s_resampled,f_s);
%-------------------------------------------------carrier signal--------------------------------------------------------------------
t = linspace(0, (length(y_t_resampled)-1)/f_s_resampled, length(y_t_resampled)); %(x2-x1)/(n-1) = 1/5*fc
carrier = cos(2*pi*f_c*t).';
%-------------------------------------------------DSBSC modulation--------------------------------------------------------------------
y_t_dsbsc = carrier.*y_t_resampled;
y_f_dsbsc = fftshift(fft(y_t_dsbsc));
x_t_dsbsc = t;
x_f_dsbsc = linspace(-f_s_resampled/2,f_s_resampled/2,length(y_f_dsbsc));
%-------------------------------------------------plotting--------------------------------------------------------------------
%plot_signal(x_t_dsbsc,y_t_dsbsc,'DSBSC modulated signal in time domain', 'Time', 'Value');
%plot_signal(x_f_dsbsc,real(y_f_dsbsc),'DSBSC modulated signal in frequency domain','Frequency','Value');
%-------------------------------------------------DSBTC modulation--------------------------------------------------------------------
modulation_index = 0.5;
A = max(abs(y_t_resampled));
y_t_dsbtc = A.*(1+modulation_index.*y_t_resampled).*carrier;
y_f_dsbtc = fftshift(fft(y_t_dsbtc));
x_t_dsbtc = t;
x_f_dsbtc = linspace(-f_s_resampled/2,f_s_resampled/2,length(y_f_dsbtc));
%-------------------------------------------------plotting--------------------------------------------------------------------
%plot_signal(x_t_dsbtc,y_t_dsbtc,'DSBTC modulated signal in time domain', 'Time', 'Value');
%plot_signal(x_f_dsbtc,real(y_f_dsbtc),'DSBTC modulated signal in frequency domain','Frequency','Value');
%-------------------------------------------------demodulation--------------------------------------------------------------------
