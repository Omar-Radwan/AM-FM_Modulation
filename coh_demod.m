function [yt_demod, yf_demod] = coh_demod(st,fs_bef,fs_aft,is_snr,snr,fc,phase,f_filter)
if is_snr == 1
    st = awgn(st, snr);
end

t = linspace(0,length(st)/fs_bef, length(st)); %(x2-x1)/(n-1) = 1/5*fc, linspace(x1,x2,n)
carrier_t = cos(2*pi*fc*t+phase).';
filter = generate_filter(length(st),fs_bef,f_filter); %filter in frequency domain
tmp = st.*carrier_t; %m(t)*c(t)
tmp = real(fftshift(fft(tmp))).*filter;  
tmp = ifft(ifftshift(tmp));
yt_demod = resample(tmp,fs_aft,fs_bef); %mutiply by carrier,filter and resample
yf_demod = fftshift(fft(yt_demod));
end