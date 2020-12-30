function [yt_demod,yf_demod] = env_demod(st,fs_cur,fs_res,is_snr,snr)
    
    if is_snr == 1
        st = awgn(st, snr);
    end
    
    yt_demod = resample(abs(hilbert(st)),fs_res,fs_cur); %envelope detector and resample 
    yf_demod = fftshift(fft(yt_demod));
end