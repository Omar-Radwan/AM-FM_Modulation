function filter = generate_filter(signal_length,fs,f_filter)
    filter = ones(signal_length,1);
    f = linspace(-fs/2,fs/2,signal_length);
    for i = 1: signal_length
        if abs(f(i))>f_filter
            filter(i)=0;
        end
    end
end