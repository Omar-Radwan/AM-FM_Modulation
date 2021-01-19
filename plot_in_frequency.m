function  plot_in_frequency(yf,fs,title_label)
    f = linspace(-fs/2,fs/2,length(yf));
    figure()
    plot(f,yf./fs);
    title(title_label);
    xlabel('Frequency'); 
    ylabel('Value'); 
end