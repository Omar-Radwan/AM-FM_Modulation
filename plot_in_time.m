function  plot_in_time(yt,fs,title_label)
    t = linspace(0, length(yt)/fs, length(yt));
    figure()
    plot(t,yt);
    title(title_label);
    xlabel('Time'); 
    ylabel('Value'); 
end