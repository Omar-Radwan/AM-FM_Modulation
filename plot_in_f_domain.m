function [ret,x_axis] = plot_in_f_domain(y,Fs,is_plot,title_label,x_label,y_label)
    ret = fftshift(fft(y));
    x_axis = linspace(-Fs/2,Fs/2,length(ret));
    if is_plot==1
        figure();
        plot(x_axis,real(ret));
        title(title_label)
        xlabel(x_label) 
        ylabel(y_label) 
    end
end
