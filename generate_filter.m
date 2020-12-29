function filter = generate_filter(len,x_axis,f_filter)
    filter = ones(len,1);
    for i = 1: len
        if x_axis(i)<-f_filter || x_axis(i)>f_filter
            filter(i)=0;
        end
    end
end