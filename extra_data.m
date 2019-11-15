function [extra_t,extra_f] = extra_data( t_data,f_data,t_min,t_max)

    Ndata=length(t_data);
    Nmin=Ndata;
    for i=Ndata:-1:1
        if t_data(i)>=t_min
            Nmin=i;
        end
    end
    Nmax=Nmin;
    for i=Nmin:Ndata
        if t_data(i)<=t_max
            Nmax=i;
        end
    end
    extra_t=t_data(Nmin:Nmax);
    extra_f=f_data(Nmin:Nmax);

end