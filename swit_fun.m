function[origin_fun] = swit_fun(X,Y,time)
    global config

    switch config.mode
        case 'r'
            origin_fun = -X./Y;
        case 'p'
            origin_fun = atan(Y./X)*180/pi;
        case 'a'
            amp = sqrt(X.^2+Y.^2);
            max_amp = GetAmp(time,amp);
            origin_fun = amp/max_amp;
    end
end