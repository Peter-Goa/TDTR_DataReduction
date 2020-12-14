function[origin_fun] = swit_fun(X,Y,time,mode)

    switch mode
        case 'r'
            origin_fun = abs(-X./Y);
        case 'p'
            origin_fun = atan(Y./X)*180/pi;
        case 'a'
            amp = sqrt(X.^2+Y.^2);
            max_amp = GetAmp(time,amp);
            origin_fun = amp/max_amp;
        case 'x'
            X_temp = abs(X);
            X_temp = X_temp - min(X_temp);
            max_x = max(X_temp);
            origin_fun = X_temp/max_x;
        case 'amp'
            amp = sqrt(X.^2+Y.^2);
            amp = amp - min(amp);
            max_amp = max(amp);
            origin_fun = amp/max_amp;
    end
end