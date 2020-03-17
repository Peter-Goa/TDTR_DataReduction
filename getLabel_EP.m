function [para_name, format_f] = getLabel_EP(fit_para)
    length_para = size(fit_para, 1);
    para_name = '';
    format_f = '';
    for index_para = 1:1:length_para
        switch fit_para(index_para, 1)
            case 1 % electron-phonon coupling factor
                para_name = [para_name, 'G [1E16 Wm-3K-1]'];
            case 2 % absorpted incident laser fluence
                para_name = [para_name, 'F [Jm-2]'];
            case 3 % delay in electron thermailization after pulse absorption
                para_name = [para_name, 'tth [ps]'];
            case 4 % constant coefficient 1 of Drude model
                para_name = [para_name, 'Aee [1E7 K-2s-1]'];
            case 5 % constant coefficient 2 of Drude model
                para_name = [para_name, 'Bep [1E11 K-1s-1]'];
        end
        para_name = [para_name, '   '];
        format_f = [format_f '%f  '];
    end
    format_f = [format_f '\r\n'];
end