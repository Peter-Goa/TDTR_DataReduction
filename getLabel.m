function [para_name, format_f] = getLabel(LayerName, fit_para)
% Assistant function to generate the layer label.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% GetLabel
% Author: RL
% Date: Dec. 4, 2019
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
length_para = size(fit_para,1);
para_name = '';
LayerName = ['Laser', LayerName];
for index_para = 1:1:length_para
    para_name = [para_name, LayerName{fit_para(index_para,1)+1}, '_'];
    if fit_para(index_para,1) == 0
        switch fit_para(index_para, 2)
            case 1
                para_name = [para_name, 'Dpump [um]'];
            case 2
                para_name = [para_name, 'Dprob [um]'];
        end
    else
        switch fit_para(index_para, 2)
            case 1
                para_name = [para_name, 'kz [W/mK]'];
            case 2
                para_name = [para_name, 'kr [W/mK]'];
            case 3
                para_name = [para_name, 'vhc [MJ/m^3K]'];
            case 4
                para_name = [para_name, 'd [nm]'];
            case 5
                para_name = [para_name, 'G [MW/m^2K]'];
        end
    end
    para_name = [para_name, '   ']; 
end

format_f = '';
for index = 1:1:length_para
    format_f = [format_f '%f  '];
end
format_f = [format_f '\r\n'];