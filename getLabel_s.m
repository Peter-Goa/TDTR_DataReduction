function [para_name] = getLabel_s(LayerName, fit_para)
% Assistant function to generate the layer label.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% GetLabel
% Author: RL
% Date: Dec. 25, 2019
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
length_para = size(fit_para,1);
para_name = cell(1, length_para);
LayerName = ['f', 'Laser', LayerName];
for index_para = 1:1:length_para
    para_name{index_para} = [para_name{index_para}, LayerName{fit_para(index_para,1)+2}, ' '];
    if fit_para(index_para,1) == 0
        switch fit_para(index_para, 2)
            case 1
                para_name{index_para} = [para_name{index_para}, 'Dpump'];
            case 2
                para_name{index_para} = [para_name{index_para}, 'Dprob'];
        end
    elseif fit_para(index_para,1) == -1
        para_name{index_para} = [para_name{index_para}, 'mod'];
    else
        switch fit_para(index_para, 2)
            case 1
                para_name{index_para} = [para_name{index_para}, 'kz'];
            case 2
                para_name{index_para} = [para_name{index_para}, 'kr'];
            case 3
                para_name{index_para} = [para_name{index_para}, 'vhc'];
            case 4
                para_name{index_para} = [para_name{index_para}, 'd'];
            case 5
                para_name{index_para} = [para_name{index_para}, 'G'];
        end
    end
end