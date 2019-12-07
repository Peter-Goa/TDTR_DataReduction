function [beta_set] = getBetaSet(beta, config)
% An assistant function used to get beta for each data
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% getBetaSet
% get beta for each data
% Author: RL
% Date: Dec. 5, 2019
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
Ndata = length(config.Data);
NCom = size(config.commonVal,1);
beta_set = cell(1, Ndata);
beta_set_temp = cell(1, Ndata);
for index = 1:1:Ndata
    NVars{index} = size(config.Data{index}.fit_para,1);
    beta_set{index} = zeros(1, NVars{index});
end
NVars_temp = NVars;
for index_1 = 1:1:NCom
    isFirst = 1;
    for index_2 = 1:1:Ndata
        if config.commonVal(index_1, index_2) ~= 0
            if isFirst == 1
                isFirst = 0;
            else
                NVars_temp{index_2} = NVars_temp{index_2} - 1;
            end
        end
    end
end
NVars_index = zeros(1,Ndata);
for index = 1:1:Ndata
    NVars_index(index+1) = NVars_index(index) + NVars_temp{index};
    beta_set_temp{index} = beta((NVars_index(index)+1):NVars_index(index+1));
end
beta_set{1} = beta_set_temp{1};
for index = 2:1:Ndata
    for index_1 = 1:1:NCom
        isFirst = 1;
        for index_2 = 1:1:index
            if config.commonVal(index_1, index_2) ~= 0
                if isFirst == 1
                    isFirst = 0;
                    Data_index = index_2;
                    Value_index = config.commonVal(index_1, index_2);
                else
                    if index_2 == index
                        beta_set{index}(config.commonVal(index_1, index_2))=beta_set{Data_index}(Value_index);
                    end
                end
            end
        end
    end
    if NVars_temp{index} ~= 0
        for index_1 = 1:1:NVars_temp{index}
            for index_2 = 1:1:NVars{index}
                if beta_set{index}(index_2) == 0
                    beta_set{index}(index_2) = beta_set_temp{index}(index_1);
                end
            end
        end
    end
end

end