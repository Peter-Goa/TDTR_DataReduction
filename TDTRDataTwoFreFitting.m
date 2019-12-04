function [Result] = TDTRDataTwoFreFitting(raw_data_set, config)
% Fit the MultiFrequency TDTR Data to get wanted parameter
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Fit the MultiFrequency TDTR Data to get wanted parameter
% Fit the TDTR Data to get wanted parameter
% Reference: main_TDTR_171120 in TDTR_Iwamoto_171120
% Author: RL
% Date: Dec. 4, 2019
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

Ndata = length(config.Data);
% get cost function for each data
isError = 0;
tau_data = cell(1, Ndata);
fun_data = cell(1, Ndata);
Result = cell(1, Ndata);
NVars = cell(1, Ndata);
value_0 = cell(1, Ndata);
value_0_unit = cell(1, Ndata);
lb = cell(1,Ndata);
ub = cell(1,Ndata);
func = cell(1,Ndata);
for index = 1:1:Ndata
    [tau_data{index}, fun_data{index}] = dealRawData(raw_data_set{index}, config.Data{index}.tau, config.Data{index}.mode);
    Result{index}.dealed_data.tau = tau_data{index}*1E9;
    Result{index}.dealed_data.fun = fun_data{index};
    NVars{index} = size(config.Data{index}.fit_para,1);
    % 0 means no error rised during calculation and 1(not 0) means some errors rised.
    value_0{index} = config.Data{index}.fit_para(:, 3)';
    value_0_unit{index} = value_0;
    for index_1 = 1:1:NVars{index}
        switch config.Data{index}.fit_para(index_1,2)
            case 1
                value_0_unit{index}(index_1) = value_0{index}(index_1);
            case 2
                value_0_unit{index}(index_1) = value_0{index}(index_1);
            case 3
                value_0_unit{index}(index_1) = value_0{index}(index_1)/1E6;
            case 4
                value_0_unit{index}(index_1) = value_0{index}(index_1)*1E9;
            case 5
                value_0_unit{index}(index_1) = value_0{index}(index_1)/1E6;
        end
    end
    lb{index} = config.Data{index}.fit_para(:, 4)'./value_0;
    ub{index} = config.Data{index}.fit_para(:, 5)'./value_0;
    %func = @(beta) getDevOfT_P(k_0, beta(1)*ky_1, beta(2)*kxy_1, k_2, Cv_0, beta(3)*Cv_1, Cv_2, d, freq, b, beta(4)*R, l, T_P_Exp);
    func{index} = @(beta) Costfunction_assist(beta, tau_data{index}, fun_data{index}, config.Data{index});
end
options = gaoptimset('Display','iter','UseParallel', true,'Generations',config.iteration,'TolCon',1E-9);
%%%%%%%%%%%%%%here%%%%%%%%%%%%%%%

end