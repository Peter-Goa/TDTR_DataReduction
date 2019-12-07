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
NCom = size(config.commonVal,1);
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
    value_0_unit{index} = value_0{index};
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
    lb{index} = config.Data{index}.fit_para(:, 4)'./value_0{index};
    ub{index} = config.Data{index}.fit_para(:, 5)'./value_0{index};
    func_set{index} = @(beta) Costfunction_assist(beta, tau_data{index}, fun_data{index}, config.Data{index});
end
options = gaoptimset('Display','iter','UseParallel', true,'Generations',config.iteration,'TolCon',1E-9);

% To combine funcs to get a total func
% Calculate the total lb and ub
value_0_temp = value_0;
value_0_unit_temp = value_0_unit;
lb_temp = lb;
ub_temp = ub;
NVars_temp = NVars;
lb_total = [];
ub_total = [];
value_0_total = [];
value_0_unit_total = [];
omit_index = cell(1, Ndata);
for index = 1:1:Ndata
    omit_index{index} = [];
end
for index_1 = 1:1:NCom
    isFirst = 1;
    for index_2 = 1:1:Ndata
        if config.commonVal(index_1, index_2) ~= 0
            if isFirst == 1
                isFirst = 0;
            else
                omit_index{index_2} = [omit_index{index_2} config.commonVal(index_1, index_2)];
                NVars_temp{index_2} = NVars_temp{index_2} - 1;
            end
        end
    end
end
for index = 1:1:Ndata
    value_0_temp{index}(omit_index{index}) = [];
    value_0_unit_temp{index}(omit_index{index}) = [];
    lb_temp{index}(omit_index{index}) = [];
    ub_temp{index}(omit_index{index}) = [];
end
NVars_total = 0;
for index = 1:1:Ndata
    lb_total = [lb_total lb_temp{index}];
    ub_total = [ub_total ub_temp{index}];
    % Calculate the total NVars
    NVars_total = NVars_total + NVars_temp{index};
end
% Calculate func_total
NVars_index = zeros(1,Ndata);
for index = 1:1:Ndata
    NVars_index(index+1) = NVars_index(index) + NVars_temp{index};
end
func = @(beta) CostfunctionTwoFreq_assist(beta, tau_data, fun_data, config);
beta = ga(func, NVars_total, [], [], [], [], lb_total, ub_total, [], options);
beta_set = getBetaSet(beta, config);
for index = 1:1:Ndata
    Loss = func_set{index}(beta_set{index});
    Result{index}.StdDev = Loss;
    Result{index}.isError = 0;
    Result{index}.fittingValue = beta_set{index}.*value_0_unit{index};
    Result{index}.theory_data.tau = tau_data{index}*1E9;
    Result{index}.theory_data.fun = TheoryFun_assist(beta_set{index},tau_data{index}, config.Data{index});  
end