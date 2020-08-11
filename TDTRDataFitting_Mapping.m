function [Result] = TDTRDataFitting(raw_data, config)
% Fit the TDTR Data to get wanted parameter
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% TDTRDataFitting
% Fit the TDTR Data to get wanted parameter
% Reference: main_TDTR_171120 in TDTR_Iwamoto_171120
% Author: RL
% Date: Nov. 14, 2019
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

    if config.modification_mode == 3
        [tau_data, fun_data] = dealRawData(raw_data, config.tau, config.mode, config.ZeroPointMode);    
    else
        [tau_data, fun_data] = dealRawData_Mapping(raw_data, config.tau, config.mode, config.modification_mode, config.modification_value);
    end
    Result.dealed_data.tau = tau_data*1E9;
    Result.dealed_data.fun = fun_data;
    NVars = size(config.fit_para,1);
    % 0 means no error rised during calculation and 1(not 0) means some errors rised.
    isError = 0;
    value_0 = config.fit_para(:, 3)';
    value_0_unit = value_0;
    for index = 1:1:NVars
        switch config.fit_para(index,2)
            case 1
                value_0_unit(index) = value_0(index);
            case 2
                value_0_unit(index) = value_0(index);
            case 3
                value_0_unit(index) = value_0(index)/1E6;
            case 4
                value_0_unit(index) = value_0(index)*1E9;
            case 5
                value_0_unit(index) = value_0(index)/1E6;
        end
    end
    lb = config.fit_para(:, 4)'./value_0;
    ub = config.fit_para(:, 5)'./value_0;
    %func = @(beta) getDevOfT_P(k_0, beta(1)*ky_1, beta(2)*kxy_1, k_2, Cv_0, beta(3)*Cv_1, Cv_2, d, freq, b, beta(4)*R, l, T_P_Exp);
    func = @(beta) Costfunction_assist(beta, tau_data, fun_data, config);

    %Genetic algorithm
%    options = gaoptimset('Display','final','UseParallel', false,'Generations',config.iteration,'TolCon',1E-9);
    options = gaoptimset('Display','iter','UseParallel', true,'Generations',config.iteration,'TolCon',1E-9);

%     try
      beta = ga(func, NVars, [], [], [], [], lb, ub, [], options);
      Loss = func(beta);
%     catch
%         disp('We found an error during using ga function!');
%         isError = 1;% 0 means no error rised during calculation and 1(not 0) means 
%         % some errors rised. 
%         Loss = 0;
%         beta = zeros(1,NVars);
%     end    
    Result.StdDev = Loss;
    Result.isError = isError;
    Result.fittingValue = beta.*value_0_unit;
    Result.theory_data.tau = tau_data*1E9;
    Result.theory_data.fun = TheoryFun_assist(beta,tau_data, config);
end