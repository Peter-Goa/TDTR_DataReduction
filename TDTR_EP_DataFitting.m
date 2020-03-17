function [Result] = TDTR_EP_DataFitting(raw_data, config)
% Deal with the data for electron-phonon coupling measurement 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Deal with the data for electron-phonon coupling measurement with
% two-temperature Model and Drude Model
% Reference:
% 1. Re-examining Electron-Fermi Relaxation in Gold Films With a Nonlinear
% Thermoreflectance Model
% 2. Influence of intraband transitions on the electron thermalreflectance
% response of metals
% 3. Effects of electron-boundary scattering on changes in
% thermoreflectance in thin metal films undergoing intraband excitations
% Components of output
% 1. Dealed data: dealed_data.tau, dealed_data.fun
% 2. Theory data: theory_data.tau, theory_data.fun
% 3. Fitting result: fittingValue, StdDev
% Author: RL
% Date: March 17, 2020
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
    [tau_data, fun_data] = dealRawData(raw_data, config.tau, config.mode, config.ZeroPointMode);
    Result.dealed_data.tau = tau_data*1E12; % ps
    Result.dealed_data.fun = fun_data;
    NVars = size(config.fit_para,1);
    % 0 means no error rised during calculation and 1(not 0) means some errors rised.
    isError = 0;
    value_0 = config.fit_para(:, 2)';
    value_0_unit = value_0;    
    for index = 1:1:NVars
        switch config.fit_para(index,1)
            case 1
                value_0_unit(index) = value_0(index)/1E16;
            case 2
                value_0_unit(index) = value_0(index);
            case 3
                value_0_unit(index) = value_0(index)*1E12;
            case 4
                value_0_unit(index) = value_0(index)/1E7;
            case 5
                value_0_unit(index) = value_0(index)/1E11;
        end
    end
    lb = config.fit_para(:, 3)'./value_0;
    ub = config.fit_para(:, 4)'./value_0;
    %func = @(beta) getDevOfT_P(k_0, beta(1)*ky_1, beta(2)*kxy_1, k_2, Cv_0, beta(3)*Cv_1, Cv_2, d, freq, b, beta(4)*R, l, T_P_Exp);
    func = @(beta) Costfunction_EP_assist(beta, tau_data, fun_data, config);

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
    Result.theory_data.tau = tau_data*1E12; %ps
    Result.theory_data.fun = TheoryFun_EP_assist(beta,tau_data, config);    
end