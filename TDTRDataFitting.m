function [Result] = TDTRDataFitting(raw_data)
% Fit the TDTR Data to get wanted parameter
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% TDTRDataFitting
% Fit the TDTR Data to get wanted parameter
% Reference: main_TDTR_171120 in TDTR_Iwamoto_171120
% Author: RL
% Date: Nov. 14, 2019
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

    global config;
    % the fitting method is a matlab built-in function named ganatic algorithm
    % do some processing on the raw data
    tau_raw = raw_data(:,1)*1E-9;
    X_raw = raw_data(:,2);
    Y_raw = raw_data(:,3);
    % make the time whose X value is maximum is zero
    [~, tau_zero_index] = max(X_raw);
    tau_zero = tau_raw(tau_zero_index);
    tau_raw = tau_raw - tau_zero;
    % shift phase to make Y value have minimum skip at time 0
    options_phase = optimset('Display','iter','MaxFunEvals',200,'TolFun',1e-14,'TolX',1e-14);
    phase_shift_rad = 0.1;
    phase_shift_rad = fminsearch(@(phase_shift_rad) SetPhaseRatio(tau_raw, X_raw, Y_raw, phase_shift_rad, tau_zero_index), phase_shift_rad, options_phase);
    X_fixed = X_raw.*cos(phase_shift_rad)-Y_raw.*sin(phase_shift_rad);
    Y_fixed = X_raw.*sin(phase_shift_rad)+Y_raw.*cos(phase_shift_rad);
    Ndata=length(tau_raw);
    Nmin=Ndata;
    for i=Ndata:-1:1
        if tau_raw(i)>=config.tau(1)
            Nmin=i;
        end
    end
    Nmax=Nmin;
    for i=Nmin:Ndata
        if tau_raw(i)<=config.tau(2)
            Nmax=i;
        end
    end

    %------ extract min to max part from raw_data -----------
    tau_data = tau_raw(Nmin:Nmax,1)';
    X_data = X_fixed(Nmin:Nmax,1)';
    Y_data = Y_fixed(Nmin:Nmax,1)';
    fun_data = swit_fun(X_data,Y_data,tau_data, config);
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