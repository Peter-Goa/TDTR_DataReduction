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
    global cal_para;
    % the time at which the value will be used to normalize amplitude of data [s]
    cal_para.norm_time = 1.00E-10;
    % modulation frequency [Hz]
    cal_para.omega_0 = 2*pi*config.f_mod;
    % laser reputation frequency
    cal_para.omega_s=2*pi*80.21*10^6;
    % the k in equation 3.27, which is used to consider the accumulation
    % effects
    kmax_n = 15000;
    cal_para.k_n = (-kmax_n:kmax_n);
    cal_para.omega = cal_para.omega_0+cal_para.k_n*cal_para.omega_s;
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
    tau_data=tau_raw(Nmin:Nmax,1)';
    X_data=X_fixed(Nmin:Nmax,1)';
    Y_data=Y_fixed(Nmin:Nmax,1)';
    fun_data = swit_fun(X_data,Y_data,tau_data);
    Result.dealed_data.tau = tau_data*1E9;
    Result.dealed_data.fun = fun_data;
    NVars = size(config.fit_para,1);
    % 0 means no error rised during calculation and 1(not 0) means some errors rised.
    isError = 0;
    value_0 = config.fit_para(:, 3)';
    lb = config.fit_para(:, 4)'./value_0;
    ub = config.fit_para(:, 5)'./value_0;
    %func = @(beta) getDevOfT_P(k_0, beta(1)*ky_1, beta(2)*kxy_1, k_2, Cv_0, beta(3)*Cv_1, Cv_2, d, freq, b, beta(4)*R, l, T_P_Exp);
    func = @(beta) Costfunction_assist(beta, tau_data, fun_data);

    %Genetic algorithm
    options = gaoptimset('Display','final','UseParallel', true,'Generations',config.iteration,'TolCon',1E-9);

    try
    beta = ga(func, NVars, [], [], [], [], lb, ub, [], options);
    Loss = func(beta);
    catch
        disp('We found an error during using ga function! Then we will skip this task and continue with next one.');
        isError = 1;% 0 means no error rised during calculation and 1(not 0) means 
        % some errors rised. 
        Loss = 0;
        beta = zeros(1,NVar);
    end    
    Result.StdDev = Loss;
    Result.isError = isError;
    Result.fittingValue = beta.*value_0;
    Result.theory_data.tau = tau_data;
    Result.theory_data.fun = TheoryFun_assist(beta,tau_data);
end

function [cost] = Costfunction_assist(beta, tau_data, fun_data)
    global config;
    value_0 = config.fit_para(:, 3)';
    NVars = size(config.fit_para,1);
    kz = config.kz;
    kr = config.kr;
    G = config.G;
    vhc = config.vhc;
    d = config.d;
    for index = 1:1:NVars
        switch config.fit_para(index,2)
            case 1
                kz(config.fit_para(index,1)) = value_0(index)*beta(index);
            case 2
                kr(config.fit_para(index,1)) = value_0(index)*beta(index);
            case 3
                vhc(config.fit_para(index,1)) = value_0(index)*beta(index);
            case 4
                d(config.fit_para(index,1)) = value_0(index)*beta(index);
            case 5
                G(config.fit_para(index,1)) = value_0(index)*beta(index);
        end
    end
    NLayer = size(config.kz, 2);
    for index = 1:1:NLayer
        if config.iso(index) == 0
            kr(index) = kz(index);
        end
    end
    [cost] = Costfunction(kz,kr,G,d,vhc,config.w,tau_data,fun_data);
end

function [func] = TheoryFun_assist(beta, tau_data)
    global config;
    value_0 = config.fit_para(:, 3)';
    NVars = size(config.fit_para,1);
    kz = config.kz;
    kr = config.kr;
    G = config.G;
    vhc = config.vhc;
    d = config.d;
    for index = 1:1:NVars
        switch config.fit_para(index,2)
            case 1
                kz(config.fit_para(index,1)) = value_0(index)*beta(index);
            case 2
                kr(config.fit_para(index,1)) = value_0(index)*beta(index);
            case 3
                vhc(config.fit_para(index,1)) = value_0(index)*beta(index);
            case 4
                d(config.fit_para(index,1)) = value_0(index)*beta(index);
            case 5
                G(config.fit_para(index,1)) = value_0(index)*beta(index);
        end
    end
    NLayer = size(config.kz, 2);
    for index = 1:1:NLayer
        if config.iso(index) == 1
            kr(index) = kz(index);
        end
    end
    func = TheoryData(kz,kr,G,d,vhc,config.w,tau_data);
end
