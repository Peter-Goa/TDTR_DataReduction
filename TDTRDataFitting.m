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
%------->Here




% the fitting method is a matlab built-in function named ganatic algorithm
% do some processing on the raw data
tau_raw = raw_data(:,1)*1E-9;
tau_raw_ns = raw_data(:,1);
X_raw = raw_data(:,2);
Y_raw = raw_data(:,3);
fun_raw = swit_fun(X_raw,Y_raw,tau_raw);
multiple_control = [0,0];
multiple_UB = [3,4];
multiple_LB = [-2,-4];
% 0 means the calculate is rejected and 1(not 0) means the calcualte is pass 
isPass = 0;
% 0 means no error rised during calculation and 1(not 0) means some errors rised.
isError = 0; 

% the method which decides whether the result pass or not concludes several 
% steps:
% 1. the acceptable range of the result is [low_band*1.3, high_band*0.9]
% 2. if the result did not pass, the center value should multiply 5 or 1/5
% 3. the range used in ganatic algorithm is [center_value*0.3,center_value*3]
while isPass == 0
    % ********************************************************************
    % Loss function
    % this section should be rewrited if you want to calculate other value
    % in your code!!!
    % ********************************************************************
    ky_1 = Meas_Conf.ky_1*(5^multiple_control(1));
    R = Meas_Conf.R*(5^multiple_control(2));
    %func = @(beta) getDevOfT_P(k_0, beta(1)*ky_1, beta(2)*kxy_1, k_2, Cv_0, beta(3)*Cv_1, Cv_2, d, freq, b, beta(4)*R, l, T_P_Exp);
    func = @(beta) getDevOfT_P(k_0, beta(1)*ky_1, kxy_1, k_2, Cv_0, Cv_1, Cv_2, d_1, d_2, freq, b, beta(2)*R, l, T_P_Exp);

    %Genetic algorithm
    options = gaoptimset('Display','final','UseParallel', true,'Generations',config.iteration,'TolCon',1E-9);
    lb = [0.3 0.3];
    ub = [3 3];
    NVars = 2;
    try
    beta = ga(func, NVars, [], [], [], [], lb, ub, [], options);
    catch
        disp('We found an error during using ga function! Then we will skip this task and continue with next one.');
        isError = 1;% 0 means no error rised during calculation and 1(not 0) means 
        % some errors rised. 
        isPass = 0;
        Loss = 0;
        beta = [0 0];
        break
    end
    if (beta(1)>0.4)&&(beta(1)<2.7)&&(beta(2)>0.4)&&(beta(2)<2.7)
        isPass = 1;
    else
        if beta(1)<0.4
            multiple_control(1) = multiple_control(1) - 1;
        end
        if beta(1)>2.7
            multiple_control(1) = multiple_control(1)+1;
        end
        if beta(2)<0.4
            multiple_control(2) = multiple_control(2)-1;
        end
        if beta(2)>2.7
            multiple_control(2) = multiple_control(2)+1;
        end
        isPass = 0;
    end
    if (multiple_control(1)<multiple_LB(1))||(multiple_control(1)>multiple_UB(1))||(multiple_control(2)<multiple_LB(2))||(multiple_control(2)>multiple_UB(2))
        isPass = 0;
        disp('The method can not find an eligible solution!!!');
        break
    end
    Loss = func(beta);
end
result.Loss = Loss;
result.isPass = isPass;
result.isError = isError;
result.value = [beta(1)*ky_1,beta(2)*R];