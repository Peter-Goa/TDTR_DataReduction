function [tau_data, fun_data] = dealRawData(raw_data, tau, mode)
% Do some processing on the raw data
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% dealRawData
% Do some processing on the raw data
% Reference: main_TDTR_171120 in TDTR_Iwamoto_171120
% Author: RL
% Date: Dec. 4, 2019
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

    % do some processing on the raw data
    tau_raw = raw_data(:,1)*1E-9;
    X_raw = raw_data(:,2);
    Y_raw = raw_data(:,3);
    % make the time whose X value is maximum is zero
    % It's possible the X is nagetie, so we add abs to make sure the programe run correctly 
    [~, tau_zero_index] = max(abs(X_raw));
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
        if tau_raw(i)>=tau(1)
            Nmin=i;
        end
    end
    Nmax=Nmin;
    for i=Nmin:Ndata
        if tau_raw(i)<=tau(2)
            Nmax=i;
        end
    end

    %------ extract min to max part from raw_data -----------
    tau_data = tau_raw(Nmin:Nmax,1)';
    X_data = X_fixed(Nmin:Nmax,1)';
    Y_data = Y_fixed(Nmin:Nmax,1)';
    fun_data = swit_fun(X_data,Y_data,tau_data, mode);
end