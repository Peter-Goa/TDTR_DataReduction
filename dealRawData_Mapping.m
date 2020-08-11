function [tau_data, fun_data] = dealRawData_Mapping(raw_data, tau, mode, modification_mode, modification_value)
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
    
    if modification_mode ~= 0
        % some method to find the zero point of time
        tau_raw = tau_raw - modification_value(1);
        % shift phase to make Y value have minimum skip at time 0
        phase_shift_rad = modification_value(2);
        X_fixed = X_raw.*cos(phase_shift_rad)-Y_raw.*sin(phase_shift_rad);
        Y_fixed = X_raw.*sin(phase_shift_rad)+Y_raw.*cos(phase_shift_rad);
    else
        X_fixed = X_raw;
        Y_fixed = Y_raw;
    end
       
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