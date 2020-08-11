function [Delta_tau, Delta_phase] = dealTDTRRawData(raw_data, ZeroPointMode)
% Do some processing on the raw data
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% dealRawData
% Do some processing on the raw data
% Reference: main_TDTR_171120 in TDTR_Iwamoto_171120
% Author: RL
% Date: Aug. 5, 2020
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

    % do some processing on the raw data
    tau_raw = raw_data(:,1)*1E-9;
    X_raw = raw_data(:,2);
    Y_raw = raw_data(:,3);
    % some method to find the zero point of time
    switch ZeroPointMode
        case 1
            % make the time whose dX value is maximum is zero
            % It's possible the dX is nagetie, so we add abs to make sure the programe run correctly
            d_X_raw = (X_raw(2 : end) - X_raw(1 : end-1))./(tau_raw(2 : end) - tau_raw(1 : end-1));
            [~, tau_zero_index] = max(abs(d_X_raw));
            tau_zero_index = tau_zero_index + 1;
        otherwise
            % make the time whose X value is maximum is zero
            % It's possible the X is nagetie, so we add abs to make sure the programe run correctly
            [~, tau_zero_index] = max(abs(X_raw));
    end
    
    Delta_tau = tau_raw(tau_zero_index);
    
    % shift phase to make Y value have minimum skip at time 0
    options_phase = optimset('Display','iter','MaxFunEvals',200,'TolFun',1e-14,'TolX',1e-14);
    phase_shift_rad = 0.1;
    Delta_phase = fminsearch(@(phase_shift_rad) SetPhaseRatio(tau_raw, X_raw, Y_raw, phase_shift_rad, tau_zero_index), phase_shift_rad, options_phase);
