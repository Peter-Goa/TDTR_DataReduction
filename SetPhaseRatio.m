function [res] = SetPhaseRatio(tau_raw, X_raw, Y_raw, phase_shift, tau_zero_index)
% the function is used to calculate the difference between average of Y in
% the range of tau_range_min to zero and the average of Y in the range of 0
% to tau_range_max

tau_range_min = -20*10^-12;
tau_range_max = 80*10^-12;

Y_fixed = X_raw.*sin(phase_shift)+Y_raw.*cos(phase_shift);

t_data_minus = tau_raw(1:tau_zero_index-2);
t_data_plus = tau_raw(tau_zero:end);

Y_data_minus = Y_fixed(1:tau_zero_index-2);
Y_data_plus = Y_fixed(tau_zero_index:end);

[~,Y_omit_data_minus]=extra_data(t_data_minus,Y_data_minus,tau_range_min,tau_range_max);
[~,Y_omit_data_plus]=extra_data(t_data_plus,Y_data_plus,tau_range_min,tau_range_max);

res=abs(mean(Y_omit_data_plus)-mean(Y_omit_data_minus));

end

