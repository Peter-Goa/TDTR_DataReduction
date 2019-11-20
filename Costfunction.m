function [cost] = Costfunction(kz,kr,G,d,vhc,w,tau_data,fun_data, config)
% Cost function using least squares method in the thesis of Aaron Jerome Schmidt (equation 3.27).
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Costfunction
% Cost function using least squares method in the thesis of Aaron Jerome Schmidt (equation 3.27).
% Reference: TDTR_fit in TDTR_Iwamoto_171120
% Author: RL
% Date: Nov. 11, 2019
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

fun = TheoryData(kz,kr,G,d,vhc,w,tau_data, config);

%%% Sum of squares of difference between theory and data
res = (fun - fun_data).^2;
cost = sqrt(sum(res)/length(res));
end