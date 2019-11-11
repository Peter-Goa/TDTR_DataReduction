function [cost] = Costfunction(kz,kr,G,d,vhc,w,tau_data,fun_data)
% Cost function using least squares method in the thesis of Aaron Jerome Schmidt (equation 3.27).
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Costfunction
% Cost function using least squares method in the thesis of Aaron Jerome Schmidt (equation 3.27).
% Reference: TDTR_fit in TDTR_Iwamoto_171120
% Author: RL
% Date: Nov. 11, 2019
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

%%%%%%  global constant number
global omega_s
global k_n

H = HeatTranferModel_1(kz,kr,G,d,vhc,w);
Z = H*exp(1i*omega_s*k_n.'*tau_data); 

X = real(Z);
Y = imag(Z);

fun = swit_fun(X,Y,tau_data);

%%% Sum of squares of difference between theory and data
res = (fun - fun_data).^2;
cost=sqrt(sum(res)/length(res));
end