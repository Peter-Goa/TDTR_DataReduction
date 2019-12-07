function [cost] = CostfunctionTwoFreq_assist(beta, tau_data, fun_data, config)
% Calculate the cost for a special beta
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% CostfunctionTwoFreq_assist
% Calculate the cost for a special beta
% Reference: main_TDTR_171120 in TDTR_Iwamoto_171120
% Author: RL
% Date: Dec. 7, 2019
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
beta_set = getBetaSet(beta, config);
cost = 0;
Ndata = length(config.Data);
for index =1:1:Ndata
    cost = cost + Costfunction_assist(beta_set{index}, tau_data{index}, fun_data{index}, config.Data{index});
end