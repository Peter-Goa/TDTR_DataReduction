function [cost] = Costfunction_EP(Param, tau_data, fun_data, config)
    fun = TheoryData_EP(Param, tau_data, config);
    
    %%% Sum of squares of difference between theory and data
    res = (fun - fun_data).^2;
    mediumValue = (max(fun_data)+min(fun_data))/2;
    cost = sqrt(sum(res)/length(res))/mediumValue;
end