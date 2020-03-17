function [func] = TheoryFun_EP_assist(beta, tau_data, config)
    value_0 = config.fit_para(:, 2)';
    NVars = size(config.fit_para,1);
    Param.R_gamma = config.R_gamma;
    Param.R_G = config.R_G;
    Param.R_CL = config.R_CL;
    Param.R_F = config.R_F;
    Param.R_tp = config.R_tp;
    Param.R_d = config.R_d;
    Param.R_T0 = config.R_T0;
    Param.R_tth = config.R_tth;
    Param.R_c = config.R_c;
    Param.R_omega = config.R_omega;
    Param.C_ns = config.C_ns;
    Param.R_omegap = config.R_omegap;
    Param.R_Aee = config.R_Aee;
    Param.R_Bep = config.R_Bep;
    for index = 1:1:NVars
        switch config.fit_para(index,1)
            case 1
                Param.R_G = value_0(index)*beta(index);
            case 2
                Param.R_F = value_0(index)*beta(index);
            case 3
                Param.R_tth = value_0(index)*beta(index);
            case 4
                Param.R_Aee = value_0(index)*beta(index);
            case 5
                Param.R_Bep = value_0(index)*beta(index);
        end
    end
    func = TheoryData_EP(Param, tau_data, config);
end
