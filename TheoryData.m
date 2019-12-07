function [func] = TheoryData(kz,kr,G,d,vhc,w,tau_data, config)
    % modulation frequency [Hz]
    cal_para.omega_0 = 2*pi*config.f_mod;
    % laser reputation frequency
    cal_para.omega_s=2*pi*80.21*10^6;
    % the k in equation 3.27, which is used to consider the accumulation
    % effects
    kmax_n = 5000;
    cal_para.k_n = (-kmax_n:kmax_n);
    cal_para.omega = cal_para.omega_0+cal_para.k_n*cal_para.omega_s;
    H = HeatTranferModel_1(kz,kr,G,d,vhc,w, config, cal_para);
    Z = H*exp(1i*cal_para.omega_s*cal_para.k_n'*tau_data); 

    X = real(Z);
    Y = imag(Z);

    func = swit_fun(X,Y,tau_data, config.mode);
end