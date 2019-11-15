function [func] = TheoryData(kz,kr,G,d,vhc,w,tau_data)
    global cal_para;

    H = HeatTranferModel_1(kz,kr,G,d,vhc,w);
    Z = H*exp(1i*cal_para.omega_s*cal_para.k_n.'*tau_data); 

    X = real(Z);
    Y = imag(Z);

    func = swit_fun(X,Y,tau_data);
end