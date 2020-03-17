function [fun] = TheoryData_EP(Param, tau_data, config)
    y0 = [Param.R_T0, Param.R_T0];
    [t, y] = ode45(@(t,y) ttm(t,y, Param), [0, 2*(config.tau(2)-config.tau(1))], y0);
    Te = y(:,1); 
    Tl = y(:,2);
    Te = [Param.R_T0; Te];
    Tl = [Param.R_T0; Tl];
    omega_tau = Param.R_Aee*Te.^2 + Param.R_Bep*Tl;
    C_epsilon = 1 - Param.R_omegap^2./(Param.R_omega*(Param.R_omega+1i*omega_tau));
    R_eps_r = real(C_epsilon);
    R_eps_i = imag(C_epsilon);
    R_nf_r = (1/(sqrt(2)))*((R_eps_r.^2+R_eps_i.^2).^(1/2)+R_eps_r).^(1/2);
    R_nf_i = (1/(sqrt(2)))*((R_eps_r.^2+R_eps_i.^2).^(1/2)-R_eps_r).^(1/2);
    C_nf = R_nf_r + 1i*R_nf_i;
    C_delta = Param.R_omega*Param.R_d*C_nf/Param.R_c;
    M11 = cos(C_delta);
    M12 = (-1i./C_nf).*sin(C_delta);
    M21 = -1i*C_nf.*sin(C_delta);
    M22 = cos(C_delta);
    C_r = ((M11+Param.C_ns*M12)-(M21+Param.C_ns*M22))./((M11+Param.C_ns*M12)+(M21+Param.C_ns*M22));
    R_Rf = conj(C_r).*C_r;
    Rf_T0 = R_Rf(1);
    R_Rf = R_Rf(2:end);
    R_deltaRf = abs((R_Rf - Rf_T0)/Rf_T0);
    [R_deltaRf_max, index] = max(R_deltaRf);
    R_deltaRf_Norm = R_deltaRf/R_deltaRf_max;
    t = t - t(index);
    fun = interp1(t, R_deltaRf_Norm, tau_data, 'spline');
end


function dy = ttm(t,y, Param)
    a = (1/(Param.R_gamma*y(1)))*(-Param.R_G*(y(1)-y(2))+(Param.R_F/(Param.R_d*(Param.R_tp+Param.R_tth)))*exp(-2.77*((t-2*(Param.R_tp+Param.R_tth))/(Param.R_tp+Param.R_tth))^2));
    b = (1/(Param.R_CL))*(Param.R_G*(y(1)-y(2)));
    dy = [a;b];
end
