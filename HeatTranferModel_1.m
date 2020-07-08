function [H] = HeatTranferModel_1(kz,kr,G,d,vhc,w, config, cal_para)
% Heat tranfer model by the method in the thesis of Aaron Jerome Schmidt (equation 3.50).

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% HeatTranferModel_1
% Using the method in the thesis of Aaron Jerome Schmidt (equation 3.50).
% Reference: norm_TDTR in TDTR_Iwamoto_171120
% Author: RL
% Date: Nov. 11, 2019
% Modification 1: Add two-direction heat tranfer model in the thesis of
% Pankaj Bansilal Kaul
% Data: July 2nd, 2020
% Author： RL
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

%Global Variable
omega = cal_para.omega;            % Matrix (1 X omegasize).
two_way = config.two_way;          % Mark. 0, Equation 3.50; 1, Unknown. 
radial_mode = config.radial_mode;  % Mark. 0, one dimensional heat transfer; 1, two dimentional heat transfer

alpha = kz./vhc;            % thermal diffusivity.
klimit = 10/sqrt(w*w');     % upper limit of k
dk = 1E5;                 % delta k of hankel transform
k = (0:dk:klimit)';         % variabel of hankel transform
k2 = k.*k;     
omegasize = length(omega);
ksize = length(k);
Nlayer=length(kz);
if two_way == 0
    q = zeros(ksize,omegasize,Nlayer);
    b = zeros(ksize,omegasize,Nlayer);
    c = zeros(ksize,omegasize,Nlayer);
elseif two_way == 1
% Using Feldman's model for the two-direction heat transfer. The detail of this model can be found
% in the thesis of PANKAJ BANSILAL KAUL
    eta = zeros(ksize,omegasize,Nlayer);
    gamma = zeros(ksize,omegasize,Nlayer);
    V_h_plus = zeros(ksize,omegasize);
    V_h_minus = zeros(ksize,omegasize);
    W_h_plus = zeros(ksize,omegasize);
    W_h_minus = zeros(ksize,omegasize);
    V_m_plus = zeros(ksize,omegasize);
    V_m_minus = zeros(ksize,omegasize);
    
else
% Using the method mentioned in Schmidt's thesis for the two-direction heat
% transfer.
    heat_Layer = config.heat(1,1);
    heat_dis = config.heat(1,2);
    Nlayer_a = heat_Layer;
    Nlayer_b = Nlayer - heat_Layer + 1;
    
    kz_a = zeros(1,Nlayer_a);
    kz_b = zeros(1,Nlayer_b);
    kr_a = zeros(1,Nlayer_a);
    kr_b = zeros(1,Nlayer_b);
    G_a = zeros(1,Nlayer_a - 1);
    G_b = zeros(1,Nlayer_b - 1);
    d_a = zeros(1,Nlayer_a);
    d_b = zeros(1,Nlayer_b);
    vhc_a = zeros(1,Nlayer_a);
    vhc_b = zeros(1,Nlayer_b);
    index_new = 1;
    for index = heat_Layer:-1:1
        kz_a(index_new) = kz(index);
        kr_a(index_new) = kr(index);
        if index > 1
            G_a(index_new) = G(index-1);
        end
        d_a(index_new) = d(index);
        vhc_a(index_new) = vhc(index);
        index_new = index_new + 1;
    end
    d_a(1) = heat_dis;
    index_new = 1;
    for index = heat_Layer:1:Nlayer
        kz_b(index_new) = kz(index);
        kr_b(index_new) = kr(index);
        if index < Nlayer
            G_b(index_new) = G(index);
        end
        d_b(index_new) = d(index);
        vhc_b(index_new) = vhc(index);
        index_new = index_new + 1;
    end
    d_b(1) = d(heat_Layer) - heat_dis;
    
    alpha_a = kz_a./vhc_a; 
    alpha_b = kz_b./vhc_b; 
    
    q_a = zeros(ksize,omegasize,Nlayer_a);
    b_a = zeros(ksize,omegasize,Nlayer_a);
    c_a = zeros(ksize,omegasize,Nlayer_a);
    q_b = zeros(ksize,omegasize,Nlayer_b);
    b_b = zeros(ksize,omegasize,Nlayer_b);
    c_b = zeros(ksize,omegasize,Nlayer_b);
end

if two_way == 0
    if radial_mode == 1
        for index = 1:1:Nlayer
            % Equation 3.44
            q(:,:,index) = sqrt(ones(ksize,1)*(1i*omega/alpha(index))+(kr(index)/kz(index))*k2*ones(1,omegasize));
        end
    else
        for index = 1:1:Nlayer
            % Equation 3.34
            q(:,:,index) = sqrt(ones(ksize,1)*(1i*omega/alpha(index)));
        end
    end
    for index = 1:1:Nlayer
        % Equation 3.44
        b(:,:,index) = -1./(kz(index)*q(:,:,index)).*tanh(q(:,:,index)*d(index));
        c(:,:,index) = -kz(index)*q(:,:,index).*tanh(q(:,:,index)*d(index));
    end
    Ctemp1 = c(:,:,Nlayer);
    Dtemp1 = ones(ksize,omegasize);
    for index = Nlayer-1:-1:1
        Ctemp2 = Ctemp1 + (-Ctemp1/G(index)+Dtemp1).*c(:,:,index);
        Dtemp2 = Ctemp1.*b(:,:,index)-Ctemp1/G(index)+Dtemp1;
        Ctemp1 = Ctemp2;
        Dtemp1 = Dtemp2;
    end
    kterm=k*ones(1,omegasize);
    expterm=exp(-k2*(w*w')/8)*ones(1,omegasize);
    fun=-kterm.*(-Dtemp1./Ctemp1).*expterm;
    H = trapz(fun)*dk;
    
elseif two_way == 1
    heat_Layer = config.heat(1,1);
    measure_Layer = config.heat(2,1);
    heat_dis = config.heat(1,2);
    measure_dis = config.heat(2,2);
    
    % Calculate eta and gamma
    for index = 1:1:Nlayer
        % Equation 142
        eta(:,:,index) = sqrt(ones(ksize,1)*(1i*omega/alpha(index))+(kr(index)/kz(index))*k2*ones(1,omegasize));
        gamma(:,:,index) = kz(index)*eta(:,:,index);
    end
    
    % Initialization
    V_h_plus(:,:) = 0.5*(1+gamma(:,:,1)./gamma(:,:,2)+gamma(:,:,1)/G(1));
    V_h_minus(:,:) = 0.5*(1-gamma(:,:,1)./gamma(:,:,2)+gamma(:,:,1)/G(1));
    W_h_plus(:,:) = 0.5*(1-gamma(:,:,Nlayer)./gamma(:,:,Nlayer-1)+gamma(:,:,Nlayer)/G(Nlayer-1));
    W_h_minus(:,:) = 0.5*(1+gamma(:,:,Nlayer)./gamma(:,:,Nlayer-1)+gamma(:,:,Nlayer)/G(Nlayer-1));
    V_m_plus(:,:) = 0.5*(1+gamma(:,:,1)./gamma(:,:,2)+gamma(:,:,1)/G(1));
    V_m_minus(:,:) = 0.5*(1-gamma(:,:,1)./gamma(:,:,2)+gamma(:,:,1)/G(1));
    
    % Calculating V_h
    for index = 2:1:(heat_Layer-1)
        V_h_plus(:,:) = 0.5*((1+gamma(:,:,index)./gamma(:,:,index+1)+gamma(:,:,index)/G(index)).*(V_h_plus(:,:).*exp(eta(:,:,index)*d(index))) + (1-gamma(:,:,index)./gamma(:,:,index+1)-gamma(:,:,index)/G(index)).*(V_h_minus(:,:).*exp(-eta(:,:,index)*d(index))));
        V_h_minus(:,:) = 0.5*((1-gamma(:,:,index)./gamma(:,:,index+1)+gamma(:,:,index)/G(index)).*(V_h_plus(:,:).*exp(eta(:,:,index)*d(index))) + (1+gamma(:,:,index)./gamma(:,:,index+1)-gamma(:,:,index)/G(index)).*(V_h_minus(:,:).*exp(-eta(:,:,index)*d(index))));
    end
    V_h_plus(:,:) = V_h_plus(:,:).*exp(eta(:,:,heat_Layer)*heat_dis);
    V_h_minus(:,:) = V_h_minus(:,:).*exp(-eta(:,:,heat_Layer)*heat_dis);
    
    % Calculationg W_h
    for index = (Nlayer-1):-1:(heat_Layer+1)
        %W_h_plus(:,:) = 0.5*((1+gamma(:,:,index)./gamma(:,:,index-1)-gamma(:,:,index)/G(index-1)).*(W_h_plus(:,:).*exp(-eta(:,:,index)*d(index))) + (1-gamma(:,:,index)./gamma(:,:,index-1)+gamma(:,:,index)/G(index-1)).*(W_h_minus(:,:).*exp(eta(:,:,index)*d(index))));
        %W_h_minus(:,:) = 0.5*((1-gamma(:,:,index)./gamma(:,:,index-1)-gamma(:,:,index)/G(index-1)).*(W_h_plus(:,:).*exp(-eta(:,:,index)*d(index))) + (1+gamma(:,:,index)./gamma(:,:,index-1)+gamma(:,:,index)/G(index-1)).*(W_h_minus(:,:).*exp(eta(:,:,index)*d(index))));
        W_h_plus(:,:) = 0.5*((1+gamma(:,:,index)./gamma(:,:,index-1)-gamma(:,:,index)/G(index-1)).*(W_h_plus(:,:).*exp(-2*eta(:,:,index)*d(index))) + (1-gamma(:,:,index)./gamma(:,:,index-1)+gamma(:,:,index)/G(index-1)).*(W_h_minus(:,:)));
        W_h_minus(:,:) = 0.5*((1-gamma(:,:,index)./gamma(:,:,index-1)-gamma(:,:,index)/G(index-1)).*(W_h_plus(:,:).*exp(-2*eta(:,:,index)*d(index))) + (1+gamma(:,:,index)./gamma(:,:,index-1)+gamma(:,:,index)/G(index-1)).*(W_h_minus(:,:)));
    end
%    W_h_plus(:,:) = W_h_plus(:,:).*exp(-eta(:,:,heat_Layer)*(d(heat_Layer) - heat_dis));
%    W_h_minus(:,:) = W_h_minus(:,:).*exp(eta(:,:,heat_Layer)*(d(heat_Layer) - heat_dis));
    W_h_plus(:,:) = W_h_plus(:,:).*exp(-2*eta(:,:,heat_Layer)*(d(heat_Layer) - heat_dis));
    W_h_minus(:,:) = W_h_minus(:,:);
    
    % Calculating V_m
    for index = 2:1:(measure_Layer-1)
        V_m_plus(:,:) = 0.5*((1+gamma(:,:,index)./gamma(:,:,index+1)+gamma(:,:,index)/G(index)).*(V_m_plus(:,:).*exp(eta(:,:,index)*d(index))) + (1-gamma(:,:,index)./gamma(:,:,index+1)-gamma(:,:,index)/G(index)).*(V_m_minus(:,:).*exp(-eta(:,:,index)*d(index))));
        V_m_minus(:,:) = 0.5*((1-gamma(:,:,index)./gamma(:,:,index+1)+gamma(:,:,index)/G(index)).*(V_m_plus(:,:).*exp(eta(:,:,index)*d(index))) + (1+gamma(:,:,index)./gamma(:,:,index+1)-gamma(:,:,index)/G(index)).*(V_m_minus(:,:).*exp(-eta(:,:,index)*d(index))));
    end
    V_m_plus(:,:) = V_m_plus(:,:).*exp(eta(:,:,measure_Layer)*measure_dis);
    V_m_minus(:,:) = V_m_minus(:,:).*exp(-eta(:,:,measure_Layer)*measure_dis);
    
    % Calculating H(omega)
    kterm=k*ones(1,omegasize);
    expterm=exp(-k2*(w*w')/8)*ones(1,omegasize);
    fun=-kterm.*((1./(2*gamma(:,:,heat_Layer))).*(W_h_plus(:,:)+W_h_minus(:,:)).*(V_m_plus(:,:)+V_m_minus(:,:))./(V_h_plus(:,:).*W_h_minus(:,:)-V_h_minus(:,:).*W_h_plus(:,:))).*expterm;
    H = trapz(fun)*dk;
    
else
    % for A side 
    for index = 1:1:Nlayer_a
        % Equation 3.44
        q_a(:,:,index) = sqrt(ones(ksize,1)*(1i*omega/alpha_a(index))+(kr_a(index)/kz_a(index))*k2*ones(1,omegasize));
    end
    for index = 1:1:Nlayer_a
        % Equation 3.44
        b_a(:,:,index) = -1./(kz_a(index)*q_a(:,:,index)).*tanh(q_a(:,:,index)*d_a(index));
        c_a(:,:,index) = -kz_a(index)*q_a(:,:,index).*tanh(q_a(:,:,index)*d_a(index));
    end
    Ctemp1_a = c_a(:,:,Nlayer_a);
    Dtemp1_a = ones(ksize,omegasize);
    for index = Nlayer_a-1:-1:1
        Ctemp2_a = Ctemp1_a + (-Ctemp1_a/G_a(index)+Dtemp1_a).*c_a(:,:,index);
        Dtemp2_a = Ctemp1_a.*b_a(:,:,index)-Ctemp1_a/G_a(index)+Dtemp1_a;
        Ctemp1_a = Ctemp2_a;
        Dtemp1_a = Dtemp2_a;
    end
    % for B side 
    for index = 1:1:Nlayer_b
        % Equation 3.44
        q_b(:,:,index) = sqrt(ones(ksize,1)*(1i*omega/alpha_b(index))+(kr_b(index)/kz_b(index))*k2*ones(1,omegasize));
    end
    for index = 1:1:Nlayer_b
        % Equation 3.44
        b_b(:,:,index) = -1./(kz_b(index)*q_b(:,:,index)).*tanh(q_b(:,:,index)*d_b(index));
        c_b(:,:,index) = -kz_b(index)*q_b(:,:,index).*tanh(q_b(:,:,index)*d_b(index));
    end
    Ctemp1_b = c_b(:,:,Nlayer_b);
    Dtemp1_b = ones(ksize,omegasize);
    for index = Nlayer_b-1:-1:1
        Ctemp2_b = Ctemp1_b + (-Ctemp1_b/G_b(index)+Dtemp1_b).*c_b(:,:,index);
        Dtemp2_b = Ctemp1_b.*b_b(:,:,index)-Ctemp1_b/G_b(index)+Dtemp1_b;
        Ctemp1_b = Ctemp2_b;
        Dtemp1_b = Dtemp2_b;
    end
    kterm=k*ones(1,omegasize);
    expterm=exp(-k2*(w*w')/8)*ones(1,omegasize);
    fun=-kterm.*(-Dtemp1_a.*Dtemp1_b./(Dtemp1_a.*Ctemp1_b + Dtemp1_b.*Ctemp1_a)).*expterm;
    H = trapz(fun)*dk;    
end