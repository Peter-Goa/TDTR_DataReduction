function [H] = HeatTranferModel_1(kz,kr,G,d,vhc,w)
% Heat tranfer model by the method in the thesis of Aaron Jerome Schmidt (equation 3.50).

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% HeatTranferModel_1
% Using the method in the thesis of Aaron Jerome Schmidt (equation 3.50).
% Reference: norm_TDTR in TDTR_Iwamoto_171120
% Author: RL
% Date: Nov. 11, 2019
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

%Global Variable
global omega            % Matrix (1 X omegasize).
global two_way          % Mark. 0, Equation 3.50; 1, Unknown. 
global radial_mode      % Mark. 0, one dimensional heat transfer; 1, two dimentional heat transfer

alpha = kz./vhc;            % thermal diffusivity.
klimit = 10/sqrt(w*w');     % upper limit of k
dk = 10000;                 % delta k of hankel transform
k = (0:dk:klimit)';         % variabel of hankel transform
k2 = k.*k;     
omegasize = length(omega);
ksize = length(k);
Nlayer=length(kz);
q = zeros(ksize,omegasize,Nlayer);
b = zeros(ksize,omegasize,Nlayer);
c = zeros(ksize,omegasize,Nlayer);

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
        Ctemp2 = Ctemp1 + (-Ctemp/G(index)+Dtemp).*c(:,:,index);
        Dtemp2 = Ctemp1.*b(:,:,index)-Ctemp1/G(index)+Dtemp1;
        Ctemp1 = Ctemp2;
        Dtemp1 = Dtemp2;
    end
    kterm=k*ones(1,omegasize);
    expterm=exp(-k2*(w*w')/8)*ones(1,omegasize);
    fun=-kterm.*(-Dtemp1./Ctemp1).*expterm;
    H = trapz(fun)*dk;
else
    disp('Sorry. Two way model does not exist.');
    return
end