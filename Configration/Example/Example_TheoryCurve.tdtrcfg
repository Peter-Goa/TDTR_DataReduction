global config

kal=238;    val=2.42;
kpdds=15; vpdds=3.38;
kmgo=41;  vmgo=3.4;

% cross-plane thermal conductivity of each layer [W/mK]
config.kz = [kal, kpdds, kmgo];

% in-plane thermal conductivity of each layer [W/mK]
config.kr = [kal, kpdds, kmgo];

% Interfacial thermal conductance between each layer[W/m^2K]
config.G = [100, 100]*1E6;

% volume specific heat of each layer [J/m^3K]
config.vhc = [val, vpdds vmgo]*1E6;

% thickness of each layer [m]
config.d = [23.8, 434.1, 500E3]*1E-9;

% diameter of laser beam ([pump probe])[m]
config.w = [15.5 10.9]/2*1E-6;

% 0 - the laser applying directly on metal transducer on the top of sample
% 1 - the applying form glass side
config.two_way = 0;

% 0 - one-dimentional heat tranfer
% 1 - two-dimentional heat tranfer
config.radial_mode = 1;

% modulation frequence of pump laser (1.111, 4.721, 11.05 MHz) [Hz]
config.f_mod = 11.05E6;

% the range of time during which the data will be used for calculation ([tau_min, tau_max])[s]
config.tau = [0.4, 8]*1E-9;

% the data type. 'r' - ratio of x and y (x/y); 'a' - normalized amplitude (sqrt(x^2+y^2)); 'p' - phase (arctan(y/x))
config.mode = 'r';

% the mode to draw the theory curve
config.TheoryCurve_mode = 1;