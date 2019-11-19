global config

kal=238;    val=2.42;
kpdds=15; vpdds=3.38;
kmgo=41;  vmgo=3.4;
% name of each layer
config.LayerName = {'Aluminum', 'PDDS', 'MgO'};

% cross-plane thermal conductivity of each layer [W/mK]
config.kz = [kal, kpdds, kmgo];

% in-plane thermal conductivity of each layer [W/mK]
config.kr = [kal, kpdds, kmgo];

% isotropy or anisotropy for each layer. 1 - isotropy; 0 - anisotropy
% if a material is isotropy, please use the kz for fitting
config.iso = [1, 1, 1];

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

% 0 - just process a single file whose path is SourcePath; 1 - process all of the files in the SourcePath
config.folder_mode = 1;

% parameters which will be fitting. It is a Npara matrix*5. 5 lines are Layer, Type, GaussValue, MinValue, MaxValue.
% Npara is the number of fitting parameter. Type can be one of the list{1 - 'kz', 2 - 'kr', 3 - 'vhc', 4 - 'd', 5 - 'G'}.
% If the material is isotropic, 'kz' rather than 'kr' should be used to get the thermal conductivity.
config.fit_para = [
1, 3, 2.4E6, 2E6, 3E6
1, 5, 100E6, 30E6, 500E6
2, 1, 15, 3, 25
2, 5, 100E6, 30E6, 500E6
2, 3, 3.4E6, 2.5E6, 4.5E6
3, 1, 45, 30, 60
3, 3, 3.5E6, 2E6, 5E6
];  

% the maximum iteration of the fitting process
config.iteration = 100;

% do fitting or don't do fitting. 0 - don't do fitting; 1 - do fitting
config.fitting_mode = 1;