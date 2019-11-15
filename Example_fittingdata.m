global config
% name of each layer
config.LayerName = {'Aluminum', 'PMMA', 'Silicon'};

% cross-plane thermal conductivity of each layer [W/mK]
config.kz = [240, 0.3, 140];

% in-plane thermal conductivity of each layer [W/mK]
config.kr = [240, 0.3, 140];

% isotropy or anisotropy for each layer. 0 - isotropy; 1 - anisotropy
% if a material is isotropy, please use the kz for fitting
config.iso = [0, 0, 0];

% Interfacial thermal conductance between each layer[W/m^2K]
config.G = [100, 100]*1E6;

% volume specific heat of each layer [J/m^3K]
config.vhc = [2.42, 1.4 1.62]*1E6;

% thickness of each layer [m]
config.d = [80, 30, 500E3]*1E-9;

% diameter of laser beam ([pump probe])[m]
config.w = [28.3 10.9]/2*1E-6;

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
config.fit_para = [
1, 3, 2.4E6, 2E6, 3E6
2, 1, 0.33, 1, 0.1
];  

% the maximum iteration of the fitting process
config.iteration = 3000;

% do fitting or don't do fitting. 0 - don't do fitting; 1 - do fitting
config.fitting_mode = 1;