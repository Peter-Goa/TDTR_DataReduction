function [] = TDTRDataProcess()
% Main function for the TDTR data processing.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% TDTRDataProcess
% Main function for the TDTR data processing.
% Reference:
% 1. main_TDTR_171120 in TDTR_Iwamoto_171120
% 2. The thesis of Aaron Jerome Schmidt
% Author: RL
% Date: Nov. 12, 2019
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

global config

clc,clear,clf
set(0,'defaultAxesFontSize',30)
set(0,'defaultTextFontSize',30)
set(0,'defaultAxesFontName','Helvetica')
set(0,'defaultTextFontName','Helvetica')

% the configuration file
[filename, path]=uigetfile('*.txt','Pick a configuration file');
if isequal(filename,0) || isequal(pathname,0)
    disp('User pressed cancel')
    return
end
config_file = fullfile(filename, path);
run(config_file);

if config.fitting_mode == 1
    if config.folder_mode == 1
        
    else
        
    end
    % the folder or file path of the data
    SourcePath = '';

    % the folder in which the result file will be
    OutputPath = '';
end
%